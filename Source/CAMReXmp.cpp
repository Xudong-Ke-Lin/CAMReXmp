#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

/*-------------------------------------------------------------
 * DECLARING STATIC MEMBER DATA
 * -----------------------------------------------------------*/

int      CAMReXmp::verbose         = 0;
Real     CAMReXmp::cfl             = 0.9; // Default value - can be overwritten in settings file
int      CAMReXmp::do_reflux       = 1;  

int      CAMReXmp::NUM_STATE       = 16;  // Number of variables in the state
int      CAMReXmp::NUM_GROW        = 3;  // number of ghost cells

int CAMReXmp::NUM_STATE_FLUID = 10;
int CAMReXmp::NUM_STATE_MAXWELL = CAMReXmp::NUM_STATE - NUM_STATE_FLUID;

Vector<BCRec> CAMReXmp::bc;
Vector<BCRec> CAMReXmp::bc_EM;

int  CAMReXmp::max_order = 2;
int CAMReXmp::max_fmg_iter = 0;
Real CAMReXmp::soln_tol = 1e-10;
Real CAMReXmp::tau = 1.0;
MultiFab CAMReXmp::acoef;
std::array<MultiFab, BL_SPACEDIM> CAMReXmp::bcoeffs;

std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_lobc_X;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_hibc_X;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_lobc_Y;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_hibc_Y;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_lobc_Z;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_hibc_Z;

advanceFunctionPointer CAMReXmp::advanceWithChosenUpdateOrder = &CAMReXmp::StrangFirst;
RKFunctionPointer CAMReXmp::RKWithChosenUpdateOrder = &CAMReXmp::RK1;
sourceFunctionPointer CAMReXmp::sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateEX;
MaxwellFunctionPointer CAMReXmp::MaxwellSolverWithChosenMethod = &CAMReXmp::hyperbolicMaxwellRK;

int CAMReXmp::advanceOrder = 1;
int CAMReXmp::RKOrder = 1;
std::string CAMReXmp::MaxwellMethod = "HYP";
std::string CAMReXmp::sourceMethod = "EX";

//
//Default constructor.  Builds invalid object.
//
CAMReXmp::CAMReXmp ()
{
  // Flux registers store fluxes at patch boundaries to ensure fluxes are conservative between AMR levels
  flux_reg = 0;
}

//
//The basic constructor.
//
CAMReXmp::CAMReXmp (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
  :
  AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
  // Flux registers are only required if AMR is actually used, and if flux fix up is being done (recommended)
  flux_reg = 0;
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
  }
}

//
//The destructor.
//
CAMReXmp::~CAMReXmp () 
{
    delete flux_reg;
}

//
//Restart from a checkpoint file.
//
// AMReX can save simultion state such
// that if the code crashes, it can be restarted, with different
// settings files parameters if necessary (e.g. to output about the
// point of the crash).
//
void
CAMReXmp::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)
{
  AmrLevel::restart(papa,is,bReadSpecial);
  
  BL_ASSERT(flux_reg == 0);
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
  }
  
}

//
// Write a checkpoint file - format is handled automatically by AMReX
void 
CAMReXmp::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
}

//
//Write a plotfile to specified directory - format is handled automatically by AMReX.
//
void
CAMReXmp::writePlotFile (const std::string& dir,
	 	            std::ostream&      os,
                            VisMF::How         how)
{
  AmrLevel::writePlotFile (dir,os,how);
}
//
//Define data descriptors.
//
// This is how the variables in a simulation are defined.  In the case
// of the advection equation, a single variable, phi, is defined.
//
void
CAMReXmp::variableSetUp ()
{
  BL_ASSERT(desc_lst.size() == 0);

  // A function which contains all processing of the settings file,
  // setting up initial data, choice of numerical methods and
  // boundary conditions
  read_params();
  
  const int storedGhostZones = 0;
    
  // Setting up a container for a variable, or vector of variables:
  // Phi_Type: Enumerator for this variable type
  // IndexType::TheCellType(): AMReX can support cell-centred and vertex-centred variables (cell centred here)
  // StateDescriptor::Point: Data can be a point in time, or an interval over time (point here)
  // storedGhostZones: Ghost zones can be stored (e.g. for output).  Generally set to zero.
  // NUM_STATE: Number of variables in the variable vector (1 in the case of advection equation)
  // cell_cons_interp: Controls interpolation between levels - cons_interp is good for finite volume
  desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,storedGhostZones,NUM_STATE+2,
			 &cell_cons_interp);
  // from IAMR (https://github.com/AMReX-Codes/IAMR/blob/main/Source/NS_setup.cpp):
  // StateDescriptor::Interval and &node_bilinear_interp
  // from documentation (https://amrex-codes.github.io/amrex/doxygen/AMReX__Interpolater_8cpp.html)
  // &face_linear_interp or &face_divfree_interp
  // can use .ixType() to check the type, xface will give (N,C) -> Nodal in x and center in y
  IndexType xface(IntVect{AMREX_D_DECL(1,0,0)});
  IndexType yface(IntVect{AMREX_D_DECL(0,1,0)});
  desc_lst.addDescriptor(EM_X_Type,xface,
			 StateDescriptor::Point,storedGhostZones,6,
			 &cell_cons_interp);
  desc_lst.addDescriptor(EM_Y_Type,yface,
			 StateDescriptor::Point,storedGhostZones,6,
			 &cell_cons_interp);
  
  /*//Set up boundary conditions, all boundaries can be set
  //independently, including for individual variables, but lo (left) and hi (right) are useful ways to
  //store them, for consistent access notation for the boundary
  //locations
  int lo_bc[amrex::SpaceDim];
  int hi_bc[amrex::SpaceDim];
  // AMReX has pre-set BCs, including periodic (int_dir) and transmissive (foextrap)
  for (int i = 0; i < amrex::SpaceDim; ++i) {
    lo_bc[i] = hi_bc[i] = BCType::int_dir;   // periodic boundaries
  }


  // Object for storing all the boundary conditions
  BCRec bc(lo_bc, hi_bc);*/

  // added by 2020D

  bc.resize(NUM_STATE+2);
  bc_EM.resize(6);
  
  ParmParse pp;
  pp.query("test",test);
  
  if (test=="BrioWu")
    {
      for (int idim = 0; idim < amrex::SpaceDim; ++idim)
	{
	  for (int n = 0; n < NUM_STATE+2; ++n)
	    {
	      bc[n].setLo(idim, BCType::foextrap); // transmissive
	      bc[n].setHi(idim, BCType::foextrap);
	    }
	}

    } else if (test=="OT")
    {
      for (int i = 0; i < amrex::SpaceDim; ++i)
	{
	  for (int n = 0; n < NUM_STATE+2; ++n)
	    {
	      bc[n].setLo(i, BCType::int_dir); // interior 
	      bc[n].setHi(i, BCType::int_dir);
	    }
	}
    } else if (test=="Harris_sheet")
    {
      for (int n = 0; n < NUM_STATE+2; ++n)
	{
	  bc[n].setLo(0, BCType::int_dir);
	  bc[n].setHi(0, BCType::int_dir);
	}
      // conducting wall boundary
      bc[RHO_I].setLo(1, BCType::foextrap);      
      bc[RHO_E].setLo(1, BCType::foextrap);      

      bc[MOMX_I].setLo(1, BCType::foextrap);
      bc[MOMY_I].setLo(1, BCType::reflect_odd);
      bc[MOMZ_I].setLo(1, BCType::foextrap);
      bc[MOMX_E].setLo(1, BCType::foextrap);
      bc[MOMY_E].setLo(1, BCType::reflect_odd);
      bc[MOMZ_E].setLo(1, BCType::foextrap);

      bc[ENER_I].setLo(1, BCType::foextrap);
      bc[ENER_E].setLo(1, BCType::foextrap);
      
      bc[BX].setLo(1, BCType::foextrap);
      bc[BY].setLo(1, BCType::reflect_odd);
      bc[BZ].setLo(1, BCType::foextrap);

      bc[EX].setLo(1, BCType::reflect_odd);
      bc[EY].setLo(1, BCType::foextrap);
      bc[EZ].setLo(1, BCType::reflect_odd);

      bc[RHO_I].setHi(1, BCType::foextrap);
      bc[RHO_E].setHi(1, BCType::foextrap);

      bc[MOMX_I].setHi(1, BCType::foextrap);
      bc[MOMY_I].setHi(1, BCType::reflect_odd);
      bc[MOMZ_I].setHi(1, BCType::foextrap);
      bc[MOMX_E].setHi(1, BCType::foextrap);
      bc[MOMY_E].setHi(1, BCType::reflect_odd);
      bc[MOMZ_E].setHi(1, BCType::foextrap);

      bc[ENER_I].setHi(1, BCType::foextrap);
      bc[ENER_E].setHi(1, BCType::foextrap);
      
      bc[BX].setHi(1, BCType::foextrap);
      bc[BY].setHi(1, BCType::reflect_odd);
      bc[BZ].setHi(1, BCType::foextrap);

      bc[EX].setHi(1, BCType::reflect_odd);
      bc[EY].setHi(1, BCType::foextrap);
      bc[EZ].setHi(1, BCType::reflect_odd);

      bc[DIVB].setLo(1, BCType::foextrap);
      bc[DIVB].setHi(1, BCType::foextrap);
      bc[DIVE].setLo(1, BCType::foextrap);
      bc[DIVE].setHi(1, BCType::foextrap);
      
    } else if (test=="zpinch2d")
    {
      /*// top and bottom periodic boundaries
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(1, BCType::int_dir);
	  bc[n].setHi(1, BCType::int_dir);
	}
      
      // left use reflective boundaries
      // reflective even for most variables
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(0, BCType::reflect_even);
	  bc[n].setHi(0, BCType::foextrap);
	}
      // reflective odd for vectors in radial and azimuthal direction
      bc[MOMX_I].setLo(0, BCType::reflect_odd);
      bc[MOMY_I].setLo(0, BCType::reflect_odd);
      bc[MOMX_E].setLo(0, BCType::reflect_odd);	  
      bc[MOMY_E].setLo(0, BCType::reflect_odd);	  
      bc[BX].setLo(0, BCType::reflect_odd);
      bc[BY].setLo(0, BCType::reflect_odd);
      bc[EX].setLo(0, BCType::reflect_odd);
      bc[EY].setLo(0, BCType::reflect_odd);

      // right use conducting wall boundary
      bc[RHO_I].setHi(0, BCType::foextrap);
      bc[RHO_E].setHi(0, BCType::foextrap);

      bc[MOMX_I].setHi(0, BCType::reflect_odd);
      bc[MOMY_I].setHi(0, BCType::foextrap);
      bc[MOMZ_I].setHi(0, BCType::foextrap);
      bc[MOMX_E].setHi(0, BCType::reflect_odd);
      bc[MOMY_E].setHi(0, BCType::foextrap);
      bc[MOMZ_E].setHi(0, BCType::foextrap);

      bc[ENER_I].setHi(0, BCType::foextrap);
      bc[ENER_E].setHi(0, BCType::foextrap);
      
      bc[BX].setHi(0, BCType::reflect_odd);
      bc[BY].setHi(0, BCType::foextrap);
      bc[BZ].setHi(0, BCType::foextrap);

      bc[EX].setHi(0, BCType::foextrap);
      bc[EY].setHi(0, BCType::reflect_odd);
      bc[EZ].setHi(0, BCType::reflect_odd);
      */
      // top and bottom periodic boundaries
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(0, BCType::int_dir);
	  bc[n].setHi(0, BCType::int_dir);
	}
      
      // left use reflective boundaries
      // reflective even for most variables
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(1, BCType::reflect_even);
	  bc[n].setHi(1, BCType::foextrap);
	}
      // reflective odd for vectors in radial and azimuthal direction
      bc[MOMY_I].setLo(1, BCType::reflect_odd);
      bc[MOMZ_I].setLo(1, BCType::reflect_odd);
      bc[MOMY_E].setLo(1, BCType::reflect_odd);	  
      bc[MOMZ_E].setLo(1, BCType::reflect_odd);	  
      bc[BY].setLo(1, BCType::reflect_odd);
      bc[BZ].setLo(1, BCType::reflect_odd);
      bc[EY].setLo(1, BCType::reflect_odd);
      bc[EZ].setLo(1, BCType::reflect_odd);

      // right use conducting wall boundary
      bc[RHO_I].setHi(1, BCType::foextrap);
      bc[RHO_E].setHi(1, BCType::foextrap);

      bc[MOMY_I].setHi(1, BCType::reflect_odd);
      bc[MOMZ_I].setHi(1, BCType::foextrap);
      bc[MOMX_I].setHi(1, BCType::foextrap);
      bc[MOMY_E].setHi(1, BCType::reflect_odd);
      bc[MOMZ_E].setHi(1, BCType::foextrap);
      bc[MOMX_E].setHi(1, BCType::foextrap);

      bc[ENER_I].setHi(1, BCType::foextrap);
      bc[ENER_E].setHi(1, BCType::foextrap);
      
      bc[BY].setHi(1, BCType::reflect_odd);
      bc[BZ].setHi(1, BCType::foextrap);
      bc[BX].setHi(1, BCType::foextrap);

      bc[EY].setHi(1, BCType::foextrap);
      bc[EZ].setHi(1, BCType::reflect_odd);
      bc[EX].setHi(1, BCType::reflect_odd);
      
    } else if (test=="convergence")
    {
      for (int i = 0; i < amrex::SpaceDim; ++i)
	{
	  for (int n = 0; n < NUM_STATE; ++n)
	    {
	      bc[n].setLo(i, BCType::int_dir); // interior                                                                 
	      bc[n].setHi(i, BCType::int_dir);
	    }
	}
    } else
    {
      amrex::Abort("Please enter valid test in inputs file");
    }

  // bc array for the electromagnetic field
  for (int n=BX_LOCAL; n<=EZ_LOCAL; n++)
    {
      bc_EM[n] = bc[n+NUM_STATE_FLUID];
    }
  // Set up variable-specific information; needs to be done for each variable in NUM_STATE
  // Phi_Type: Enumerator for the variable type being set
  // 0: Position of the variable in the variable vector.  Single variable for advection.
  // phi: Name of the variable - appears in output to identify what is being plotted
  // bc: Boundary condition object for this variable (defined above)
  // BndryFunc: Function for setting boundary conditions.  For basic BCs, AMReX can handle these automatically
  desc_lst.setComponent(Phi_Type, RHO_I, "density_i", bc[RHO_I], 
			StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, MOMX_I, "momx_i", bc[MOMX_I],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, MOMY_I, "momy_i", bc[MOMY_I],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, MOMZ_I, "momz_i", bc[MOMZ_I],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, ENER_I, "energy_i", bc[ENER_I],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, RHO_E, "density_e", bc[RHO_E],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, MOMX_E, "momx_e", bc[MOMX_E],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, MOMY_E, "momy_e", bc[MOMY_E],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, MOMZ_E, "momz_e", bc[MOMZ_E],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, ENER_E, "energy_e", bc[ENER_E],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, BX, "magx", bc[BX],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, BY, "magy", bc[BY],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, BZ, "magz", bc[BZ],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, EX, "elecx", bc[EX],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, EY, "elecy", bc[EY],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, EZ, "elecz", bc[EZ],
			StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, DIVB, "divBerror", bc[DIVB],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(Phi_Type, DIVE, "divEerror", bc[DIVE],
			StateDescriptor::BndryFunc(nullfill));

  // face-cenctred primary variables in 2D
  desc_lst.setComponent(EM_X_Type, BX_LOCAL, "magxFace", bc[BX],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, BY_LOCAL, "magyFace", bc[BY],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, EX_LOCAL, "elecxFace", bc[EX],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, EY_LOCAL, "elecyFace", bc[EY],
                        StateDescriptor::BndryFunc(nullfill));
  // face-centred moments of the primary variables in 2D
  // used for the divergence-free (second or higher order) reconstruction
  desc_lst.setComponent(EM_X_Type, BY_LOCAL, "magyFaceX", bc[BY],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, BZ_LOCAL, "magzFaceX", bc[BZ],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, BX_LOCAL, "magxFaceY", bc[BX],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, BZ_LOCAL, "magzFaceY", bc[BZ],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, EY_LOCAL, "elecyFaceX", bc[EY],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, EZ_LOCAL, "eleczFaceX", bc[EZ],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, EX_LOCAL, "elecxFaceY", bc[EX],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, EZ_LOCAL, "eleczFaceY", bc[EZ],
                        StateDescriptor::BndryFunc(nullfill));
  
}

//
//Cleanup data descriptors at end of run.
//
void
CAMReXmp::variableCleanUp () 
{
    desc_lst.clear();
}


//
//Initialize data on this level from another CAMReXmp (during regrid).
// These are standard AMReX commands which are unlikely to need altering
//
void
CAMReXmp::init (AmrLevel &old)
{
  
  CAMReXmp* oldlev = (CAMReXmp*) &old;
  //
  // Create new grid data by fillpatching from old.
  //
  Real dt_new    = parent->dtLevel(level);
  Real cur_time  = oldlev->state[Phi_Type].curTime();
  Real prev_time = oldlev->state[Phi_Type].prevTime();
  Real dt_old    = cur_time - prev_time;
  setTimeLevel(cur_time,dt_old,dt_new);
  
  MultiFab& S_new = get_new_data(Phi_Type);
  
  const int zeroGhosts = 0;
  // FillPatch takes the data from the first argument (which contains
  // all patches at a refinement level) and fills (copies) the
  // appropriate data onto the patch specified by the second argument:
  // old: Source data
  // S_new: destination data
  // zeroGhosts: If this is non-zero, ghost zones could be filled too - not needed for init routines
  // cur_time: AMReX can attempt interpolation if a different time is specified - not recommended for advection eq.
  // Phi_Type: Specify the type of data being set
  // 0: This is the first data index that is to be copied
  // NUM_STATE: This is the number of states to be copied
  FillPatch(old, S_new, zeroGhosts, cur_time, Phi_Type, 0, NUM_STATE);
  
  // Note: In this example above, the all states in Phi_Type (which is
  // only 1 to start with) are being copied.  However, the FillPatch
  // command could be used to create a velocity vector from a
  // primitive variable vector.  In this case, the `0' argument is
  // replaced with the position of the first velocity component in the
  // primitive variable vector, and the NUM_STATE arguement with the
  // dimensionality - this argument is the number of variables that
  // are being filled/copied, and NOT the position of the final
  // component in e.g. the primitive variable vector.
}

//
//Initialize data on this level after regri_EM[0]dding if old level did not previously exist
// These are standard AMReX commands which are unlikely to need altering
//
void
CAMReXmp::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    
    // See first init function for documentation
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Advance grids at this level in time.
//  This function is the one that actually calls the flux functions.
//
Real
CAMReXmp::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{
  
  MultiFab& S_mm = get_new_data(Phi_Type);
  
  // This ensures that all data computed last time step is moved from
  // `new' data to `old data' - this should not need changing If more
  // than one type of data were declared in variableSetUp(), then the
  // loop ensures that all of it is updated appropriately
  for (int k = 0; k < NUM_STATE_TYPE; k++) {
    state[k].allocOldData();
    state[k].swapTimeLevels(dt);
  }

  // S_new is the MultiFab that will be operated upon to update the data
  MultiFab& S_new = get_new_data(Phi_Type);

  // Set up a multifab that will contain the electromagnetic fields
  MultiFab& S_EM_X_new = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_new = get_new_data(EM_Y_Type);
  
  const Real prev_time = state[Phi_Type].prevTime();
  const Real cur_time = state[Phi_Type].curTime();
  const Real ctr_time = 0.5*(prev_time + cur_time);

  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
 
  //
  // Get pointers to Flux registers, or set pointer to zero if not there.
  //
  FluxRegister *fine    = 0;
  FluxRegister *current = 0;
    
  int finest_level = parent->finestLevel();

  // If we are not on the finest level, fluxes may need correcting
  // from those from finer levels.  To start this process, we set the
  // flux register values to zero
  if (do_reflux && level < finest_level) {
    fine = &getFluxReg(level+1);
    fine->setVal(0.0);
  }

  // If we are not on the coarsest level, the fluxes are going to be
  // used to correct those on coarser levels.  We get the appropriate
  // flux level to include our fluxes within
  if (do_reflux && level > 0)
  {
    current = &getFluxReg(level);
  }

  // Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes[amrex::SpaceDim];
  
  // Define the appropriate size for the flux MultiFab.
  // Fluxes are defined at cell faces - this is taken care of by the
  // surroundingNodes(j) command, ensuring the size of the flux
  for (int j = 0; j < amrex::SpaceDim; j++)
  {
    BoxArray ba = S_new.boxArray();
    ba.surroundingNodes(j);
    fluxes[j].define(ba, dmap, NUM_STATE, 0);
  }

  // Set up a multifab that will contain the fluxes for the electromagnetic fields
  // It will be a nodal multifab (i.e. stored at the corner)
  MultiFab fluxesEM;
  BoxArray ba = S_EM_X_new.boxArray();
  ba.surroundingNodes(1);
  const DistributionMapping& dmX = S_EM_X_new.DistributionMap();
  fluxesEM.define(ba, dmX, 6, 0);
  fluxesEM = 0.0;
  
  // State with ghost cells - this is used to compute fluxes and perform the update.
  MultiFab Sborder(grids, dmap, NUM_STATE+2, NUM_GROW);
  // See init function for details about the FillPatch function
  FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE+2);
  // Fill periodic boundaries where they exist.  More accurately, the
  // FillBoundary call will fill overlapping boundaries (with periodic
  // domains effectively being overlapping).  It also takes care of
  // AMR patch and CPU boundaries.
  Sborder.FillBoundary(geom.periodicity());

  // Fill non-periodic physical boundaries
  FillDomainBoundary(Sborder, geom, bc);
  
  // Set up a multifab that will contain the electromagnetic fields
  //MultiFab S_EM[2];
  Array<MultiFab,2> S_EM;

  S_EM[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  S_EM[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM[0], NUM_GROW, time, EM_X_Type, 0, 6);
  FillPatch(*this, S_EM[1], NUM_GROW, time, EM_Y_Type, 0, 6);
  
  S_EM[0].FillBoundary(geom.periodicity());
  S_EM[1].FillBoundary(geom.periodicity());

  //std::cout << Sborder.ixType().cellCentered() << " " << S_EM[0].ixType().cellCentered() << " " << S_EM[1].ixType().cellCentered() << std::endl;
  
  // added by 2020D 
  // Fill non-periodic physical boundaries 
  FillDomainBoundary(S_EM[0], geom, bc_EM);
  FillDomainBoundary(S_EM[1], geom, bc_EM);  

  // We use FillPatcher to do fillpatch here if we can
  //FillPatcherFill(SY, 0, 2, NUM_GROW, time, 3, 0);
  //FillPatch(*this, SY, NUM_GROW, time, 3, 0, 2);
  //SY.FillBoundary(geom.periodicity());
  //FillDomainBoundary(SY, geom, bc3);
  
  
  //FillPatcherFill(S_EM[1], 0, 6, NUM_GROW, time, EM_Y_Type, 0);
  //FillPatcherFill(Sborder, 0, NUM_STATE, NUM_GROW, time, Phi_Type, 0);
  //MFIter::allowMultipleMFIters(true);

  MultiFab S0(grids, dmap, NUM_STATE+2, NUM_GROW);
  // See init function for details about the FillPatch function
  FillPatch(*this, S0, NUM_GROW, time, Phi_Type, 0, NUM_STATE+2);
  S0.FillBoundary(geom.periodicity());
  // Fill non-periodic physical boundaries
  FillDomainBoundary(S0, geom, bc);


  //Array<MultiFab,2> Slopes;
  //for (int n = 0; n<nSlopes; n++)
  //Slopes[n].define(grids, dmap, NUM_STATE, NUM_GROW);
  MultiFab Slopes(grids, dmap, NUM_STATE_FLUID*nSlopes, NUM_GROW);
  
  //(this->*advanceWithChosenUpdateOrder)(Sborder,fluxes,dx,dt);   

  const int iOffset = 1;
  const int jOffset = ( amrex::SpaceDim > 1 ? 1 : 0);
  const int kOffset = ( amrex::SpaceDim == 3 ? 1 : 0);

  // Compute slopes for WENO reconstruction
  for (MFIter mfi(Slopes, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& slopes = Slopes.array(mfi);

      for(int k = lo.z-kOffset; k <= hi.z+kOffset; k++)
	{
	  for(int j = lo.y-jOffset; j <= hi.y+jOffset; j++)
	    {
	      for(int i = lo.x-iOffset; i <= hi.x+iOffset; i++)
		{		  
		  std::array<Vector<Real>, 2> slopesX = WENO_slopes(arr, i, j, k, 1, 0, 0);
		  for (int n = 0; n<NUM_STATE_FLUID; n++)
		    {
		      slopes(i,j,k,n) = slopesX[0][n];		      
		      slopes(i,j,k,n+NUM_STATE_FLUID) = slopesX[1][n];
		    }
#if (AMREX_SPACEDIM >= 2)
		  std::array<Vector<Real>, 2> slopesY = WENO_slopes(arr, i, j, k, 0, 1, 0); 
		  for (int n = 0; n<NUM_STATE_FLUID; n++)
		    {
		      slopes(i,j,k,n+2*NUM_STATE_FLUID) = slopesY[0][n];
		      slopes(i,j,k,n+3*NUM_STATE_FLUID) = slopesY[1][n];
		    }
		  Vector<Real> slopesCross = WENO_slopesCross(arr, slopes, i, j, k, 1, 1, 0);
		  for (int n = 0; n<NUM_STATE_FLUID; n++)
		    {
		      slopes(i,j,k,n+4*NUM_STATE_FLUID) = slopesCross[n]; 
		    }
#endif
		  
		}
	    }
	}
    }
  
  // Update cell-centred fluid variables and z-components of EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
  {

    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);

    // Loop over all the patches at this level
    for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& fluxArr = fluxes[d].array(mfi);
      const auto& arrEM = S_EM[d].array(mfi);      


      const auto& slopes = Slopes.array(mfi);
      //const auto& slopesxx = Slopes[1].array(mfi);
#if (AMREX_SPACEDIM == 1)      
      //const std::array<Array4<Real>, nSlopes> slopes = {slopesx,slopesxx};
#endif
#if (AMREX_SPACEDIM == 2)
      //const auto& slopesy = Slopes[2].array(mfi);
      //const auto& slopesyy = Slopes[3].array(mfi);
      //const auto& slopesxy = Slopes[4].array(mfi);
      //const std::array<Array4<Real>, nSlopes> slopes = {slopesx,slopesxx,slopesy,slopesyy,slopesxy};
#endif
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = MUSCL_Hancock_WENOHLLC_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
							      0, NUM_STATE_FLUID/2, dx[d], dt, d,
							      fluidFlux, HLLC);
	    Vector<Real> flux_e = MUSCL_Hancock_WENOHLLC_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
							      NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, dx[d], dt, d,
							      fluidFlux, HLLC);	    
	    //Vector<Real> flux = fluid_flux_HLLC(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);

	    for(int n=0; n<NUM_STATE_FLUID/2; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
		fluxArr(i,j,k,n+NUM_STATE_FLUID/2) = flux_e[n];
	      }
	    
	    // Fluxes for the z-components of the EM fields because it is 2D code
	    //Vector<Real> fluxEM = Maxwell_flux_Godunov(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);
	    //Vector<Real> fluxEM = MUSCL_Hancock_WENOGodunov_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);
	    /*Vector<Real> fluxEM = MUSCL_Hancock_WENOHLLC_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
	    						      NUM_STATE_FLUID, NUM_STATE_MAXWELL, dx[d], dt, d,
							      MaxwellFlux, Godunov);
	    for(int n=NUM_STATE_FLUID; n<NUM_STATE; n++)
	      {		
		fluxArr(i,j,k,n) = fluxEM[n-NUM_STATE_FLUID];		
	      }
	    */
	    /*
	    // Calculate transverse facial components from cell-centred components
	    // transverse values, e.g. By, Bz, Ey and Ez for the x-face
	    Vector<Real> transverseEM = Maxwell_transverse_comp(arr, i, j, k, iOffset, jOffset, kOffset);	    
	    arrEM(i,j,k,BX_LOCAL+(1+d)%3) = transverseEM[BX_LOCAL+(1+d)%3];
	    arrEM(i,j,k,BX_LOCAL+(2+d)%3) = transverseEM[BX_LOCAL+(2+d)%3];
	    arrEM(i,j,k,EX_LOCAL+(1+d)%3) = transverseEM[EX_LOCAL+(1+d)%3];
	    arrEM(i,j,k,EX_LOCAL+(2+d)%3) = transverseEM[EX_LOCAL+(2+d)%3];	   
	    */
	    /*
	    std::array<Vector<Real>, 2> flux_and_transverse = MUSCL_Hancock_Godunov_Maxwell_and_transverse(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);
	    fluxArr(i,j,k,BZ) = flux_and_transverse[0][BZ_LOCAL];
	    fluxArr(i,j,k,EZ) = flux_and_transverse[0][EZ_LOCAL];

	    arrEM(i,j,k,BX_LOCAL+(1+d)%3) = flux_and_transverse[1][BX_LOCAL+(1+d)%3];
	    arrEM(i,j,k,BX_LOCAL+(2+d)%3) = flux_and_transverse[1][BX_LOCAL+(2+d)%3];
	    arrEM(i,j,k,EX_LOCAL+(1+d)%3) = flux_and_transverse[1][EX_LOCAL+(1+d)%3];
	    arrEM(i,j,k,EX_LOCAL+(2+d)%3) = flux_and_transverse[1][EX_LOCAL+(2+d)%3];
	    */	
	  }
	}
      }
      
      for(int k = lo.z; k <= hi.z; k++)
      {
      	for(int j = lo.y; j <= hi.y; j++)
      	{
      	  for(int i = lo.x; i <= hi.x; i++)
      	  {
      	    // Update fluid variables
      	    //for(int n=0; n<NUM_STATE; n++)
      	    for(int n=0; n<NUM_STATE_FLUID; n++)
              {
      		// Conservative update formula
      		arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
              }
      	    /*
      	    // Update cell-centred z-components becuause it is 2D code
      	    arr(i,j,k,BZ) = arr(i,j,k,BZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, BZ) - fluxArr(i,j,k,BZ));
      	    arr(i,j,k,EZ) = arr(i,j,k,EZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, EZ) - fluxArr(i,j,k,EZ));
      	    */
      	    // Initialise to zero
      	    arr(i,j,k,DIVB) = 0.0;
      	    arr(i,j,k,DIVE) = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
      	  }
      	}
      }
    }
      
    // We need to compute boundary conditions again after each update
    Sborder.FillBoundary(geom.periodicity());
    //S_EM[0].FillBoundary(geom.periodicity());
    //S_EM[1].FillBoundary(geom.periodicity());
     
    // added by 2020D 
    // Fill non-periodic physical boundaries
    FillDomainBoundary(Sborder, geom, bc);
    //FillDomainBoundary(S_EM[0], geom, bc_EM);
    //FillDomainBoundary(S_EM[1], geom, bc_EM);
    
    // The fluxes now need scaling for the reflux command.
    // This scaling is by the size of the boundary through which the flux passes, e.g. the x-flux needs scaling by the dy, dz and dt
    /*if(do_reflux)
    {
      Real scaleFactor = dt;
      for(int scaledir = 0; scaledir < amrex::SpaceDim; ++scaledir)
      {
	// Fluxes don't need scaling by dx[d]
	if(scaledir == d)
	{
	  continue;
	}
	scaleFactor *= dx[scaledir];
      }
      // The mult function automatically multiplies entries in a multifab by a scalar
      // scaleFactor: The scalar to multiply by
      // 0: The first data index in the multifab to multiply
      // NUM_STATE:  The total number of data indices that will be multiplied
      fluxes[d].mult(scaleFactor, 0, NUM_STATE);
      }*/
  }
  
  /*
  // Unsplit hyperbolic update  
  for (int d = 0; d < amrex::SpaceDim ; d++)   
  {

    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);

    // Loop over all the patches at this level
    for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& fluxArr = fluxes[d].array(mfi);
      const auto& arrEM = S_EM[d].array(mfi);
      //const auto& fluxArrEM = fluxesEM.array(mfi);

      //amrex::Print() << fluxArr(10,5,0,RHO_I) << " " << fluxArr(10,5,0,RHO_E) << std::endl;

      
      //std::cout << " " << fluxArr << " " << arr << " " << fluxArrEM << " " << arrEM << " " << lo.x << " " << hi.x << std::endl;      
      //amrex::Abort();
      
  
      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    for(int n=0; n<NUM_STATE_FLUID; n++)
              {
		// Conservative update formula
		arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
              }
	    // Update cell-centred z-components becuause it is 2D code
	    //arr(i,j,k,BZ) = arr(i,j,k,BZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, BZ) - fluxArr(i,j,k,BZ));
	    //arr(i,j,k,EZ) = arr(i,j,k,EZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, EZ) - fluxArr(i,j,k,EZ));
	    
	    // Initialise to zero
	    arr(i,j,k,DIVB) = 0.0;
	    arr(i,j,k,DIVE) = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
	  }
	}
      }    
    }
    
    // We need to compute boundary conditions again after each update
    Sborder.FillBoundary(geom.periodicity());
     
    // added by 2020D 
    // Fill non-periodic physical boundaries                      
    FillDomainBoundary(Sborder, geom, bc);
  }  
  */  


  // At this point the current at the faces are updated in fluxes
  // Sborder has updated fluid variables and has cell-centred EM variables
  // S_EM has primary facial variables already defined at t=t_old and reconstructed transverse components
  // divergence free update including the source terms for the electric field
  //hyperbolicMaxwellSolverDivFree(S_EM,fluxesEM,Sborder,fluxes,dx,dt);
  
  /*// Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes0[amrex::SpaceDim];
  for (int j = 0; j < amrex::SpaceDim; j++)
  {
    BoxArray ba = S_new.boxArray();
    ba.surroundingNodes(j);
    fluxes0[j].define(ba, dmap, NUM_STATE, 0);
  }
  
  // State with ghost cells - this is used to compute fluxes and perform the update.
  MultiFab S0(grids, dmap, NUM_STATE+2, NUM_GROW);
  // See init function for details about the FillPatch function
  FillPatch(*this, S0, NUM_GROW, time, Phi_Type, 0, NUM_STATE+2);
  S0.FillBoundary(geom.periodicity());
  // Fill non-periodic physical boundaries
  FillDomainBoundary(S0, geom, bc);
  hyperbolicMaxwellSolverDivFreeIneff(S_EM,fluxesEM,Sborder,fluxes,S0,fluxes0,dx,dt);
  */

  // // Compute cell-centred EM fields from face-centred
  // for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
  //   {
  //     const Box& bx = mfi.tilebox();
      
  //     const Dim3 lo = lbound(bx);
  //     const Dim3 hi = ubound(bx);

  //     const Dim3 hiDomain = ubound(geom.Domain());

  //     Array4<Real> arr = Sborder.array(mfi);
  //     Array4<Real> arrEMX = S_EM[0].array(mfi);
  //     Array4<Real> arrEMY = S_EM[1].array(mfi);

  //     for(int k = lo.z; k <= hi.z; k++)
  //     {
  //       for(int j = lo.y; j <= hi.y; j++)
  //       {
  //         for(int i = lo.x; i <= hi.x; i++)
  //         {
  // 	    arr(i,j,k,BX) = 0.5*(arrEMX(i,j,k,BX_LOCAL)+arrEMX(i+1,j,k,BX_LOCAL));
  // 	    arr(i,j,k,BY) = 0.5*(arrEMY(i,j,k,BY_LOCAL)+arrEMY(i,j+1,k,BY_LOCAL));
  // 	    arr(i,j,k,EX) = 0.5*(arrEMX(i,j,k,EX_LOCAL)+arrEMX(i+1,j,k,EX_LOCAL));
  // 	    arr(i,j,k,EY) = 0.5*(arrEMY(i,j,k,EY_LOCAL)+arrEMY(i,j+1,k,EY_LOCAL));

  // 	    Real ax = (arrEMX(i+1,j,k,EX_LOCAL)-arrEMX(i,j,k,EX_LOCAL))/dx[0];
  // 	    Real by = (arrEMY(i,j+1,k,EY_LOCAL)-arrEMY(i,j,k,EY_LOCAL))/dx[1];
  // 	    Real cz = 0.0;

  // 	    //if (std::abs(ax+by+cz-1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))) > 1e-13)
  // 	    //std::cout << i << " " << j << " " << ax+by+cz << " " << 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E)) << std::endl;
  // 	    /*
  // 	    if (i==1 && bc[MOMX_I].lo(0) == BCType::reflect_odd)
  // 	      arr(0,j,k,EX_LOCAL) = arr(1,j,k,EX_LOCAL);
  // 	    if (i==hiDomain.x && bc[MOMX_I].hi(0) == BCType::reflect_odd)
  // 	      arr(hiDomain.x,j,k,EX_LOCAL) = arr(hiDomain.x-1,j,k,EX_LOCAL);
  // 	    if (j==1 && bc[MOMY_I].lo(1) == BCType::reflect_odd)
  // 	      arr(i,0,k,EY_LOCAL) = arr(i,1,k,EY_LOCAL);
  // 	    if (j==hiDomain.y && bc[MOMY_I].hi(1) == BCType::reflect_odd)
  // 	      arr(i,hiDomain.y,k,EY_LOCAL) = arr(i,hiDomain.y-1,k,EY_LOCAL);
  // 	    */	    
  // 	    // Ez source terms
  // 	    //arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
  // 	  }
  // 	}
  //     }      
  //   }
  // // We need to compute boundary conditions again after each update       
  // Sborder.FillBoundary(geom.periodicity());
	    
  // // Fill non-periodic physical boundaries          
  // FillDomainBoundary(Sborder, geom, bc);
  
  // At this point the current at the faces are updated in fluxes
  // Sborder has updated fluid variables and has cell-centred EM variables
  // S_EM has primary facial variables already defined at t=t_old and reconstructed transverse components
  // S0 has non-updated fluid variables used for charge density reconstruction
  // divergence free update and reconstruction including the source terms for the electric field
  // also includes computation of cell-centred EM fields from face-centred by using coeff. that give higher order
  //hyperbolicMaxwellSolverDivFreeBalsara(S_EM,fluxesEM,Sborder,fluxes,S0,dx,dt);
  MaxwellSolverDivFreeWENO(S_EM,fluxesEM,Sborder,fluxes,S0,dx,dt);

  // Source term update
  // Loop over all the patches at this level         
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      Array4<Real> arr = Sborder.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
      {
        for(int j = lo.y; j <= hi.y; j++)
        {
          for(int i = lo.x; i <= hi.x; i++)
          {
	    sourceUpdateANEX(arr, i, j, k, dt);
	  }
	}
      }      
    }
  // We need to compute boundary conditions again after each update       
  Sborder.FillBoundary(geom.periodicity());
	    
  // Fill non-periodic physical boundaries
  FillDomainBoundary(Sborder, geom, bc);
  /*
  // Update source terms for face-centred EM fields (ANEX)
  for (int d = 0; d < amrex::SpaceDim ; d++)   
  {

    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);
    
    // Loop over all the patches at this level
    for (MFIter mfi(S_EM[d], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      //const auto& arr = Sborder.array(mfi);
      const auto& arrEM = S_EM[d].array(mfi);
      //const auto& fluxArrEM = fluxesEM.array(mfi);
      const auto& fluxArr = fluxes[d].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    // source terms
	    Real currentFace = r_i*fluxArr(i,j,k,RHO_I) + r_e*fluxArr(i,j,k,RHO_E);
	    arrEM(i,j,k,EX_LOCAL+d) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;
	    //Vector<Real> transverseEM = Maxwell_transverse_comp(arr, i, j, k, iOffset, jOffset, kOffset);
	    //arrEM(i,j,k,EX_LOCAL+d) = transverseEM[EX_LOCAL+d];	    	    
	    
	  }
	}
      }      
    }         
  }
  
  // Compute cell-centred EM fields from face-centred
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      Array4<Real> arr = Sborder.array(mfi);
      Array4<Real> arrEMX = S_EM[0].array(mfi);
      Array4<Real> arrEMY = S_EM[1].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
      {
        for(int j = lo.y; j <= hi.y; j++)
        {
          for(int i = lo.x; i <= hi.x; i++)
          {
	    arr(i,j,k,BX) = 0.5*(arrEMX(i,j,k,BX_LOCAL)+arrEMX(i+1,j,k,BX_LOCAL));
	    arr(i,j,k,BY) = 0.5*(arrEMY(i,j,k,BY_LOCAL)+arrEMY(i,j+1,k,BY_LOCAL));
	    arr(i,j,k,EX) = 0.5*(arrEMX(i,j,k,EX_LOCAL)+arrEMX(i+1,j,k,EX_LOCAL));
	    arr(i,j,k,EY) = 0.5*(arrEMY(i,j,k,EY_LOCAL)+arrEMY(i,j+1,k,EY_LOCAL));
	  }
	}
      }      
    }
  */
  // Compute divergence errors
  for (int d = 0; d < amrex::SpaceDim ; d++)
    {
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
      {
	const Box& bx = mfi.tilebox();
	  
	const Dim3 lo = lbound(bx);
	const Dim3 hi = ubound(bx);
	
	Array4<Real> arr = Sborder.array(mfi);
	Array4<Real> arrEM = S_EM[d].array(mfi);	
	
	//for(int k = lo.z; k <= hi.z+kOffset; k++)
	for(int k = lo.z; k <= hi.z; k++)
	{
	  //for(int j = lo.y; j <= hi.y+jOffset; j++)
	  for(int j = lo.y; j <= hi.y; j++)
	  {
	    //for(int i = lo.x; i <= hi.x+iOffset; i++)
	    for(int i = lo.x; i <= hi.x; i++)
	      {
		arr(i,j,k,DIVB) += (arrEM(i+iOffset,j+jOffset,k,BX_LOCAL+d)-arrEM(i,j,k,BX_LOCAL+d))/dx[d];
		// substract because divEerror was set to be the charge density
		arr(i,j,k,DIVE) -= (arrEM(i+iOffset,j+jOffset,k,EX_LOCAL+d)-arrEM(i,j,k,EX_LOCAL+d))/dx[d];
	      }
	  }
	}      
      }
    }

  // We need to compute boundary conditions again after each update       
  Sborder.FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries 
  FillDomainBoundary(Sborder, geom, bc);
  
  /*for (int d = 0; d < amrex::SpaceDim ; d++)   
  {

    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);

    int d_EM = (d==0) ? 1 : 0;
    
    // Loop over all the patches at this level
    for (MFIter mfi(S_EM[d_EM], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arrEM = S_EM[d_EM].array(mfi);
      const auto& fluxArrEM = fluxesEM.array(mfi);
      const auto& fluxArr = fluxes[d_EM].array(mfi);

      //std::cout << fluxArr << " " << fluxArrEM << " " << arrEM << " " << lo.x << " " << hi.x << " " << lo.y << " " << hi.y << std::endl;      
      // amrex::Abort();
      
      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {

	    // source terms
	    Vector<Real> currentFace = calculateCurrentFace(fluxArr, i, j, k);	  
	    arrEM(i,j,k,EX_LOCAL+(1+d)%3) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace[(1+d)%3];
	    arrEM(i,j,k,EX_LOCAL+(2+d)%3) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace[(2+d)%3];   
	  }
	}
      }
      
    }
      
    // We need to compute boundary conditions again after each update
    S_EM[0].FillBoundary(geom.periodicity());
    S_EM[1].FillBoundary(geom.periodicity());
     
    // added by 2020D 
    // Fill non-periodic physical boundaries                          
    FillDomainBoundary(S_EM[0], geom, bc_EM);
    FillDomainBoundary(S_EM[1], geom, bc_EM);
    
  }
  */

  // The updated data is now copied to the S_new multifab.  This means
  // it is now accessible through the get_new_data command, and AMReX
  // can automatically interpolate or extrapolate between layers etc.
  // S_new: Destination
  // Sborder: Source
  // Third entry: Starting variable in the source array to be copied (the zeroth variable in this case)
  // Fourth entry: Starting variable in the destination array to receive the copy (again zeroth here)
  // NUM_STATE: Total number of variables being copied
  // Sixth entry: Number of ghost cells to be included in the copy (zero in this case, since only real
  //              data is needed for S_new)
  MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE+2, 0);
  MultiFab::Copy(S_EM_X_new, S_EM[0], 0, 0, 6, 0);
  MultiFab::Copy(S_EM_Y_new, S_EM[1], 0, 0, 6, 0); 
  // Refluxing at patch boundaries.  Amrex automatically does this
  // where needed, but you need to state a few things to make sure it
  // happens correctly:
  // FineAdd: If we are not on the coarsest level, the fluxes at this level will form part of the correction
  //          to a coarse level
  // CrseInit:  If we are not the finest level, the fluxes at patch boundaries need correcting.  Since we
  //            know that the coarse level happens first, we initialise the boundary fluxes through this
  //            function, and subsequently FineAdd will modify things ready for the correction
  // Both functions have the same arguments:
  // First: Name of the flux MultiFab (this is done dimension-by-dimension
  // Second: Direction, to ensure the correct vertices are being corrected
  // Third: Source component - the first entry of the flux MultiFab that is to be copied (it is possible that
  //        some variables will not need refluxing, or will be computed elsewhere (not in this example though)
  // Fourth: Destinatinon component - the first entry of the flux register that this call to FineAdd sends to
  // Fifth: NUM_STATE - number of states being added to the flux register
  // Sixth: Multiplier - in general, the least accurate (coarsest) flux is subtracted (-1) and the most
  //        accurate (finest) flux is added (+1)
  if (do_reflux) {
    if (current) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
	current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
    }
    if (fine) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
	fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
    }
  }

  return dt;
}

//
//Estimate time step.
// This function is called by all of the other time step functions in AMReX, and is the only one that should
// need modifying
//
Real
CAMReXmp::estTimeStep (Real)
{
  // This is just a dummy value to start with 
  Real dt_est  = 1.0e+20;
  
  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
  const Real cur_time = state[Phi_Type].curTime();
  const MultiFab& S_new = get_new_data(Phi_Type);
  
  // This should not really be hard coded
  
  // State with ghost cells - this is used to compute fluxes and perform the update.
  MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
  // See init function for details about the FillPatch function
  FillPatch(*this, Sborder, NUM_GROW, state[Phi_Type].curTime(), Phi_Type, 0, NUM_STATE);
  Sborder.FillBoundary(geom.periodicity());
  // added by 2020D            
  // Fill non-periodic physical boundaries     
  FillDomainBoundary(Sborder, geom, bc);

  /*Vector<Real> dim;
  if (geom.Coord()==0)
    dim = {0,1};
  else if (geom.Coord()==1)
  dim = {0,2};*/
  /*ParmParse pp;
  std::string coor;
  pp.query("coor",coor);  
  if (coor=="xz")
    dim = {0,2};
  else
    dim = {0,1};*/

  MFIter::allowMultipleMFIters(true);

  for(unsigned int d = 0; d < amrex::SpaceDim; ++d)
  {
    Vector<Real> c_array {0.0};
    Vector<Real> omega_pe_array {0.0};
    Vector<Real> omega_ce_array {0.0};
    
    for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = Sborder.array(mfi);
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  // compute fastest speeds
		  Real v_x_e = arr(i,j,k,MOMX_E+d)/arr(i,j,k,RHO_E);		  
		  Real c_e = get_speed({arr(i,j,k,RHO_E),arr(i,j,k,MOMX_E),arr(i,j,k,MOMY_E),arr(i,j,k,MOMZ_E),arr(i,j,k,ENER_E)});
		  if (MaxwellMethod=="IM")
		    c_array.push_back(std::abs(v_x_e)+c_e);
		  else
		    c_array.push_back(v_x_e+c_e);
		  // compute electron frequencies
		  omega_pe_array.push_back(std::sqrt(m*m*arr(i,j,k,RHO_E)));
		  Real B = get_magnitude(arr(i,j,k,BX),arr(i,j,k,BY),arr(i,j,k,BZ));		  
		  omega_ce_array.push_back(m*B);
		}
	    }
	}
    }
    c_h = *std::max_element(c_array.begin(), c_array.end());
    Real omega_pe = *std::max_element(omega_pe_array.begin(), omega_pe_array.end());
    Real omega_ce = *std::max_element(omega_ce_array.begin(), omega_ce_array.end());

    // for implicit Maxwell solver use maximum fluid velocity
    if (MaxwellMethod=="IM")
      dt_est = std::min(dt_est, dx[d]/c_h);
    // for hyperbolic Maxwell solver use speed of light
    else
      dt_est = std::min(dt_est, dx[d]/std::max(c_h, c));
      //dt_est = std::min(dt_est, dx[d]/c_h);
      // 0.5 means subcycling two times
      //dt_est = std::min(dt_est, dx[d]/std::max(c_h, c/2.0));

    // dt also needs to resolve plasma and cyclotron frequencies
    dt_est = std::min(dt_est, 0.5*std::min(omega_pe,omega_ce));

  }
  
  
  // Ensure that we really do have the minimum across all processors
  ParallelDescriptor::ReduceRealMin(dt_est);
  ParallelDescriptor::ReduceRealMax(c_h);

  dt_est *= cfl;
  /*
  if (c_h>c)
    amrex::Abort("Fluid velocity is higher than speed of light!");
  
  if (dt_est<1e-8)
    amrex::Abort("Too small time step!");
  */
  if (verbose) {
    amrex::Print() << "CAMReXmp::estTimeStep at level " << level 
		   << ":  dt_est = " << dt_est << std::endl;
  }
  
  return dt_est;
}
//Compute initial time step.
//
Real
CAMReXmp::initialTimeStep ()
{
  return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
CAMReXmp::computeInitialDt (int                   finest_level,
	  	               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
  //
  // Grids have been constructed, compute dt for all levels.
  //
  // AMReX's AMR Level mode assumes that the time step only needs
  // calculating on the coarsest level - all subsequent time steps are
  // reduced by the refinement factor
  if (level > 0)
    return;

  // Initial guess
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor   *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_level[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[Phi_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

//
//Compute new `dt'.
//
void
CAMReXmp::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
  //
  // We are at the end of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  //
  if (level > 0)
    return;

  // Although we only compute the time step on the finest level, we
  // need to take information from all levels into account.  The
  // sharpest features may be smeared out on coarse levels, so not
  // using finer levels could cause instability
  for (int i = 0; i <= finest_level; i++)
  {
    CAMReXmp& adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  // A couple of things are implemented to ensure that time step's
  // don't suddenly grow by a lot, as this could lead to errors - for
  // sensible mesh refinement choices, these shouldn't really change
  // anything
  if (post_regrid_flag == 1) 
  {
    //
    // Limit dt's by pre-regrid dt
    //
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i],dt_level[i]);
    }
  }
  else 
  {
    //
    // Limit dt's by change_max * old dt
    //
    static Real change_max = 1.1;
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
    }
  }
    
  //
  // Find the minimum over all levels
  //
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_min[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[Phi_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }
  
  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

//
//Do work after timestep().
// If something has to wait until all processors have done their advance function, the post_timestep function
// is the place to put it.  Refluxing and averaging down are two standard examples for AMR
//
void
CAMReXmp::post_timestep (int iteration)
{
  //
  // Integration cycle on fine level grids is complete
  // do post_timestep stuff here.
  //
  int finest_level = parent->finestLevel();
  
  if (do_reflux && level < finest_level)
    reflux();
  
  if (level < finest_level)
    avgDown();
  
}

//
//Do work after regrid().
// Nothing normally needs doing here, but if something was calculated on a per-patch basis, new patches might
// this to be calcuated immediately
//
void
CAMReXmp::post_regrid (int lbase, int new_finest)
{

}

//
//Do work after a restart().
// Similar to post_regrid, nothing normally needs doing here
//
void
CAMReXmp::post_restart() 
{

}

//
//Do work after init().
// Once new patches have been initialised, work may need to be done to ensure consistency, for example,
// averaging down - though for linear interpolation, this probably won't change anything
//
void
CAMReXmp::post_init (Real stop_time)
{
  if (level > 0)
    return;
  //
  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  //
  int finest_level = parent->finestLevel();
  for (int k = finest_level-1; k>= 0; k--)
    getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//  Determine which parts of the domain need refinement
//
void
CAMReXmp::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
  const Real* dx        = geom.CellSize();
  const Real* prob_lo   = geom.ProbLo();

  MultiFab& S_new = get_new_data(Phi_Type);
  
  Vector<int> itags;
	
  for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
  {
    const Box&  tilebx  = mfi.tilebox();

    // An AMReX construction, effectively a boolean array which is true in positions that are valid for refinement
    TagBox&     tagfab  = tags[mfi];

    // Traditionally, a lot of the array-based operations in AMReX happened in Fortran.  The standard template
    // for these is short and easy to read, flagging on values or gradients (first order calculation)
    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
    // So we are going to get a temporary integer array.
    tagfab.get_itags(itags, tilebx);
	    
    // data pointer and index space
    int*        tptr    = itags.dataPtr();
    const int*  tlo     = tilebx.loVect();
    const int*  thi     = tilebx.hiVect();

    // Various macros exist to convert the C++ data structures to Fortran
    state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
		BL_TO_FORTRAN_3D(S_new[mfi]),
		&tagval, &clearval, 
		AMREX_ARLIM_3D(tilebx.loVect()), AMREX_ARLIM_3D(tilebx.hiVect()), 
		AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &level);
    //
    // Now update the tags in the TagBox.
    //
    tagfab.tags_and_untags(itags, tilebx);
  }
}

//
// This function reads the settings file
//
void
CAMReXmp::read_params ()
{
  // Make sure that this is only done once
  static bool done = false;

  if (done) return;

  done = true;

  // A ParmParse object allows settings, with the correct prefix, to be read in from the settings file
  // The prefix can help identify what a settings parameter is used for
  // AMReX has some default ParmParse names, amr and geometry are two commonly needed ones
  ParmParse pp("adv");   

  // ParmParse has two options; query and get.  Query will only alter
  // a parameter if it can be found (if these aren't in the settings
  // file, then the values at the top of this file will be used).  Get
  // will throw an error if the parameter is not found in the settings
  // file.
  pp.query("v",verbose);
  pp.query("cfl",cfl);
  pp.query("do_reflux",do_reflux);

  // Vector variables can be read in; these require e.g.\ pp.queryarr
  // and pp.getarr, so that the ParmParse object knows to look for
  // more than one variable

  // Geometries can be Cartesian, cylindrical or spherical - some
  // functions (e.g. divergence in linear solvers) are coded with this
  // geometric dependency
  Geometry const* gg = AMReX::top()->getDefaultGeometry();

  // This tutorial code only supports Cartesian coordinates.
  /*if (! gg->IsCartesian()) {
    amrex::Abort("Please set geom.coord_sys = 0");
    }*/

  /*// This tutorial code only supports periodic boundaries.
  // The periodicity is read from the settings file in AMReX source code, but can be accessed here
  if (! gg->isAllPeriodic()) {
    amrex::Abort("Please set geometry.is_periodic = 1 1 1");
    }*/

  //
  // read tagging parameters from probin file
  //
  // Tradtionally, the inputs file with ParmParse functionality is handled by C++.  However, a Fortran settings
  // file, by default named probin, can also supply variables.  Mostly used for mesh refinement (tagging) critera
  std::string probin_file("probin");

  ParmParse ppa("amr");
  ppa.query("probin_file",probin_file);

  int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  // use a fortran routine to
  // read in tagging parameters from probin file
  get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

  ParmParse ppn("num");
  
  ppn.get("Strang", advanceOrder);
  amrex::Print() << "Reading Strang splitting order: " << std::endl; 
  if (advanceOrder==1){
    amrex::Print() << "1st order Strang splitting" << std::endl;
    advanceWithChosenUpdateOrder = &CAMReXmp::StrangFirst;
  }
  else if (advanceOrder==2){
    amrex::Print() << "2nd order Strang splitting" << std::endl;
    advanceWithChosenUpdateOrder = &CAMReXmp::StrangSecond;
  }
  else
    amrex::Abort("Please specify a valid Strang splitting order: 1 or 2");
  
  ppn.get("RK", RKOrder);
  amrex::Print() << "Reading RK order: " << std::endl;
  if (RKOrder==1){
    amrex::Print() << "1st order RK splitting" << std::endl;
    RKWithChosenUpdateOrder = &CAMReXmp::RK1;
    tau = 1.0;
  }
  else if (RKOrder==2){
    amrex::Print() << "2nd order RK splitting" << std::endl;
    RKWithChosenUpdateOrder = &CAMReXmp::RK2;
    tau = 0.5;
  }
  else
    amrex::Abort("Please specify a valid RK order: 1 or 2");
  
  ppn.get("Maxwell",MaxwellMethod);
  amrex::Print() << "Reading Maxwell method: " << std::endl;
  if (MaxwellMethod=="HYP"){
    amrex::Print() << "HYP Maxwell treatment" << std::endl;
    MaxwellSolverWithChosenMethod = &CAMReXmp::hyperbolicMaxwellRK;
  }
  else if (MaxwellMethod=="IM"){
    amrex::Print() << "IM Maxwell treatment" << std::endl;
    MaxwellSolverWithChosenMethod = &CAMReXmp::implicitMaxwellRK;
  }
  else
    amrex::Abort("Please specify a valid Maxwell method: HYP or IM");

  ppn.get("source",sourceMethod);
  amrex::Print() << "Reading source method: " << std::endl;
  if (sourceMethod=="EX"){
    amrex::Print() << "EX source treatment" << std::endl;
    sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateEX;
  }
  else if (sourceMethod=="IM"){
    amrex::Print() << "IM source treatment" << std::endl;
    sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateIM;
  }
  else if (sourceMethod=="ANEX"){
    amrex::Print() << "ANEX source treatment" << std::endl;
    sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateANEX;
  }
  else
    amrex::Abort("Please specify a valid source method: EX, IM or ANEX");
}

//
// AMReX has an inbuilt reflux command, but we still have the freedom
// to decide what goes into it (for example, which variables are
// actually refluxed).  This also gives a little flexibility as to
// where flux registers are stored.  In this example, they are stored
// on levels [1,fine] but not level 0.  
//
void
CAMReXmp::reflux ()
{
  BL_ASSERT(level<parent->finestLevel());

  const Real strt = amrex::second();

  // Call the reflux command with the appropriate data.  Because there
  // are no flux registers on the coarse level, they start from the
  // first level.  But the coarse level to the (n-1)^th are the ones
  // that need refluxing, hence the `level+1'.  
  getFluxReg(level+1).Reflux(get_new_data(Phi_Type),1.0,0,0,NUM_STATE,geom);
  
  if (verbose)
  {
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    Real      end    = amrex::second() - strt;
    
    ParallelDescriptor::ReduceRealMax(end,IOProc);
    
    amrex::Print() << "CAMReXmp::reflux() at level " << level 
		   << " : time = " << end << std::endl;
  }
}

//
// Generic function for averaging down - in this case it just makes sure it doesn't happen on the finest level
//
void
CAMReXmp::avgDown ()
{
  if (level == parent->finestLevel())
  {
    return;
  }
  // Can select which variables averaging down will happen on - only one to choose from in this case!
  avgDown(Phi_Type);
}

//
// Setting up the call to the AMReX-implemented average down function
//
void
CAMReXmp::avgDown (int state_indx)
{
  // For safety, again make sure this only happens if a finer level exists
  if (level == parent->finestLevel()) return;

  // You can access data at other refinement levels, use this to
  // specify your current data, and the finer data that is to be
  // averaged down
  CAMReXmp& fine_lev = getLevel(level+1);
  MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
  MultiFab&  S_crse   = get_new_data(state_indx);

  // Call the AMReX average down function:
  // S_fine: Multifab with the fine data to be averaged down
  // S_crse: Multifab with the coarse data to receive the fine data where necessary
  // fine_lev.geom:  Geometric information (cell size etc.) for the fine level
  // geom: Geometric information for the coarse level (i.e. this level)
  // 0: First variable to be averaged (as not all variables need averaging down
  // S_fine.nComp(): Number of variables to average - this can be computed automatically from a multifab
  // refRatio: The refinement ratio between this level and the finer level
  amrex::average_down(S_fine,S_crse,
		      fine_lev.geom,geom,
		      0,S_fine.nComp(),parent->refRatio(level));
}
/*void
CAMReXmp::FillDomainBoundaryFC(Array<MultiFab,AMREX_SPACEDIM>& phi, const Vector<BCRec>& bc)
{
  for (int d = 0; d < AMREX_SPACEDIM; d++)
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
      
      for (MFIter mfi(phi[d], true); mfi.isValid(); ++mfi)
	{
	  
	  const Dim3 lo = lbound(geom.Domain());
	  const Dim3 hi = ubound(geom.Domain());

	  // Indexable arrays for the data, and the directional flux
	  // Based on the corner-centred definition of the flux array, the
	  // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	  const auto& arr = phi[d].array(mfi);

	  std::cout << arr << std::endl;
	  
	  for (int n = 0; n < arr.nComp(); n++)
	    {
	      for (int ng = 1; ng <= NUM_GROW; ng++)
		{
		  if (bc[n].lo(d) == BCType::int_dir){
		    arr(0-ng*iOffset,0-ng*jOffset,
			0-ng*kOffset,n) = arr(hi.x-ng*iOffset,
					      hi.y-ng*jOffset,hi.z-ng*kOffset,n);
		  } else if (bc[n].lo(d) == BCType::foextrap){
		    arr(0-ng*iOffset,0-ng*jOffset,
			0-ng*kOffset,n) = arr(0,0,0,n);
		  } else if (bc[n].lo(d) == BCType::reflect_even){
		    arr(0-ng*iOffset,0-ng*jOffset,
			0-ng*kOffset,n) = arr(0+ng*iOffset,
					      0+ng*jOffset,0+ng*kOffset,n);
		  } else if (bc[n].lo(d) == BCType::reflect_odd){
		    arr(0-ng*iOffset,0-ng*jOffset,
			0-ng*kOffset,n) = -arr(0+ng*iOffset,
					       0+ng*jOffset,0+ng*kOffset,n);
		  } else{
		    amrex::Abort("FC BC not implemented yet!");
		  }

		  if (bc[n].hi(d) == BCType::int_dir){
		    arr(hi.x+ng*iOffset,hi.y+ng*jOffset,
			hi.z+ng*kOffset,n) = arr(0+ng*iOffset,
						 0+ng*jOffset,0+ng*kOffset,n);
		  } else if (bc[n].hi(d) == BCType::foextrap){
		    arr(hi.x+ng*iOffset,hi.y+ng*jOffset,
			hi.z+ng*kOffset,n) = arr(hi.x,hi.y,hi.z,n);
		  } else if (bc[n].hi(d) == BCType::reflect_even){
		    arr(hi.x+ng*iOffset,hi.y+ng*jOffset,
			hi.z+ng*kOffset,n) = arr(hi.x-ng*iOffset,
						 hi.y-ng*jOffset,hi.z-ng*kOffset,n);
		  } else if (bc[n].hi(d) == BCType::reflect_odd){
		    arr(hi.x+ng*iOffset,hi.y+ng*jOffset,
			hi.z+ng*kOffset,n) = -arr(hi.x-ng*iOffset,
						  hi.y-ng*jOffset,hi.z-ng*kOffset,n);
		  } else{
		    amrex::Abort("FC BC not implemented yet!");
		  }
		}
	    }
	}
    }
}
*/
