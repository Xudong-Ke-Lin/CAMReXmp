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

int      CAMReXmp::NUM_STATE       = 18;  // Number of variables in the state
int      CAMReXmp::NUM_GROW        = 3;  // number of ghost cells

int CAMReXmp::NUM_STATE_FLUID = 10;
int CAMReXmp::NUM_STATE_MAXWELL = CAMReXmp::NUM_STATE - NUM_STATE_FLUID;

Vector<BCRec> CAMReXmp::bc;
Vector<BCRec> CAMReXmp::bc_EM;

int  CAMReXmp::max_order = 2;
int CAMReXmp::max_fmg_iter = 0;
Real CAMReXmp::soln_tol = 1e-10;

std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_lobc_X;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_hibc_X;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_lobc_Y;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_hibc_Y;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_lobc_Z;
std::array<LinOpBCType, AMREX_SPACEDIM> CAMReXmp::mlmg_hibc_Z;

sourceFunctionPointer CAMReXmp::sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateIMMidpoint;
fluidSpaceFunctionPointer CAMReXmp::fluidSolverWithChosenSpaceOrder = &CAMReXmp::fluidSolverTVD;
#if (AMREX_SPACEDIM >= 2)
MaxwellSpaceFunctionPointer CAMReXmp::MaxwellSolverWithChosenSpaceOrder = &CAMReXmp::MaxwellSolverFVTDTVD;
#endif

int CAMReXmp::fluidOrder = 2;
int CAMReXmp::MaxwellOrder = 2;
std::string CAMReXmp::MaxwellTimeMethod = "EX";
std::string CAMReXmp::MaxwellDivMethod = "FVTD";
std::string CAMReXmp::sourceMethod = "IM";
int CAMReXmp::projectionStep = 10;

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
			 StateDescriptor::Point,storedGhostZones,NUM_STATE,
			 &cell_cons_interp);
  // from IAMR (https://github.com/AMReX-Codes/IAMR/blob/main/Source/NS_setup.cpp):
  // StateDescriptor::Interval and &node_bilinear_interp
  // from documentation (https://amrex-codes.github.io/amrex/doxygen/AMReX__Interpolater_8cpp.html)
  // &face_linear_interp or &face_divfree_interp
  // can use .ixType() to check the type, xface will give (N,C) -> Nodal in x and center in y
#if (AMREX_SPACEDIM >= 2)
  IndexType xface(IntVect{AMREX_D_DECL(1,0,0)});
  IndexType yface(IntVect{AMREX_D_DECL(0,1,0)});
  IndexType edge(IntVect{AMREX_D_DECL(1,1,0)});
  desc_lst.addDescriptor(EM_X_Type,xface,
			 StateDescriptor::Point,storedGhostZones,6,
			 &face_linear_interp);
  desc_lst.addDescriptor(EM_Y_Type,yface,
			 StateDescriptor::Point,storedGhostZones,6,
			 &face_linear_interp);
  desc_lst.addDescriptor(EM_XY_Type,edge,
			 StateDescriptor::Point,storedGhostZones,6,
			 &node_bilinear_interp);
#endif
  
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

  bc.resize(NUM_STATE);
#if (AMREX_SPACEDIM >= 2)  
  bc_EM.resize(6);
#endif
  
  ParmParse pp;
  pp.query("test",test);
  
  if (test=="BrioWu" || test=="BrioWu2DCart" || (test.find("Toro") != std::string::npos) || test=="cylExp"
      || test=="LiuLax"  || test=="LiuLax2")
    {
      for (int idim = 0; idim < amrex::SpaceDim; ++idim)
	{
	  for (int n = 0; n < NUM_STATE; ++n)
	    {
	      bc[n].setLo(idim, BCType::foextrap); // transmissive
	      bc[n].setHi(idim, BCType::foextrap);
	    }
	}
      // Euler problem in 1D
    } else if (test=="BrioWu1DCyl")
    {

      // top and bottom are periodic (z-axis)
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(1, BCType::int_dir);
	  bc[n].setHi(1, BCType::int_dir);
	}
      
      // left use geometric boundaries
      // reflective even for most variables      
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(0, BCType::reflect_even);
	}
      // reflective odd for vectors in radial and azimuthal direction
      // these are x and z respectively
      bc[MOMX_I].setLo(0, BCType::reflect_odd);
      bc[MOMZ_I].setLo(0, BCType::reflect_odd);
      bc[MOMX_E].setLo(0, BCType::reflect_odd);
      bc[MOMZ_E].setLo(0, BCType::reflect_odd);
      bc[BX].setLo(0, BCType::reflect_odd);
      bc[BZ].setLo(0, BCType::reflect_odd);
      bc[EX].setLo(0, BCType::reflect_odd);	  
      bc[EZ].setLo(0, BCType::reflect_odd);
      
      // right use transmissive
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setHi(0, BCType::foextrap);
	}
    } else if (test=="Colella")
    {
      for (int idim = 0; idim < amrex::SpaceDim; ++idim)
	{
	  bc[0].setLo(idim, BCType::foextrap);
	  bc[0].setHi(idim, BCType::foextrap);
	  bc[1].setLo(idim, BCType::reflect_odd);
	  bc[1].setHi(idim, BCType::reflect_odd);
	  bc[2].setLo(idim, BCType::foextrap);
	  bc[2].setHi(idim, BCType::foextrap);
	  bc[3].setLo(idim, BCType::foextrap);
	  bc[3].setHi(idim, BCType::foextrap);
	  bc[4].setLo(idim, BCType::foextrap);
	  bc[4].setHi(idim, BCType::foextrap);
	}
      // sperical explosion in cylindrical coordinates
    } else if (test=="sphExp")
    {

      // top and bottom are transmissive (z-axis)
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(1, BCType::foextrap);
	  bc[n].setHi(1, BCType::foextrap);
	}
      
      // left use geometric boundaries
      // reflective even for most variables      
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(0, BCType::reflect_even);
	  // transmissive on the right
	  bc[n].setHi(0, BCType::foextrap);
	}
      // reflective odd for vectors in radial and azimuthal direction
      // these are x and z respectively
      bc[MOMX_I].setLo(0, BCType::reflect_odd);
      bc[MOMZ_I].setLo(0, BCType::reflect_odd);
      bc[MOMX_E].setLo(0, BCType::reflect_odd);	  
      bc[MOMZ_E].setLo(0, BCType::reflect_odd);	  
      
      // EM problem in 2D
    } else if (test=="EMwave" || test=="EMwave1d")
    {
      for (int idim = 0; idim < amrex::SpaceDim; ++idim)
	{
	  for (int n = 0; n < NUM_STATE; ++n)
	    {
	      bc[n].setLo(idim, BCType::int_dir); // interior
	      bc[n].setHi(idim, BCType::int_dir);
	    }
	}
    } else if (test=="EMwaveTM")
    {

      // top and bottom are periodic (z-axis)
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(1, BCType::int_dir);
	  bc[n].setHi(1, BCType::int_dir);
	}
      
      // left use geometric boundaries
      // reflective even for most variables      
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(0, BCType::reflect_even);
	}
      // reflective odd for vectors in radial and azimuthal direction
      // these are x and z respectively
      bc[BX].setLo(0, BCType::reflect_odd);
      bc[BZ].setLo(0, BCType::reflect_odd);
      //bc[BZ].setLo(0, BCType::reflect_even);
      bc[EX].setLo(0, BCType::reflect_odd);	  
      bc[EZ].setLo(0, BCType::reflect_odd);
      //bc[EZ].setLo(0, BCType::reflect_even);	  
      
      // right use conducting wall
      // foextrap for most variables      
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setHi(0, BCType::foextrap);
	}
      // reflective odd for normal magnetic field and tangencials electric fields
      // these are Br (Bx) and, Ez (Ey) and Etheta (Ez) respectively
      bc[BX].setHi(0, BCType::reflect_odd);
      bc[EY].setHi(0, BCType::reflect_odd);	  
      bc[EZ].setHi(0, BCType::reflect_odd);	  
            
    } else if (test=="EMwaveTE1d")
    {

      // left use geometric boundaries
      // reflective even for most variables      
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(0, BCType::reflect_even);
	}
      // reflective odd for vectors in radial and azimuthal direction
      // these are x and z respectively
      bc[BX].setLo(0, BCType::reflect_odd);
      bc[BZ].setLo(0, BCType::reflect_odd);
      bc[EX].setLo(0, BCType::reflect_odd);	  
      bc[EZ].setLo(0, BCType::reflect_odd);	  
      
      // right use conducting wall
      // foextrap for most variables      
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setHi(0, BCType::foextrap);
	}
      // reflective odd for normal magnetic field and tangencials electric fields
      // these are Br (Bx) and, Ez (Ey) and Etheta (Ez) respectively
      bc[BX].setHi(0, BCType::reflect_odd);
      bc[EY].setHi(0, BCType::reflect_odd);	  
      bc[EZ].setHi(0, BCType::reflect_odd);	  
            
    } else if (test=="gaussianEM")
    {
      for (int idim = 0; idim < amrex::SpaceDim; ++idim)
	{
	  for (int n = 0; n < NUM_STATE; ++n)
	    {
	      bc[n].setLo(idim, BCType::foextrap);
	      bc[n].setHi(idim, BCType::foextrap);
	    }
	}      
    } else if (test=="OT" || test=="OTideal")
    {
      for (int i = 0; i < amrex::SpaceDim; ++i)
	{
	  for (int n = 0; n < NUM_STATE; ++n)
	    {
	      bc[n].setLo(i, BCType::int_dir); // interior 
	      bc[n].setHi(i, BCType::int_dir);
	    }
	}
    } else if (test=="Harris_sheet_full")
    {
      for (int i = 0; i < amrex::SpaceDim; ++i)
	{
	  for (int n = 0; n < NUM_STATE; ++n)
	    {
	      bc[n].setLo(i, BCType::int_dir); // interior 
	      bc[n].setHi(i, BCType::int_dir);
	    }
	}
    } else if (test=="Harris_sheet")
    {
      for (int n = 0; n < NUM_STATE; ++n)
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
      bc[DIVE].setLo(1, BCType::reflect_odd);
      bc[DIVE].setHi(1, BCType::reflect_odd);

    } else if (test=="blast")
    {
      for (int idim = 0; idim < amrex::SpaceDim; ++idim)
	{
	  for (int n = 0; n < NUM_STATE; ++n)
	    {
	      bc[n].setLo(idim, BCType::foextrap);
	      bc[n].setHi(idim, BCType::foextrap);
	    }
	}      
    } else if (test=="zpinch1d" || test=="zpinch2d" || test=="zpinch2dTrue")
    {

      // top and bottom are periodic (z-axis)
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(1, BCType::int_dir);
	  bc[n].setHi(1, BCType::int_dir);
	}
      
      // left use geometric boundaries
      // reflective even for most variables      
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setLo(0, BCType::reflect_even);
	}
      // reflective odd for vectors in radial and azimuthal direction
      // these are x and z respectively
      bc[MOMX_I].setLo(0, BCType::reflect_odd);
      bc[MOMZ_I].setLo(0, BCType::reflect_odd);
      bc[MOMX_E].setLo(0, BCType::reflect_odd);	  
      bc[MOMZ_E].setLo(0, BCType::reflect_odd);	  
      bc[BX].setLo(0, BCType::reflect_odd);
      bc[BZ].setLo(0, BCType::reflect_odd);
      bc[EX].setLo(0, BCType::reflect_odd);	  
      bc[EZ].setLo(0, BCType::reflect_odd);
      
      // right use conducting wall
      // foextrap for most variables      
      for (int n = 0; n < NUM_STATE; ++n)
	{
	  bc[n].setHi(0, BCType::foextrap);
	}
      // reflective odd for normal magnetic field and tangencials electric fields
      // these are Br (Bx) and, Ez (Ey) and Etheta (Ez) respectively
      bc[MOMX_I].setHi(0, BCType::reflect_odd);
      bc[MOMX_E].setHi(0, BCType::reflect_odd);
      bc[BX].setHi(0, BCType::reflect_odd);
      bc[EY].setHi(0, BCType::reflect_odd);	  
      bc[EZ].setHi(0, BCType::reflect_odd);	  
                  
    } else if (test=="convergence" || test=="convergence2D")
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
#if (AMREX_SPACEDIM >= 2)
  // bc array for the electromagnetic field
  for (int n=BX_LOCAL; n<=EZ_LOCAL; n++)
    {
      bc_EM[n] = bc[n+NUM_STATE_FLUID];
    }
#endif
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

#if (AMREX_SPACEDIM >= 2)  
  // face-cenctred primary variables in 2D
  //desc_lst.setComponent(EM_X_Type, BX_LOCAL, "magxFace", bc_EM[BX_LOCAL],
  //                      StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, BX_LOCAL, "magxFace", bc_EM[BX_LOCAL],
			StateDescriptor::BndryFunc(nullfill));
  //desc_lst.setComponent(EM_Y_Type, BY_LOCAL, "magyFace", bc_EM[BY_LOCAL],
  //                      StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, BY_LOCAL, "magyFace", bc_EM[BY_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, EX_LOCAL, "elecxFace", bc_EM[EX_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, EY_LOCAL, "elecyFace", bc_EM[EY_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  // face-centred moments of the primary variables in 2D
  // used for the divergence-free (second or higher order) reconstruction
  //desc_lst.setComponent(EM_X_Type, BY_LOCAL, "magyFaceX", bc_EM[BY_LOCAL],
  //                      StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, BY_LOCAL, "magyFaceY", bc_EM[BY_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, BZ_LOCAL, "magzFaceX", bc_EM[BZ_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  //desc_lst.setComponent(EM_Y_Type, BX_LOCAL, "magxFaceY", bc_EM[BX_LOCAL],
  //                      StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, BX_LOCAL, "magxFaceX", bc_EM[BX_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, BZ_LOCAL, "magzFaceY", bc_EM[BZ_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, EY_LOCAL, "elecyFaceX", bc_EM[EY_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_X_Type, EZ_LOCAL, "eleczFaceX", bc_EM[EZ_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, EX_LOCAL, "elecxFaceY", bc_EM[EX_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_Y_Type, EZ_LOCAL, "eleczFaceY", bc_EM[EZ_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));

  desc_lst.setComponent(EM_XY_Type, BX_LOCAL, "magxEdge", bc_EM[BX_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_XY_Type, BY_LOCAL, "magyEdge", bc_EM[BY_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_XY_Type, BZ_LOCAL, "magzEdge", bc_EM[BZ_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_XY_Type, EX_LOCAL, "elecxEdge", bc_EM[EX_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_XY_Type, EY_LOCAL, "elecyEdge", bc_EM[EY_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
  desc_lst.setComponent(EM_XY_Type, EZ_LOCAL, "eleczEdge", bc_EM[EZ_LOCAL],
                        StateDescriptor::BndryFunc(nullfill));
#endif
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
//Initialize data on this level after regridding if old level did not previously exist
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
    //MultiFab& S_EM_X_new = get_new_data(EM_X_Type);
    //MultiFab& S_EM_Y_new = get_new_data(EM_Y_Type);

    // See first init function for documentation
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
    //FillCoarsePatch(S_EM_X_new, 0, cur_time, EM_X_Type, 0, 6);
    //FillCoarsePatch(S_EM_Y_new, 0, cur_time, EM_Y_Type, 0, 6);
    
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
  StrangSecond(dx, dt, time);

  MultiFab SNew(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, SNew, NUM_GROW, time+dt, Phi_Type, 0, NUM_STATE);
#if (AMREX_SPACEDIM >= 2)  
  Array<MultiFab,AMREX_SPACEDIM> S_EMNew;
  S_EMNew[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EMNew[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
  S_EMNew[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EMNew[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);  
#endif
  
  ParmParse pp;
  Real stop_time;
  pp.query("stop_time",stop_time);
  // conditional used in case of hyperbolic divergence cleaning, because DIVB and DIVE will be divergence cleaning parameters
  if (cur_time>=stop_time)
    {
      amrex::Print() << "Last time step: calculating divergence errors" << std::endl;
      // Loop over all the patches at this level
      for (MFIter mfi(SNew, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();

	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  // Indexable arrays for the data, and the directional flux
	  // Based on the vertex-centred definition of the flux array, the
	  // data array runs from e.g. [0,N] and the flux array from [0,N+1]
	  const auto& arr = SNew.array(mfi);
        
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {	    	    
		      // Initialise to zero
		      arr(i,j,k,DIVB) = 0.0;
		      arr(i,j,k,DIVE) = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
		    }
		}
	    }    
	}
      // Compute divergence errors
      for (int d = 0; d < amrex::SpaceDim ; d++)
	{
	  const int iOffset = ( d == 0 ? 1 : 0);
	  const int jOffset = ( d == 1 ? 1 : 0);
	  const int kOffset = ( d == 2 ? 1 : 0);

	  for (MFIter mfi(SNew, true); mfi.isValid(); ++mfi)
	    {
	      const Box& bx = mfi.tilebox();
	  
	      const Dim3 lo = lbound(bx);
	      const Dim3 hi = ubound(bx);
	
	      Array4<Real> arr = SNew.array(mfi);
#if (AMREX_SPACEDIM >= 2)	      
	      //Array4<Real> arrEM = S_EMNew[d].array(mfi);	
#endif
	      
	      //for(int k = lo.z; k <= hi.z+kOffset; k++)
	      for(int k = lo.z; k <= hi.z; k++)
		{
		  //for(int j = lo.y; j <= hi.y+jOffset; j++)
		  for(int j = lo.y; j <= hi.y; j++)
		    {
		      //for(int i = lo.x; i <= hi.x+iOffset; i++)
		      for(int i = lo.x; i <= hi.x; i++)
			{
			  //arr(i,j,k,DIVB) += (arrEM(i+iOffset,j+jOffset,k,BX_LOCAL+d)-arrEM(i,j,k,BX_LOCAL+d))/ dx[d];
			  arr(i,j,k,DIVB) += (arr(i+iOffset,j+jOffset,k,BX+d)-arr(i-iOffset,j-jOffset,k,BX+d))/(2.0*dx[d]);
			  // substract because divEerror was set to be the charge density
			  //arr(i,j,k,DIVE) -= (arrEM(i+iOffset,j+jOffset,k,EX_LOCAL+d)-arrEM(i,j,k,EX_LOCAL+d))/dx[d];
			  arr(i,j,k,DIVE) -= (arr(i+iOffset,j+jOffset,k,EX+d)-arr(i-iOffset,j-jOffset,k,EX+d))/(2.0*dx[d]);
			}
		    }
		}      
	    }
	}
      // copy divergence errors
      MultiFab::Copy(S_new, SNew, DIVB, DIVB, 2, 0);
    }
  
  /*
  // Density errors for convergence problem
  for (MFIter mfi(S_EMNew[0], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      //const auto& arr = S_source.array(mfi);
      const auto& arrEM = S_EMNew[0].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
	      const Real y = prob_lo[1] + (double(j)+0.5) * dx[1];
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  const Real x = prob_lo[0] + (double(i)) * dx[0];
  		  // forcing term for Ex
  		  //arrEM(i,j,k,EX_LOCAL) += -dt*(2.0+std::sin(2.0*M_PI*(x-cur_time)))/(lambda_d*lambda_d*l_r);
		  arrEM(i,j,k,EX_LOCAL) += -dt*c/std::sqrt(2.0)*(2.0+std::sin(2.0*M_PI*(x+y-std::sqrt(2.0)*c*cur_time)))/(lambda_d*lambda_d*l_r);		  
  		}
  	    }
  	}      
    }
  for (MFIter mfi(S_EMNew[1], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      //const auto& arr = S_source.array(mfi);
      const auto& arrEM = S_EMNew[1].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
	      const Real y = prob_lo[1] + (double(j)) * dx[1];
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  const Real x = prob_lo[0] + (double(i)+0.5) * dx[0];
  		  // forcing term for Ey
		  arrEM(i,j,k,EY_LOCAL) += -dt*c/std::sqrt(2.0)*(2.0+std::sin(2.0*M_PI*(x+y-std::sqrt(2.0)*c*cur_time)))/(lambda_d*lambda_d*l_r);
  		}
  	    }
  	}      
    }
  for (MFIter mfi(SNew, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      //const auto& arr = S_source.array(mfi);
      const auto& arr = SNew.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
	      const Real y = prob_lo[1] + (double(j)+0.5) * dx[1];
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  const Real x = prob_lo[0] + (double(i)+0.5) * dx[0];
  		  // forcing term for Ex
  		  //arr(i,j,k,EX) += -dt*(2.0+std::sin(2.0*M_PI*(x-cur_time)))/(lambda_d*lambda_d*l_r);
		  // in 2D
		  arr(i,j,k,EX) += -dt*c/std::sqrt(2.0)*(2.0+std::sin(2.0*M_PI*(x+y-std::sqrt(2.0)*c*cur_time)))/(lambda_d*lambda_d*l_r);
		  arr(i,j,k,EY) += -dt*c/std::sqrt(2.0)*(2.0+std::sin(2.0*M_PI*(x+y-std::sqrt(2.0)*c*cur_time)))/(lambda_d*lambda_d*l_r);
		  // Real rho_exact = 2.0+std::sin(2.0*M_PI*(x-cur_time));
		  Real rho_exact = 2.0+std::sin(2.0*M_PI*(x+y-std::sqrt(2.0)*c*cur_time));
		  arr(i,j,k,DIVB) = arr(i,j,k,0)-rho_exact;
  		}
  	    }
  	}      
    }
  */
  /*
  // Bz and Ey errors for EM wave problem
  ParmParse pp;
  Real stop_time;
  pp.query("stop_time",stop_time);
  if (cur_time>=stop_time)
    {
      amrex::Print() << "Last time step: calculating EM errors" << std::endl;
      for (MFIter mfi(SNew, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);
	  
	  Array4<Real> arr = SNew.array(mfi);
	  
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  const Real y = prob_lo[1] + (double(j)+0.5) * dx[1];
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      const Real x = prob_lo[0] + (double(i)+0.5) * dx[0];
		      arr(i,j,k,DIVB) = arr(i,j,k,BZ)-std::cos(2.0*M_PI*(x+y-std::sqrt(2.0)*c*cur_time));	 
		      arr(i,j,k,DIVE) = (arr(i,j,k,EY)-c*std::cos(2.0*M_PI*(x+y-std::sqrt(2.0)*c*cur_time))/std::sqrt(2.0)); 
		    }
		}
	    }      
	}  
      // copy errors
      MultiFab::Copy(S_new, SNew, DIVB, DIVB, 2, 0);
    }
  // We need to compute boundary conditions again after each update       
  //Sborder.FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries 
  //FillDomainBoundary(Sborder, geom, bc);
  */
  // Btheta (Bz) and Er (Ex) errors for EM circular waveguide problem
  /*ParmParse pp;
  Real stop_time;
  pp.query("stop_time",stop_time);
  //if (cur_time>=stop_time)
    {
      amrex::Print() << "Last time step: calculating EM errors" << std::endl;
      for (MFIter mfi(SNew, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);
	  
	  Array4<Real> arr = SNew.array(mfi);
	  
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  const Real y = prob_lo[1] + (double(j)+0.5) * dx[1];
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      const Real x = prob_lo[0] + (double(i)+0.5) * dx[0];

		      using namespace std::complex_literals;
		      
		      // cylindrical coordinates
		      const Real rCyl = x;
		      const Real zCyl = y;

		      // parameters
		      Real xi0 = 2.40482555769577;
		      Real kr = xi0;
		      Real kz = M_PI;
		      Real omega = std::sqrt(kr*kr + kz*kz);

		      // function
		      auto Bessel0ExpFunc = std::cyl_bessel_j(0,kr*rCyl)*std::exp(1i*kz*zCyl)*std::exp(-1i*omega*cur_time);
		      auto Bessel1ExpFunc = std::cyl_bessel_j(1,kr*rCyl)*std::exp(1i*kz*zCyl)*std::exp(-1i*omega*cur_time);

		      // negative sign for theta components in cyl. coord.
		      arr(i,j,k,DIVB) = arr(i,j,k,BZ)+std::real(-1i*omega*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
		      //arr(i,j,k,DIVE) = arr(i,j,k,EX)-std::real(-1i*kz*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
		      arr(i,j,k,DIVE) = arr(i,j,k,EY)-std::real(Bessel0ExpFunc);

		    }
		}
	    }      
	}  
      // copy errors
      MultiFab::Copy(S_new, SNew, DIVB, DIVB, 2, 0);
    }
  */
  // Btheta (Bz) and Ez (Ex) exact sols for cylindrical Brio-Wu test with EM circular waveguide
  /*ParmParse pp;
  Real stop_time;
  pp.query("stop_time",stop_time);
  //if (cur_time>=stop_time)
    {
      //amrex::Print() << "Last time step: calculating EM errors" << std::endl;
      for (MFIter mfi(SNew, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);
	  
	  Array4<Real> arr = SNew.array(mfi);
	  
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  const Real y = prob_lo[1] + (double(j)+0.5) * dx[1];
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      const Real x = prob_lo[0] + (double(i)+0.5) * dx[0];

		      using namespace std::complex_literals;
		      
		      // cylindrical coordinates
		      const Real rCyl = x;

		      // parameters
		      Real xi0 = 2.40482555769577;
		      Real omega = xi0/(0.1*2.0*M_PI);

		      // function
		      auto Bessel0SinFunc = std::cyl_bessel_j(0,omega*rCyl)*std::sin(omega*c*cur_time);
		      auto Bessel1CosFunc = std::cyl_bessel_j(1,omega*rCyl)*std::cos(omega*c*cur_time);

		      // negative sign for theta components in cyl. coord.
		      arr(i,j,k,DIVB) = Bessel1CosFunc;		     
		      arr(i,j,k,DIVE) = c*Bessel0SinFunc;

		    }
		}
	    }      
	}  
      // copy errors
      MultiFab::Copy(S_new, SNew, DIVB, DIVB, 2, 0);
      //MultiFab::Copy(S_new, SNew, EX, EX, 1, 0);
    }  
  */
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
  
  //MultiFab::Copy(S_new, SNew, 0, 0, NUM_STATE, 0);
  // !!
  // Coment out for HDC
  //MultiFab::Copy(S_EM_X_new, S_EMNew[0], 0, 0, 6, 0);
#if (AMREX_SPACEDIM >= 2) 
  //MultiFab::Copy(S_EM_Y_new, S_EMNew[1], 0, 0, 6, 0);
#endif
  
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
  /*if (do_reflux) {
    if (current) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
	current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
    }
    if (fine) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
	fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
    }
  }
  */
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
  
  // State with ghost cells - this is used to compute fluxes and perform the update.
  MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
  // See init function for details about the FillPatch function
  FillPatch(*this, Sborder, NUM_GROW, state[Phi_Type].curTime(), Phi_Type, 0, NUM_STATE);
  Sborder.FillBoundary(geom.periodicity());
  // added by 2020D            
  // Fill non-periodic physical boundaries     
  FillDomainBoundary(Sborder, geom, bc);
  
  bool fluidconstraint = (fluidOrder != 0 ? true : false);
  bool EMconstraint = (MaxwellOrder != 0 ? true : false);
  
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

  // if electron mass smaller than ion mass, use electron velocities
  const int fluid = ( m>1.0 ? NUM_STATE_FLUID/2 : 0);
  
  MFIter::allowMultipleMFIters(true);

  for(unsigned int d = 0; d < amrex::SpaceDim; ++d)
  {
    Vector<Real> c_array {0.0};
    
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
		  Real v = arr(i,j,k,MOMX_I+d+fluid)/arr(i,j,k,RHO_I+fluid);
		  Vector<Real> u_i = get_data_zone(arr,i,j,k,fluid,NUM_STATE_FLUID/2);
		  Real a = get_speed(u_i);
		  c_array.push_back(std::abs(v)+a);
		}
	    }
	}
    }
    c_h = *std::max_element(c_array.begin(), c_array.end());

    // dt does not need to resolve plasma and cyclotron frequencies
    // because usually an implicit souce term update is used
    if (EMconstraint && fluidconstraint && MaxwellTimeMethod=="EX")
      dt_est = std::min(dt_est, dx[d]/std::max(c_h, c));
    else if (EMconstraint && fluidconstraint && (MaxwellTimeMethod=="EXsubcycling" || MaxwellTimeMethod=="IM"))
      dt_est = std::min(dt_est, dx[d]/c_h);
    else if (EMconstraint && !fluidconstraint)
      dt_est = std::min(dt_est, dx[d]/c);
    else
      dt_est = std::min(dt_est, dx[d]/c_h);
        
  }  
  
  // Ensure that we really do have the minimum across all processors
  ParallelDescriptor::ReduceRealMin(dt_est);
  ParallelDescriptor::ReduceRealMax(c_h);

  dt_est *= cfl;
  amrex::Print() << c_h << " " << dt_est << std::endl;
  
  if (c_h>c && EMconstraint && fluidconstraint)
    {
      amrex::Abort("Fluid velocity is higher than speed of light!");
      amrex::Print() << "Fluid velocity is higher than speed of light!" << std::endl;
      ParmParse pp;
      Real stop_time;
      pp.query("stop_time",stop_time);
      //dt_est = stop_time - cur_time;      
    }
  if (dt_est<1e-8)
    {
      amrex::Abort("Too small time step!");
      amrex::Print() <<	"Too small time step!"	<< std::endl;
      ParmParse pp;
      Real stop_time;
      pp.query("stop_time",stop_time);
      dt_est = stop_time - cur_time;
    }

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
      
  ppn.get("fluid",fluidOrder);
  amrex::Print() << "Reading fluid method: " << std::endl;
  if (fluidOrder==0){
    amrex::Print() << "no fluid solver" << std::endl;
    fluidSolverWithChosenSpaceOrder = &CAMReXmp::fluidSolverVoid;
  }
  else if (fluidOrder==2){
    amrex::Print() << "fluid 2nd order TVD" << std::endl;
    fluidSolverWithChosenSpaceOrder = &CAMReXmp::fluidSolverAdvTVD;//fluidSolverTVD;
  }
  else if (fluidOrder==3){
    amrex::Print() << "fluid 3rd order WENO" << std::endl;
    fluidSolverWithChosenSpaceOrder = &CAMReXmp::fluidSolverWENO;
  }
  else
    amrex::Abort("Please specify a valid fluid order: 0, 2 or 3");

#if (AMREX_SPACEDIM >= 2)  
  ppn.get("Maxwell",MaxwellOrder);
  amrex::Print() << "Reading Maxwell method: " << std::endl;
  if (MaxwellOrder==0){
    amrex::Print() << "no Maxwell solver" << std::endl;
    MaxwellSolverWithChosenSpaceOrder = &CAMReXmp::MaxwellSolverFVTDVoid;
  }
  else if (MaxwellOrder==2){
    amrex::Print() << "Maxwell 2nd order TVD" << std::endl;
    MaxwellSolverWithChosenSpaceOrder = &CAMReXmp::MaxwellSolverFVTDTVD;
  }
  else if (MaxwellOrder==3){
    amrex::Print() << "Maxwell 3rd order WENO" << std::endl;
    MaxwellSolverWithChosenSpaceOrder = &CAMReXmp::MaxwellSolverFVTDWENO;
  }
  else
    amrex::Abort("Please specify a valid Maxwell order: 2 or 3");
#endif
  
  ppn.get("MaxwellTimeMethod",MaxwellTimeMethod);
  amrex::Print() << "Reading Maxwell time stepping method: " << std::endl;
  if (MaxwellTimeMethod=="EX"){
    amrex::Print() << "explicit Maxwell time stepping" << std::endl;
  }
  else if (MaxwellTimeMethod=="EXsubcycling"){
    amrex::Print() << "explicit Maxwell time stepping with sub-cycling" << std::endl;
  }
  else if (MaxwellTimeMethod=="IM"){
    amrex::Print() << "implicit Maxwell time stepping" << std::endl;
  }
  else
    amrex::Abort("Please specify a valid Maxwell time stepping method: EX, EXsubcycling or IM");

  ppn.get("MaxwellDivMethod",MaxwellDivMethod);
  amrex::Print() << "Reading Maxwell divergence cleaning method: " << std::endl;
  if (MaxwellDivMethod=="FVTD"){
    amrex::Print() << "finite-volume time-domain" << std::endl;
    if (MaxwellTimeMethod!="EX" && MaxwellTimeMethod!="EXsubcycling")
      amrex::Abort("MaxwellDivMethod=FVTD only works with MaxwellTimeMethod=EX or EXsubcycling");
  }
  else if (MaxwellDivMethod=="FDTD"){
    amrex::Print() << "finite-difference time-domain" << std::endl;
    if (MaxwellTimeMethod!="IM")
      amrex::Abort("MaxwellDivMethod=FDTD only works with MaxwellTimeMethod=IM");
  }
  else if (MaxwellDivMethod=="HDC"){
    amrex::Print() << "hyperbolic divergence cleaning" << std::endl;
    amrex::Print() << "reading cb and ce" << std::endl;
    ppn.get("cb",cb);
    ppn.get("ce",ce);
    if (MaxwellTimeMethod!="EX")
      amrex::Abort("MaxwellDivMethod=HDC only works with MaxwellTimeMethod=EX");
  }
  else if (MaxwellDivMethod=="NONE"){
  }
  else
    amrex::Abort("Please specify a valid Maxwell divergence cleaning method: FVTD, FDTD or HDC");

  if (amrex::SpaceDim==1 && (MaxwellDivMethod=="FVTD" || MaxwellDivMethod=="FDTD"))
    amrex::Abort("FVTD and FDTD only work for non-1D setup");
  
  ppn.get("source",sourceMethod);
  amrex::Print() << "Reading source method: " << std::endl;
  if (sourceMethod=="no"){
    amrex::Print() << "no source treatment" << std::endl;
    sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateVoid;
  }
  /*if (sourceMethod=="EX"){
    amrex::Print() << "EX source treatment" << std::endl;
    sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateEX;
  }
  */
  else if (sourceMethod=="IM"){
    amrex::Print() << "IM source treatment" << std::endl;
    sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateIMMidpoint;
  }
  /*else if (sourceMethod=="ANEX"){
    amrex::Print() << "ANEX source treatment" << std::endl;
    sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateANEX;
    }*/
  else if (sourceMethod=="STIFF"){
    amrex::Print() << "STIFF source treatment" << std::endl;
    sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateStiff;
  }
  else if (sourceMethod=="EXACT"){
    amrex::Print() << "EXACR source treatment" << std::endl;
    sourceUpdateWithChosenMethod = &CAMReXmp::sourceUpdateExact;
  }
  else
    amrex::Abort("Please specify a valid source method: EX, IM, ANEX or STIFF");

  ppn.get("projectionStep", projectionStep);
  amrex::Print() << "Reading steps per projection: " << std::endl; 
  
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
