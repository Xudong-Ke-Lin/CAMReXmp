#ifndef _CAMReXmp_H_
#define _CAMReXmp_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>

#include <memory>
#include <iostream>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif

#include "utils.H"

#include <AMReX_BLFort.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Utility.H>
#include <AMReX_MLLinOp.H>

enum StateType { Phi_Type = 0, EM_X_Type = 1, EM_Y_Type = 2, 
		 NUM_STATE_TYPE };
enum StateVariable {
		    RHO_I = 0, MOMX_I = 1, MOMY_I = 2, MOMZ_I = 3, ENER_I = 4,
		    RHO_E = 5, MOMX_E = 6, MOMY_E = 7, MOMZ_E = 8, ENER_E = 9,
		    BX = 10, BY = 11, BZ = 12, EX = 13, EY = 14, EZ = 15,
		    DIVB = 16, DIVE = 17
};
enum StateVariableMaxwell {
			   BX_LOCAL = 0, BY_LOCAL = 1, BZ_LOCAL = 2, EX_LOCAL = 3, EY_LOCAL = 4, EZ_LOCAL = 5
};


class CAMReXmp;

typedef void (CAMReXmp::*advanceFunctionPointer)(MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const Real*, Real);
typedef void (CAMReXmp::*StrangFunctionPointer)(MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const Real*, Real);
typedef void (CAMReXmp::*RKFunctionPointer)(MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const Real*, Real,
					    void (CAMReXmp::*)(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const double*, double), int, int);
typedef void (CAMReXmp::*MaxwellFunctionPointer)(MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const Real*, Real);
typedef void (CAMReXmp::*sourceFunctionPointer)(Array4<Real>&, int, int, int, Real);

//
// AmrLevel-derived class for hyperbolic conservation equations
//

class CAMReXmp
    :
    public amrex::AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    CAMReXmp ();
    //
    //The basic constructor.
    //
    CAMReXmp (amrex::Amr&     papa,
	         int             lev,
                 const amrex::Geometry& level_geom,
                 const amrex::BoxArray& bl,
                 const amrex::DistributionMapping& dm,
                 amrex::Real            time);
    //
    //The destructor.
    //
    virtual ~CAMReXmp () override;
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (amrex::Amr&   papa,
                          std::istream& is,
			  bool          bReadSpecial = false) override;

    virtual void checkPoint (const std::string& dir,
			     std::ostream&      os,
			     amrex::VisMF::How  how = amrex::VisMF::NFiles,
			     bool               dump_old = true) override;

    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                amrex::VisMF::How  how) override;
    //
    //Define data descriptors.
    //
    static void variableSetUp ();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData () override;
    //
    //Initialize data on this level from another CAMReXmp (during regrid).
    //
    virtual void init (amrex::AmrLevel& old) override;
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () override;
    //
    //Advance grids at this level in time.
    //
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;
    //
    //Estimate time step.
    //
    amrex::Real estTimeStep (amrex::Real dt_old);
    //
    //Compute initial time step.
    //
    amrex::Real initialTimeStep ();
    //
    //Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   amrex::Vector<int>&           n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                  stop_time) override;
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               amrex::Vector<int>&           n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                  stop_time,
                               int                   post_regrid_flag) override;
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration) override;

    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase, int new_finest) override;
    //
    //Do work after a restart().
    //
    virtual void post_restart () override;
    //
    //Do work after init().
    //
    virtual void post_init (amrex::Real stop_time) override;
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           amrex::Real         time,
			   int          n_error_buf = 0, int ngrow = 0) override;

    static int  NUM_STATE;
    static int  NUM_GROW;
  
  void sourceUpdate(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt);
  void sourceUpdateEX(Array4<Real>& arr, int i, int j, int k, Real dt);
  void sourceUpdateIM(Array4<Real>& arr, int i, int j, int k, Real dt);
  void sourceUpdateANEX(Array4<Real>& arr, int i, int j, int k, Real dt);
  void cylSourceUpdate(Array4<Real>& arr, int i, int j, int k, Real dt, Real x);
  void cylSourceUpdateImplicitRK2(MultiFab& Sborder, MultiFab& Sold, const Real* dx, Real dt);

  static int NUM_STATE_FLUID;
  static int NUM_STATE_MAXWELL;

  static int advanceOrder;
  static int RKOrder;
  static std::string MaxwellMethod;
  static std::string sourceMethod;
  
protected:

    static void read_params ();

    /**
     * Function to read tagging parameters from input file.
     * See Tagging_params.cpp for implementation.
     */
    //static void get_tagging_params ();

    // inline functioned implemented below
    CAMReXmp& getLevel (int lev);
    amrex::FluxRegister& getFluxReg ();
    amrex::FluxRegister& getFluxReg (int lev);

    void reflux ();

    void avgDown ();

    void avgDown (int state_indx);

  void setDomainBC (std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
		    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc, int index);  

    //
    // The data.
    //
    amrex::FluxRegister*        flux_reg;
    //
    // Static data members.
    //
    static int          verbose;
    static amrex::Real  cfl;
    static int          do_reflux;

    static Vector<BCRec> bc;
    static Vector<BCRec> bc_EM;
  
  static advanceFunctionPointer advanceWithChosenUpdateOrder;
  static RKFunctionPointer RKWithChosenUpdateOrder;
  static MaxwellFunctionPointer MaxwellSolverWithChosenMethod;
  static sourceFunctionPointer sourceUpdateWithChosenMethod;  

  using functionPointer = void (CAMReXmp::*)(MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const Real*, Real);
  //using RKfunctionPointer = void (CAMReXmp::*)(MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const Real*, Real, functionPointer);
  //using sourcefunctionPointer = void (CAMReXmp::*)(MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const Real*, Real, functionPointer);

  // Parameters for the implicit Solver
  static int  max_order;
  static int max_fmg_iter;
  static Real soln_tol;
  static Real tau;
  static amrex::MultiFab acoef;
  static std::array<amrex::MultiFab, BL_SPACEDIM> bcoeffs;
  
  // Boundary conditions for MLABecLaplacian
  static std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_lobc_X;
  static std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_hibc_X;
  static std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_lobc_Y;
  static std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_hibc_Y;
  static std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_lobc_Z;
  static std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_hibc_Z;


private:

  void StrangFirst(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real);
  void StrangSecond(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real);
  //void RK1(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real,
  //	   std::function<void (amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real)>);
  void RK1(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real, functionPointer, int, int);
  void RK2(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real, functionPointer, int, int);
  void linearCombination(MultiFab&, MultiFab&, Real, MultiFab&, Real, int, int);
  void fluidSolver(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real);  
  void hyperbolicMaxwellRK(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real);
  void hyperbolicMaxwellSolver(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real);
  void implicitMaxwellRK(amrex::MultiFab&, amrex::MultiFab (&)[AMREX_SPACEDIM], const amrex::Real*, amrex::Real);
  void implicitMaxwellSolverSetUp();  
  void implicitMaxwellSolver(amrex::MultiFab&, const amrex::Real*, amrex::Real);
  void implicitMaxwellSolverElecFieldUpdate(amrex::MultiFab&, const amrex::MultiFab&, const amrex::MultiFab&, const amrex::Real*, amrex::Real);
  void hyperbolicMaxwellSolverDivFree(Array<MultiFab,AMREX_SPACEDIM>&, MultiFab&, MultiFab&, MultiFab (&)[AMREX_SPACEDIM], const Real*, Real);
  void hyperbolicMaxwellSolverDivFreeBalsara(Array<MultiFab,AMREX_SPACEDIM>&, MultiFab&, MultiFab&, MultiFab (&)[AMREX_SPACEDIM], MultiFab&, const Real*, Real);
  void MaxwellSolverDivFree(Array<MultiFab,AMREX_SPACEDIM>&, MultiFab&, MultiFab&, MultiFab (&)[AMREX_SPACEDIM], MultiFab&, const Real*, Real);
  void hyperbolicMaxwellSolverDivFreeSubCycle(Array<MultiFab,AMREX_SPACEDIM>&, MultiFab&, MultiFab&, MultiFab (&)[AMREX_SPACEDIM], const Real*, Real);
  void hyperbolicMaxwellSolverDivFreeSubCycleIneff(Array<MultiFab,AMREX_SPACEDIM>&, MultiFab&, MultiFab&, MultiFab (&)[AMREX_SPACEDIM], MultiFab&, MultiFab (&)[AMREX_SPACEDIM], const Real*, Real);
  void computeCurrentUpdate(amrex::MultiFab&, const amrex::MultiFab&, const amrex::Real*, amrex::Real);
  void computeAlpha(amrex::MultiFab&);
  void computeBeta(std::array<amrex::MultiFab,AMREX_SPACEDIM>&);
  void computeRhs(std::array<MultiFab, 3>&, amrex::MultiFab&, amrex::MultiFab&, const amrex::Real*, amrex::Real);
  void printComponents(amrex::MultiFab&);

  void computeRhsX(amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&, amrex::Real, amrex::Real,
		   amrex::Real, amrex::Real);
  void computeRhsY(amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&, amrex::Real, amrex::Real,
                   amrex::Real, amrex::Real);
  void computeRhsZ(amrex::MultiFab&, amrex::MultiFab&, amrex::MultiFab&, amrex::Real, amrex::Real,
		   amrex::Real, amrex::Real);


};    

/**
 * Boundary condition function that does nothing.
 * See bc_nullfill.cpp for implementation.
 */
void nullfill (amrex::Box const& bx, amrex::FArrayBox& data,
               const int dcomp, const int numcomp,
               amrex::Geometry const& geom, const amrex::Real time,
               const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
               const int scomp);



//
// Inlines.
//

inline
CAMReXmp&
CAMReXmp::getLevel (int lev)
{
    return *(CAMReXmp *) &parent->getLevel(lev);
}

inline
amrex::FluxRegister&
CAMReXmp::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
amrex::FluxRegister&
CAMReXmp::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}
#endif /*_CAMReXmp_H_*/
