#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::StrangFirst(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, const Real* dx, Real time, Real dt)
{
  //MultiFab Sprev(grids, dmap, NUM_STATE, NUM_GROW);
  /*if (geom.Coord()==1 && MaxwellMethod=="IM")
    {
      MultiFab::Copy(Sprev, Sborder, 0, 0, NUM_STATE, NUM_GROW); 
    }
  */
  //(this->*MaxwellSolverWithChosenMethod)(Sborder,fluxes,dx,dt);
  //RK1(Sborder, fluxes, dx, dt, &CAMReXmp::hyperbolicMaxwellSolver, BX, NUM_STATE_MAXWELL);
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::MaxwellSolver, BX, NUM_STATE);
  //implicitMaxwellSolver(Sborder, dx, dt);

  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::fluidSolver, 0, NUM_STATE_FLUID);
  //RK1(Sborder, fluxes, dx, dt, &CAMReXmp::fluidSolver, 0, NUM_STATE_FLUID); 
  
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::sourceUpdate, 0, NUM_STATE);
  //RK1(Sborder, fluxes, dx, dt, &CAMReXmp::sourceUpdate, 0, NUM_STATE);
  /*
  if (geom.Coord()==1 && MaxwellMethod=="IM")
    {
      cylSourceUpdateImplicitRK2(Sborder, Sprev, dx, dt);
    }
  */
  (this->*RKWithChosenUpdateOrder)(S_dest, S_source, fluxes, S_EM_dest, S_EM_source, fluxesEM, dx, time, dt);
  
  sourceUpdate(S_dest, fluxes, dx, dt);

}

void CAMReXmp::StrangSecond(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, const Real* dx, Real time, Real dt)
{
  
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, 0.5*dt, &CAMReXmp::sourceUpdate, 0, NUM_STATE);

  //RK1(Sborder, fluxes, dx, dt, &CAMReXmp::hyperbolicMaxwellSolver, BX, NUM_STATE_MAXWELL);
  //(this->*MaxwellSolverWithChosenMethod)(Sborder,fluxes,dx,dt);   
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::MaxwellSolver, BX, NUM_STATE);

  //RK1(Sborder, fluxes, dx, dt, &CAMReXmp::fluidSolver, 0, NUM_STATE_FLUID);
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::fluidSolver, 0, NUM_STATE_FLUID);
  
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, 0.5*dt, &CAMReXmp::sourceUpdate, 0, NUM_STATE);

  sourceUpdate(S_source, fluxes, dx, 0.5*dt);

  (this->*RKWithChosenUpdateOrder)(S_dest, S_source, fluxes, S_EM_dest, S_EM_source, fluxesEM, dx, time, dt);
  
  sourceUpdate(S_dest, fluxes, dx, 0.5*dt);
}
