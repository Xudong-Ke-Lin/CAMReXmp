#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::StrangFirst(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  MultiFab Sprev(grids, dmap, NUM_STATE, NUM_GROW);
  if (geom.Coord()==1 && MaxwellMethod=="IM")
    {
      MultiFab::Copy(Sprev, Sborder, 0, 0, NUM_STATE, NUM_GROW); 
    }
  
  (this->*MaxwellSolverWithChosenMethod)(Sborder,fluxes,dx,dt);
  //RK1(Sborder, fluxes, dx, dt, &CAMReXmp::hyperbolicMaxwellSolver, BX, NUM_STATE_MAXWELL);
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::MaxwellSolver, BX, NUM_STATE);
  //implicitMaxwellSolver(Sborder, dx, dt);

  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::fluidSolver, 0, NUM_STATE_FLUID);
  RK1(Sborder, fluxes, dx, dt, &CAMReXmp::fluidSolver, 0, NUM_STATE_FLUID); 

  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::sourceUpdate, 0, NUM_STATE);
  RK1(Sborder, fluxes, dx, dt, &CAMReXmp::sourceUpdate, 0, NUM_STATE);

  if (geom.Coord()==1 && MaxwellMethod=="IM")
    {
      cylSourceUpdateImplicitRK2(Sborder, Sprev, dx, dt);
    }
}

void CAMReXmp::StrangSecond(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  
  (this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, 0.5*dt, &CAMReXmp::sourceUpdate, 0, NUM_STATE);

  RK1(Sborder, fluxes, dx, dt, &CAMReXmp::hyperbolicMaxwellSolver, BX, NUM_STATE_MAXWELL);
  //(this->*MaxwellSolverWithChosenMethod)(Sborder,fluxes,dx,dt);   
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::MaxwellSolver, BX, NUM_STATE);

  RK1(Sborder, fluxes, dx, dt, &CAMReXmp::fluidSolver, 0, NUM_STATE_FLUID);
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::fluidSolver, 0, NUM_STATE_FLUID);
  
  (this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, 0.5*dt, &CAMReXmp::sourceUpdate, 0, NUM_STATE);

}
