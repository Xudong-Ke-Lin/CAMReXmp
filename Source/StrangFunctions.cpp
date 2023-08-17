#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::StrangZero(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, const Real* dx, Real time, Real dt)
{

  (this->*RKWithChosenUpdateOrder)(S_dest, S_source, fluxes, S_EM_dest, S_EM_source, fluxesEM, dx, time, dt);
 
}

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
  
  (this->*RKWithChosenUpdateOrder)(S_dest, S_source, fluxes, S_EM_dest, S_EM_source, fluxesEM, dx, time, dt);
  implicitYeeMaxwellSolver(S_EM_dest,S_EM_source,S_dest,S_source,dx,dt,time);

  return;

  sourceUpdate(S_source, fluxes, dx, 0.5*dt);

  if (sourceMethod=="IM")
    {
      // Update face-centred EM fields
      for (int d = 0; d < amrex::SpaceDim ; d++)   
	{

	  const int iOffset = ( d == 0 ? 1 : 0);
	  const int jOffset = ( d == 1 ? 1 : 0);
	  const int kOffset = ( d == 2 ? 1 : 0);
    
	  // Loop over all the patches at this level
	  for (MFIter mfi(S_EM_source[d], true); mfi.isValid(); ++mfi)
	    {
	      const Box& bx = mfi.tilebox();

	      const Dim3 lo = lbound(bx);
	      const Dim3 hi = ubound(bx);

	      // Indexable arrays for the data, and the directional flux
	      // Based on the corner-centred definition of the flux array, the
	      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	      //const auto& arr = S_dest.array(mfi);
	      const auto& arrEM = S_EM_source[d].array(mfi);
	      const auto& arr = S_source.array(mfi);

	      const Dim3 hiDomain = ubound(geom.Domain());
      
	      for(int k = lo.z; k <= hi.z; k++)
		{
		  for(int j = lo.y; j <= hi.y; j++)
		    {
		      for(int i = lo.x; i <= hi.x; i++)
			{
			  // 2D code; z-component updated using cell-centred scheme
			  arrEM(i,j,k,EX_LOCAL+d) = 0.5*(arr(i-iOffset,j-jOffset,k-kOffset,EX+d)+arr(i,j,k,EX+d));
			}
		    }
		}      
	    }
	}
      // We need to compute boundary conditions again after each update
      S_EM_source[0].FillBoundary(geom.periodicity());
      S_EM_source[1].FillBoundary(geom.periodicity());
     
      // added by 2020D 
      // Fill non-periodic physical boundaries                          
      FillDomainBoundary(S_EM_source[0], geom, bc_EM);    
      FillDomainBoundary(S_EM_source[1], geom, bc_EM);

      MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
      MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);  
      MultiFab::Copy(S_EM_X_int, S_EM_source[0], 0, 0, 6, 0);
      FillPatch(*this, S_EM_source[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
      MultiFab::Copy(S_EM_Y_int, S_EM_source[1], 0, 0, 6, 0);
      FillPatch(*this, S_EM_source[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
#endif
    }
  
  (this->*RKWithChosenUpdateOrder)(S_dest, S_source, fluxes, S_EM_dest, S_EM_source, fluxesEM, dx, time, dt);
  implicitYeeMaxwellSolver(S_EM_dest,S_EM_source,S_dest,S_source,dx,dt,time);

  sourceUpdate(S_dest, fluxes, dx, 0.5*dt);

  if (sourceMethod=="IM")
    {
      // Update face-centred EM fields
      for (int d = 0; d < amrex::SpaceDim ; d++)   
	{

	  const int iOffset = ( d == 0 ? 1 : 0);
	  const int jOffset = ( d == 1 ? 1 : 0);
	  const int kOffset = ( d == 2 ? 1 : 0);
    
	  // Loop over all the patches at this level
	  for (MFIter mfi(S_EM_dest[d], true); mfi.isValid(); ++mfi)
	    {
	      const Box& bx = mfi.tilebox();

	      const Dim3 lo = lbound(bx);
	      const Dim3 hi = ubound(bx);

	      // Indexable arrays for the data, and the directional flux
	      // Based on the corner-centred definition of the flux array, the
	      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	      //const auto& arr = S_dest.array(mfi);
	      const auto& arrEM = S_EM_dest[d].array(mfi);
	      const auto& arr = S_dest.array(mfi);

	      const Dim3 hiDomain = ubound(geom.Domain());
      
	      for(int k = lo.z; k <= hi.z; k++)
		{
		  for(int j = lo.y; j <= hi.y; j++)
		    {
		      for(int i = lo.x; i <= hi.x; i++)
			{
			  // 2D code; z-component updated using cell-centred scheme
			  arrEM(i,j,k,EX_LOCAL+d) = 0.5*(arr(i-iOffset,j-jOffset,k-kOffset,EX+d)+arr(i,j,k,EX+d));
			}
		    }
		}      
	    }
	}
      // We need to compute boundary conditions again after each update
      S_EM_dest[0].FillBoundary(geom.periodicity());
      S_EM_dest[1].FillBoundary(geom.periodicity());
     
      // added by 2020D 
      // Fill non-periodic physical boundaries                          
      FillDomainBoundary(S_EM_dest[0], geom, bc_EM);    
      FillDomainBoundary(S_EM_dest[1], geom, bc_EM);

      MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
      MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);  
      MultiFab::Copy(S_EM_X_int, S_EM_dest[0], 0, 0, 6, 0);
      FillPatch(*this, S_EM_dest[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
      MultiFab::Copy(S_EM_Y_int, S_EM_dest[1], 0, 0, 6, 0);
      FillPatch(*this, S_EM_dest[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
#endif

      //if ((parent->levelSteps(0))%1==0)// && parent->levelSteps(0)!=0)
      {
	//amrex::Print() << parent->levelSteps(0) << std::endl;
	//Projection(S_dest,S_EM_dest,dx,time+dt);      
      }
    }
}
