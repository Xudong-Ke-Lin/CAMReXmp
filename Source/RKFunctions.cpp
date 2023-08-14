#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::RK1(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, const Real* dx, Real time, Real dt)
{

  // notice that RK1 might be unstable with higher spatial reconstruction
  amrex::Abort("Add a function like MUSCLHancokFluidSolverTVD");
  
  S_dest.FillBoundary(geom.periodicity());
  FillDomainBoundary(S_dest, geom, bc);  
  S_EM_dest[0].FillBoundary(geom.periodicity());
  FillDomainBoundary(S_EM_dest[0], geom, bc_EM);
#if (AMREX_SPACEDIM >= 2)
  S_EM_dest[1].FillBoundary(geom.periodicity());
  FillDomainBoundary(S_EM_dest[1], geom, bc_EM);  
#endif    

}

void CAMReXmp::RK2(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, const Real* dx, Real time, Real dt)
{
  /*
  //plasmaSolverTVD(S_dest,S_source,fluxes,dx,dt); 
  //return;
  MultiFab S1(grids, dmap, NUM_STATE, NUM_GROW);
  plasmaSolverTVD(S1,S_source,fluxes,dx,dt);  
  
  MultiFab S2(grids, dmap, NUM_STATE, NUM_GROW);
  plasmaSolverTVD(S2,S1,fluxes,dx,dt);

  linearCombination(S_dest, S_source, 1.0/2.0, S2, 1.0/2.0, 0, NUM_STATE);
  
  S_dest.FillBoundary(geom.periodicity());
  FillDomainBoundary(S_dest, geom, bc);  
  */

  MultiFab S1(grids, dmap, NUM_STATE, NUM_GROW);
  Array<MultiFab,AMREX_SPACEDIM> S_EM1;
  S_EM1[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);  
#if (AMREX_SPACEDIM >= 2)
  S_EM1[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
#endif

  (this->*fluidSolverWithChosenOrder)(S1,S_source,fluxes,dx,dt);  
  (this->*MaxwellSolverWithChosenOrder)(S_EM1,S_EM_source,fluxesEM,S1,S_source,fluxes,dx,dt);

  MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);  
  MultiFab::Copy(S_EM_X_int, S_EM1[0], 0, 0, 6, 0);
  FillPatch(*this, S_EM1[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_int, S_EM1[1], 0, 0, 6, 0);
  FillPatch(*this, S_EM1[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
#endif
  
  MultiFab S2(grids, dmap, NUM_STATE, NUM_GROW);
  Array<MultiFab,AMREX_SPACEDIM> S_EM2;
  S_EM2[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
#if (AMREX_SPACEDIM >= 2)
  S_EM2[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
#endif

  (this->*fluidSolverWithChosenOrder)(S2,S1,fluxes,dx,dt);
  (this->*MaxwellSolverWithChosenOrder)(S_EM2,S_EM1,fluxesEM,S2,S1,fluxes,dx,dt);

  linearCombination(S_dest, S_source, 1.0/2.0, S2, 1.0/2.0, 0, NUM_STATE);
  linearCombination(S_EM_dest[0], S_EM_source[0], 1.0/2.0, S_EM2[0], 1.0/2.0, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  linearCombination(S_EM_dest[1], S_EM_source[1], 1.0/2.0, S_EM2[1], 1.0/2.0, 0, 6);
#endif
  
  S_dest.FillBoundary(geom.periodicity());
  FillDomainBoundary(S_dest, geom, bc);  
  S_EM_dest[0].FillBoundary(geom.periodicity());
  FillDomainBoundary(S_EM_dest[0], geom, bc_EM);
#if (AMREX_SPACEDIM >= 2)
  S_EM_dest[1].FillBoundary(geom.periodicity());
  FillDomainBoundary(S_EM_dest[1], geom, bc_EM);  
#endif    

}
// RK2 with subcycling
// void CAMReXmp::RK2(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, const Real* dx, Real time, Real dt)
// {
//   MultiFab S1(grids, dmap, NUM_STATE+2, NUM_GROW);
//   Array<MultiFab,AMREX_SPACEDIM> S_EM1;
//   S_EM1[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);  
// #if (AMREX_SPACEDIM >= 2)
//   S_EM1[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
// #endif

//   MultiFab S_sourceCycle(grids, dmap, NUM_STATE+2, NUM_GROW);
//   MultiFab::Copy(S_sourceCycle, S_source, 0, 0, NUM_STATE+2, NUM_GROW);
//   Array<MultiFab,AMREX_SPACEDIM> S_EM_sourceCycle;
//   S_EM_sourceCycle[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);  
// #if (AMREX_SPACEDIM >= 2)
//   S_EM_sourceCycle[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
// #endif
//   MultiFab::Copy(S_EM_sourceCycle[0], S_EM_source[0], 0, 0, 6, NUM_GROW);
//   MultiFab::Copy(S_EM_sourceCycle[1], S_EM_source[1], 0, 0, 6, NUM_GROW);

//   (this->*fluidSolverWithChosenOrder)(S1,S_source,fluxes,dx,dt);
//   Real cflEM=cfl/1.0;
//   Real dtEM = cflEM*std::min(dx[0],dx[1])/c;
//   Real dt_current = dtEM;
//   do{
//     amrex::Print() << "Updating hyperbolic part, dt=" << dtEM << ", current dt=" << dt_current << " and final dt=" << dt << std::endl;
//     (this->*MaxwellSolverWithChosenOrder)(S_EM1,S_EM_sourceCycle,fluxesEM,S1,S_sourceCycle,fluxes,dx,dtEM);

//     MultiFab::Copy(S_sourceCycle, S1, BX, BX, 6, NUM_GROW);
    
//     MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
//     MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);  
//     MultiFab::Copy(S_EM_X_int, S_EM1[0], 0, 0, 6, 0);
//     MultiFab::Copy(S_EM_Y_int, S_EM1[1], 0, 0, 6, 0);
//     FillPatch(*this, S_EM_sourceCycle[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
//     FillPatch(*this, S_EM_sourceCycle[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
    
//     if (dt_current+dtEM > dt)
//       {
// 	// if at final time step, stop
// 	if (dt_current==dt)
// 	  break;
// 	// at the final time step, use the remaining dt
// 	dtEM = std::abs(dt - dt_current);      	
//       }
//     dt_current += dtEM;
//   } while (dt_current <= dt);

//   MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
//   MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);  
//   MultiFab::Copy(S_EM_X_int, S_EM1[0], 0, 0, 6, 0);
//   MultiFab::Copy(S_EM_Y_int, S_EM1[1], 0, 0, 6, 0);
//   FillPatch(*this, S_EM1[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
//   FillPatch(*this, S_EM1[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
  
//   MultiFab S2(grids, dmap, NUM_STATE+2, NUM_GROW);
//   Array<MultiFab,AMREX_SPACEDIM> S_EM2;
//   S_EM2[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
// #if (AMREX_SPACEDIM >= 2)
//   S_EM2[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
// #endif
//   /*    
//     // Compute y-components slopes in x-direction
//   for (MFIter mfi(S_EM1[1], true); mfi.isValid(); ++mfi)
//     {
//       const Box& bx = mfi.tilebox();
      
//       const Dim3 lo = lbound(bx);
//       const Dim3 hi = ubound(bx);

//       const Dim3 hiDomain = ubound(geom.Domain());

//       const Array4<Real> arr = S_EM1[1].array(mfi);
      
//       for(int k = lo.z; k <= hi.z; k++)
// 	{
// 	  for(int j = lo.y-1; j <= hi.y+1; j++)
// 	    {
// 	      for(int i = lo.x-1; i <= hi.x+1; i++)
// 		{
// 		  if (j>=hiDomain.y-1 && i==188)
// 		    std::cout << "before " << j << " " << arr(i,j-1,k,1) << " " << arr(i,j,k,1) << " " << arr(i,j+1,k,1) << std::endl;
// 		}
// 	    }
// 	}      
//     }
//   */
//   (this->*fluidSolverWithChosenOrder)(S2,S1,fluxes,dx,dt);
//   dtEM = cflEM*std::min(dx[0],dx[1])/c;
//   dt_current = dtEM;
//   do{
//     amrex::Print() << "Updating hyperbolic part, dt=" << dtEM << ", current dt=" << dt_current << " and final dt=" << dt << std::endl;
//     (this->*MaxwellSolverWithChosenOrder)(S_EM2,S_EM1,fluxesEM,S2,S1,fluxes,dx,dtEM);

//     MultiFab::Copy(S1, S2, BX, BX, 6, NUM_GROW);
    
//     MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
//     MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);  
//     MultiFab::Copy(S_EM_X_int, S_EM2[0], 0, 0, 6, 0);
//     MultiFab::Copy(S_EM_Y_int, S_EM2[1], 0, 0, 6, 0);
//     FillPatch(*this, S_EM1[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
//     FillPatch(*this, S_EM1[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
//     /*
//     // Compute y-components slopes in x-direction
//   for (MFIter mfi(S_EM1[1], true); mfi.isValid(); ++mfi)
//     {
//       const Box& bx = mfi.tilebox();
      
//       const Dim3 lo = lbound(bx);
//       const Dim3 hi = ubound(bx);

//       const Dim3 hiDomain = ubound(geom.Domain());

//       const Array4<Real> arr = S_EM1[1].array(mfi);
      
//       for(int k = lo.z; k <= hi.z; k++)
// 	{
// 	  for(int j = lo.y-1; j <= hi.y+1; j++)
// 	    {
// 	      for(int i = lo.x-1; i <= hi.x+1; i++)
// 		{
// 		  if (j>=hiDomain.y-1 && i==188)
// 		    std::cout << "after " << j << " " << arr(i,j-1,k,1) << " " << arr(i,j,k,1) << " " << arr(i,j+1,k,1) << std::endl;
// 		}
// 	    }
// 	}      
//     }
//     */
//   if (dt_current+dtEM > dt)
//       {
// 	// if at final time step, stop
// 	if (dt_current==dt)
// 	  break;
// 	// at the final time step, use the remaining dt
// 	dtEM = std::abs(dt - dt_current);      	
//       }
//     dt_current += dtEM;
//   } while (dt_current <= dt);  

//   linearCombination(S_dest, S_source, 1.0/2.0, S2, 1.0/2.0, 0, NUM_STATE+2);
//   linearCombination(S_EM_dest[0], S_EM_source[0], 1.0/2.0, S_EM2[0], 1.0/2.0, 0, 6);
// #if (AMREX_SPACEDIM >= 2)
//   linearCombination(S_EM_dest[1], S_EM_source[1], 1.0/2.0, S_EM2[1], 1.0/2.0, 0, 6);
// #endif
  
//   S_dest.FillBoundary(geom.periodicity());
//   FillDomainBoundary(S_dest, geom, bc);  
//   S_EM_dest[0].FillBoundary(geom.periodicity());
//   FillDomainBoundary(S_EM_dest[0], geom, bc_EM);
// #if (AMREX_SPACEDIM >= 2)
//   S_EM_dest[1].FillBoundary(geom.periodicity());
//   FillDomainBoundary(S_EM_dest[1], geom, bc_EM);  
// #endif    
// }
void CAMReXmp::RK3(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, const Real* dx, Real time, Real dt)
{
  MultiFab S1(grids, dmap, NUM_STATE+2, NUM_GROW);
  Array<MultiFab,AMREX_SPACEDIM> S_EM1;
  S_EM1[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);  
#if (AMREX_SPACEDIM >= 2)
  S_EM1[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
#endif
  
  (this->*fluidSolverWithChosenOrder)(S1,S_source,fluxes,dx,dt);
  (this->*MaxwellSolverWithChosenOrder)(S_EM1,S_EM_source,fluxesEM,S1,S_source,fluxes,dx,dt);
  
  MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);  
  MultiFab::Copy(S_EM_X_int, S_EM1[0], 0, 0, 6, 0);
  FillPatch(*this, S_EM1[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_int, S_EM1[1], 0, 0, 6, 0);  
  FillPatch(*this, S_EM1[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
#endif
  
  MultiFab S2(grids, dmap, NUM_STATE+2, NUM_GROW);
  Array<MultiFab,AMREX_SPACEDIM> S_EM2;
  S_EM2[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
#if (AMREX_SPACEDIM >= 2)
  S_EM2[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
#endif

  (this->*fluidSolverWithChosenOrder)(S2,S1,fluxes,dx,dt);
  (this->*MaxwellSolverWithChosenOrder)(S_EM2,S_EM1,fluxesEM,S2,S1,fluxes,dx,dt);
  
  linearCombination(S1, S_source, 0.75, S2, 0.25, 0, NUM_STATE+2);
  linearCombination(S_EM1[0], S_EM_source[0], 0.75, S_EM2[0], 0.25, 0, 6);
#if (AMREX_SPACEDIM >= 2)  
  linearCombination(S_EM1[1], S_EM_source[1], 0.75, S_EM2[1], 0.25, 0, 6);
#endif
  
  // boundary conditions
  S1.FillBoundary(geom.periodicity());
  FillDomainBoundary(S1, geom, bc);  
  MultiFab::Copy(S_EM_X_int, S_EM1[0], 0, 0, 6, 0);
  FillPatch(*this, S_EM1[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_int, S_EM1[1], 0, 0, 6, 0);  
  FillPatch(*this, S_EM1[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
#endif
  
  (this->*fluidSolverWithChosenOrder)(S2,S1,fluxes,dx,dt);
  (this->*MaxwellSolverWithChosenOrder)(S_EM2,S_EM1,fluxesEM,S2,S1,fluxes,dx,dt);

  linearCombination(S_dest, S_source, 1.0/3.0, S2, 2.0/3.0, 0, NUM_STATE+2);
  linearCombination(S_EM_dest[0], S_EM_source[0], 1.0/3.0, S_EM2[0], 2.0/3.0, 0, 6);
#if (AMREX_SPACEDIM >= 2)  
  linearCombination(S_EM_dest[1], S_EM_source[1], 1.0/3.0, S_EM2[1], 2.0/3.0, 0, 6);
#endif
  
  S_dest.FillBoundary(geom.periodicity());
  FillDomainBoundary(S_dest, geom, bc);  
  S_EM_dest[0].FillBoundary(geom.periodicity());
  FillDomainBoundary(S_EM_dest[0], geom, bc_EM);
#if (AMREX_SPACEDIM >= 2)
  S_EM_dest[1].FillBoundary(geom.periodicity());
  FillDomainBoundary(S_EM_dest[1], geom, bc_EM);  
#endif    
}
void CAMReXmp::linearCombination(MultiFab& S_new, MultiFab& S1, Real a1, MultiFab& S2, Real a2, int idxStart, int num_var)
{
  for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      
      Array4<Real> arr = S_new.array(mfi);
      Array4<Real> arr1 = S1.array(mfi);
      Array4<Real> arr2 = S2.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  for(int n = idxStart; n < idxStart+num_var; n++)
		    arr(i,j,k,n) = a1*arr1(i,j,k,n) + a2*arr2(i,j,k,n);
		    //arr(i,j,k,n) = arr1(i,j,k,n);
		}
	    }
	}      
    }
  // We need to compute boundary conditions again after each update
  //S_new.FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries  
  //FillDomainBoundary(S_new, geom, bc);

}
