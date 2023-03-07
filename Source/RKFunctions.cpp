#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

//void CAMReXmp::RK1(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt,
//		   std::function<void (MultiFab&, MultiFab (&)[AMREX_SPACEDIM], const Real*, Real)> func)
void CAMReXmp::RK1(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt, functionPointer func, int idxStart, int num_var)
{
  (this->*func)(Sborder, fluxes, dx, dt);
}

void CAMReXmp::RK2(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt, functionPointer func, int idxStart, int num_var)
{
  // Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes1[amrex::SpaceDim];
  MultiFab fluxes2[amrex::SpaceDim];
  
  // Define the appropriate size for the flux MultiFab.
  for (int j = 0; j < amrex::SpaceDim; j++)
  {
    BoxArray ba = Sborder.boxArray();
    ba.surroundingNodes(j);
    fluxes1[j].define(ba, dmap, NUM_STATE, 0);
    fluxes2[j].define(ba, dmap, NUM_STATE, 0); 
  }

  MultiFab S0(grids, dmap, NUM_STATE, NUM_GROW);
  MultiFab S1(grids, dmap, NUM_STATE, NUM_GROW);
  MultiFab::Copy(S0, Sborder, idxStart, idxStart, num_var, NUM_GROW);
  MultiFab::Copy(S1, Sborder, idxStart, idxStart, num_var, NUM_GROW);  
  (this->*func)(S1, fluxes1, dx, dt);
  
  MultiFab S2(grids, dmap, NUM_STATE, NUM_GROW);
  MultiFab::Copy(S2, S1, idxStart, idxStart, num_var, NUM_GROW);
  //(this->*func)(S2, fluxes2, dx, dt);
  (this->*func)(S2, fluxes2, dx, dt);

  linearCombination(Sborder, S0, 0.5, S2, 0.5, idxStart, num_var);
  //MultiFab::Copy(Sborder, S1, idxStart, idxStart, num_var, NUM_GROW);
}
void CAMReXmp::hyperbolicMaxwellRK(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  //(this->*RKWithChosenUpdateOrder)(Sborder, fluxes, dx, dt, &CAMReXmp::hyperbolicMaxwellSolver, BX, NUM_STATE_MAXWELL);
  RK1(Sborder, fluxes, dx, dt, &CAMReXmp::hyperbolicMaxwellSolver, BX, NUM_STATE_MAXWELL);
}
void CAMReXmp::implicitMaxwellRK(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  implicitMaxwellSolver(Sborder, dx, dt);
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
