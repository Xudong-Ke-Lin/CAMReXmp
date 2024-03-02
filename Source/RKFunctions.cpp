#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::RK2(const Real* dx, Real dt, Real time)
{

  // get multifabs references
  MultiFab& S_new = get_new_data(Phi_Type);
#if (AMREX_SPACEDIM >= 2)  
  MultiFab& S_EM_X_new = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_new = get_new_data(EM_Y_Type);  
#endif
  
  // input states
  MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
#if (AMREX_SPACEDIM >= 2)  
  Array<MultiFab,AMREX_SPACEDIM> S_EM_input;
  S_EM_input[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW); 
  S_EM_input[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
#endif
  
  // intermediate states in RK2
  MultiFab S1(grids, dmap, NUM_STATE, NUM_GROW);
#if (AMREX_SPACEDIM >= 2)  
  Array<MultiFab,AMREX_SPACEDIM> S_EM1;
  S_EM1[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);  
  S_EM1[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
#endif  

  //////////////////////////////////////////////////////////////////////////////////////////
  // fluid or plasma solver
  //////////////////////////////////////////////////////////////////////////////////////////
  // fill initial input data
  if (fluidOrder!=0)
    {
      FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
      if (MaxwellDivMethod=="HDC")
	plasmaSolverTVD(S_input,dx,dt);
      else
	fluidSolverTVD(S_input,dx,dt);
      // fill intermediate state
      FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    }
  //////////////////////////////////////////////////////////////////////////////////////////
  // Maxwell solver for FVTD with explicit time stepping
  //////////////////////////////////////////////////////////////////////////////////////////
#if (AMREX_SPACEDIM >= 2)  
  if (MaxwellOrder!=0 && MaxwellDivMethod=="FVTD" && MaxwellTimeMethod=="EX")
    {
      // fill initial input data
      FillPatch(*this, S_EM_input[0], NUM_GROW, time, EM_X_Type, 0, 6);
      FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);
      //MaxwellSolverFVTDTVD(S_EM_input,S_input,dx,dt);
      (this->*MaxwellSolverWithChosenSpaceOrder)(S_EM_input,S_input,dx,dt);
      // fill intermediate states
      FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
      FillPatch(*this, S_EM1[0], NUM_GROW, time, EM_X_Type, 0, 6);
      FillPatch(*this, S_EM1[1], NUM_GROW, time, EM_Y_Type, 0, 6);
    }
#endif
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // fluid or plasma solver
  //////////////////////////////////////////////////////////////////////////////////////////
  if (fluidOrder!=0)
    {
      if (MaxwellDivMethod=="HDC")
	{
	  plasmaSolverTVD(S1,dx,dt);
	  linearCombination(S_new, S_new, 1.0/2.0, S_input, 1.0/2.0, 0, NUM_STATE);
	}
      else
	{
	  fluidSolverTVD(S1,dx,dt);
	  // only needs to update the fluid variables
	  linearCombination(S_new, S_new, 1.0/2.0, S_input, 1.0/2.0, 0, NUM_STATE_FLUID);
	}
    }
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Maxwell solver for FVTD with explicit time stepping
  //////////////////////////////////////////////////////////////////////////////////////////
#if (AMREX_SPACEDIM >= 2)
  if (MaxwellOrder!=0 && MaxwellDivMethod=="FVTD" && MaxwellTimeMethod=="EX")
    {
      //MaxwellSolverFVTDTVD(S_EM1,S1,dx,dt);
      (this->*MaxwellSolverWithChosenSpaceOrder)(S_EM1,S1,dx,dt);
      linearCombination(S_new, S_new, 1.0/2.0, S_input, 1.0/2.0, BX, 6);
      linearCombination(S_EM_X_new, S_EM_X_new, 1.0/2.0, S_EM_input[0], 1.0/2.0, 0, 6);
      linearCombination(S_EM_Y_new, S_EM_Y_new, 1.0/2.0, S_EM_input[1], 1.0/2.0, 0, 6);
    }
#endif  
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Maxwell solver for FVTD with explicit time stepping with sub-cycling
  //////////////////////////////////////////////////////////////////////////////////////////
#if (AMREX_SPACEDIM >= 2)  
  if (MaxwellOrder!=0 && MaxwellDivMethod=="FVTD" && MaxwellTimeMethod=="EXsubcycling")
    {
      Real dtEM = cfl*std::min(dx[0],dx[1])/c;
      //Real dtEM = 0.4*std::min(dx[0],dx[1])/c;
      Real dt_current = dtEM; 
  
      do{
	amrex::Print() << "Updating hyperbolic part, dt=" << dtEM << ", current dt=" << dt_current << " and final dt=" << dt << std::endl;

	// fill initial input data
	// note that it is using the updated fluid states for the entire loop
	FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
	FillPatch(*this, S_EM_input[0], NUM_GROW, time, EM_X_Type, 0, 6);
	FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);

	//MaxwellSolverFVTDTVD(S_EM_input,S_input,dx,dtEM);
	(this->*MaxwellSolverWithChosenSpaceOrder)(S_EM_input,S_input,dx,dtEM);
	// fill intermediate states 
	FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
	FillPatch(*this, S_EM1[0], NUM_GROW, time, EM_X_Type, 0, 6);
	FillPatch(*this, S_EM1[1], NUM_GROW, time, EM_Y_Type, 0, 6);

	//MaxwellSolverFVTDTVD(S_EM1,S1,dx,dtEM);
	(this->*MaxwellSolverWithChosenSpaceOrder)(S_EM1,S1,dx,dtEM);
	linearCombination(S_new, S_new, 1.0/2.0, S_input, 1.0/2.0, BX, 6);
	linearCombination(S_EM_X_new, S_EM_X_new, 1.0/2.0, S_EM_input[0], 1.0/2.0, 0, 6);
	linearCombination(S_EM_Y_new, S_EM_Y_new, 1.0/2.0, S_EM_input[1], 1.0/2.0, 0, 6);
    
	if (dt_current+dtEM > dt)
	  {
	    // if at final time step, stop
	    if (dt_current==dt)
	      break;
	    // at the final time step, use the remaining dt
	    dtEM = std::abs(dt - dt_current);
	    if (dtEM<1e-14)
	      break;
	  }
	dt_current += dtEM;
      } while (dt_current <= dt);
    }
#endif            
}
void CAMReXmp::RK3(const Real* dx, Real dt, Real time)
{

  // get multifabs references
  MultiFab& S_new = get_new_data(Phi_Type);
  MultiFab& S_EM_X_new = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_new = get_new_data(EM_Y_Type);  
  
  // input states
  MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
#if (AMREX_SPACEDIM >= 2)   
  Array<MultiFab,AMREX_SPACEDIM> S_EM_input;
  S_EM_input[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  S_EM_input[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
#endif
  
  // intermediate states in RK2
  MultiFab S1(grids, dmap, NUM_STATE, NUM_GROW);
  MultiFab S2(grids, dmap, NUM_STATE, NUM_GROW);
#if (AMREX_SPACEDIM >= 2)  
  Array<MultiFab,AMREX_SPACEDIM> S_EM1, S_EM2;
  S_EM1[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);  
  S_EM2[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);  
  S_EM1[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  S_EM2[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
#endif  

  //////////////////////////////////////////////////////////////////////////////////////////
  // fluid or plasma solver
  //////////////////////////////////////////////////////////////////////////////////////////
  // fill initial input data
  if (fluidOrder!=0)
    {
      FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
      if (MaxwellDivMethod=="HDC")
	//plasmaSolverTVD(S_input,dx,dt);
	amrex::Abort("Currently does not support WENO for HDC");
      else
	(this->*fluidSolverWithChosenSpaceOrder)(S_input,dx,dt);
      // fill intermediate state
      FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    }
  //////////////////////////////////////////////////////////////////////////////////////////
  // Maxwell solver for FVTD with explicit time stepping
  //////////////////////////////////////////////////////////////////////////////////////////
#if (AMREX_SPACEDIM >= 2)  
  if (MaxwellOrder!=0 && MaxwellDivMethod=="FVTD" && MaxwellTimeMethod=="EX")
    {
      // fill initial input data
      FillPatch(*this, S_EM_input[0], NUM_GROW, time, EM_X_Type, 0, 6);
      FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);
      (this->*MaxwellSolverWithChosenSpaceOrder)(S_EM_input,S_input,dx,dt);
      // fill intermediate states
      FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
      FillPatch(*this, S_EM1[0], NUM_GROW, time, EM_X_Type, 0, 6);
      FillPatch(*this, S_EM1[1], NUM_GROW, time, EM_Y_Type, 0, 6);
    }
#endif
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // fluid or plasma solver
  //////////////////////////////////////////////////////////////////////////////////////////
  if (fluidOrder!=0)
    {
      if (MaxwellDivMethod=="HDC")
	{
	  //plasmaSolverTVD(S1,dx,dt);
	}
      else
	{
	  (this->*fluidSolverWithChosenSpaceOrder)(S1,dx,dt);
	  // only needs to update the fluid variables
	  FillPatch(*this, S2, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
	  linearCombination(S_new, S_input, 0.75, S2, 0.25, 0, NUM_STATE_FLUID);
	  FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
	}
    }
  //////////////////////////////////////////////////////////////////////////////////////////
  // Maxwell solver for FVTD with explicit time stepping
  //////////////////////////////////////////////////////////////////////////////////////////
#if (AMREX_SPACEDIM >= 2)  
  if (MaxwellOrder!=0 && MaxwellDivMethod=="FVTD" && MaxwellTimeMethod=="EX")
    {
      (this->*MaxwellSolverWithChosenSpaceOrder)(S_EM1,S1,dx,dt);

      // fill intermediate states 
      FillPatch(*this, S2, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
      FillPatch(*this, S_EM2[0], NUM_GROW, time, EM_X_Type, 0, 6);
      FillPatch(*this, S_EM2[1], NUM_GROW, time, EM_Y_Type, 0, 6);
	
      linearCombination(S_new, S_input, 0.75, S2, 0.25, BX, 6);
      linearCombination(S_EM_X_new, S_EM_input[0], 0.75, S_EM2[0], 0.25, 0, 6);
      linearCombination(S_EM_Y_new, S_EM_input[1], 0.75, S_EM2[1], 0.25, 0, 6);

      // fill intermediate states 
      FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
      FillPatch(*this, S_EM1[0], NUM_GROW, time, EM_X_Type, 0, 6);
      FillPatch(*this, S_EM1[1], NUM_GROW, time, EM_Y_Type, 0, 6);
    }
#endif                  

  //////////////////////////////////////////////////////////////////////////////////////////
  // fluid or plasma solver
  //////////////////////////////////////////////////////////////////////////////////////////
  if (fluidOrder!=0)
    {
      if (MaxwellDivMethod=="HDC")
	{
	  //plasmaSolverTVD(S1,dx,dt);
	}
      else
	{
	  (this->*fluidSolverWithChosenSpaceOrder)(S1,dx,dt);
	  // only needs to update the fluid variables
	  linearCombination(S_new, S_new, 2.0/3.0, S_input, 1.0/3.0, 0, NUM_STATE_FLUID);
	}
    }
  //////////////////////////////////////////////////////////////////////////////////////////
  // Maxwell solver for FVTD with explicit time stepping
  //////////////////////////////////////////////////////////////////////////////////////////
#if (AMREX_SPACEDIM >= 2)        
  if (MaxwellOrder!=0 && MaxwellDivMethod=="FVTD" && MaxwellTimeMethod=="EX")
    {
      (this->*MaxwellSolverWithChosenSpaceOrder)(S_EM1,S1,dx,dt);

      linearCombination(S_new, S_new, 2.0/3.0, S_input, 1.0/3.0, BX, 6);
      linearCombination(S_EM_X_new, S_EM_X_new, 2.0/3.0, S_EM_input[0], 1.0/3.0, 0, 6);
      linearCombination(S_EM_Y_new, S_EM_Y_new, 2.0/3.0, S_EM_input[1], 1.0/3.0, 0, 6);
    }
#endif  
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // Maxwell solver for FVTD with explicit time stepping with sub-cycling
  //////////////////////////////////////////////////////////////////////////////////////////
#if (AMREX_SPACEDIM >= 2)  
  if (MaxwellOrder!=0 && MaxwellDivMethod=="FVTD" && MaxwellTimeMethod=="EXsubcycling")
    {

      Real dtEM = cfl*std::min(dx[0],dx[1])/c;
      //Real dtEM = 0.4*std::min(dx[0],dx[1])/c;
      Real dt_current = dtEM; 
  
      do{
	amrex::Print() << "Updating hyperbolic part, dt=" << dtEM << ", current dt=" << dt_current << " and final dt=" << dt << std::endl;

	// fill initial input data
	// note that it is using the updated fluid states for the entire loop
	FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);	
	FillPatch(*this, S_EM_input[0], NUM_GROW, time, EM_X_Type, 0, 6);
	FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);
	
	(this->*MaxwellSolverWithChosenSpaceOrder)(S_EM_input,S_input,dx,dtEM);
	// fill intermediate states 
	FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
	FillPatch(*this, S_EM1[0], NUM_GROW, time, EM_X_Type, 0, 6);
	FillPatch(*this, S_EM1[1], NUM_GROW, time, EM_Y_Type, 0, 6);

	(this->*MaxwellSolverWithChosenSpaceOrder)(S_EM1,S1,dx,dtEM);
	// fill intermediate states 
	FillPatch(*this, S2, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
	FillPatch(*this, S_EM2[0], NUM_GROW, time, EM_X_Type, 0, 6);
	FillPatch(*this, S_EM2[1], NUM_GROW, time, EM_Y_Type, 0, 6);
	
	linearCombination(S_new, S_input, 0.75, S2, 0.25, BX, 6);
	linearCombination(S_EM_X_new, S_EM_input[0], 0.75, S_EM2[0], 0.25, 0, 6);
	linearCombination(S_EM_Y_new, S_EM_input[1], 0.75, S_EM2[1], 0.25, 0, 6);
	// fill intermediate states 
	FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
	FillPatch(*this, S_EM1[0], NUM_GROW, time, EM_X_Type, 0, 6);
	FillPatch(*this, S_EM1[1], NUM_GROW, time, EM_Y_Type, 0, 6);

	(this->*MaxwellSolverWithChosenSpaceOrder)(S_EM1,S1,dx,dtEM);

	linearCombination(S_new, S_new, 2.0/3.0, S_input, 1.0/3.0, BX, 6);
	linearCombination(S_EM_X_new, S_EM_X_new, 2.0/3.0, S_EM_input[0], 1.0/3.0, 0, 6);
	linearCombination(S_EM_Y_new, S_EM_Y_new, 2.0/3.0, S_EM_input[1], 1.0/3.0, 0, 6);
    
	if (dt_current+dtEM > dt)
	  {
	    // if at final time step, stop
	    if (dt_current==dt)
	      break;
	    // at the final time step, use the remaining dt
	    dtEM = std::abs(dt - dt_current);
	    if (dtEM<1e-14)
	      break;
	  }
	dt_current += dtEM;
      } while (dt_current <= dt);
    }
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
		}
	    }
	}      
    }
  // We need to compute boundary conditions again after each update
  //S_new.FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries  
  //FillDomainBoundary(S_new, geom, bc);

}
void CAMReXmp::RK2(const Real* dx, Real dt, Real time, int start, int len,
		   std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver)
{
  // get multifabs references
  MultiFab& S_new = get_new_data(Phi_Type);

  // input states
  MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  // intermediate states in RK2
  MultiFab S1(grids, dmap, NUM_STATE, NUM_GROW);

  generalSolverTVD(S_input,dx,dt,start,len,solver);
  FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  generalSolverTVD(S1,dx,dt,start,len,solver);

  MultiFab::LinComb(S_new, 0.5, S_new, start, 0.5, S_input, start, start, len, 0);

}
