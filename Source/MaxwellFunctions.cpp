#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

// Implicit Maxwell Solver
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>

//#include <AMReX_Hypre.H>
#include <AMReX_HypreSolver.H>

// AMReX linear solvers in nodal variables
#include <AMReX_MLNodeABecLaplacian.H>

using namespace amrex;

/*void CAMReXmp::MaxwellSolver(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  (this->*MaxwellSolverWithChosenMethod)(Sborder,fluxes,dx,dt);
}
*/
//template<int size>
void CAMReXmp::hyperbolicMaxwellSolver(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  /*Vector<Real> dim;
  if (geom.Coord()==0)
    dim = {0,1};
  else if (geom.Coord()==1)
  dim = {0,2};*/
  /*ParmParse pp;
  std::string coor;
  pp.query("coor",coor);
  Vector<Real> dim;
  if (coor=="xz")
    dim = {0,2};
  else
    dim = {0,1};
  */

  MFIter::allowMultipleMFIters(true);

  int from, to, step;
  if (geom.Coord()==0){
    from = 0;
    to = amrex::SpaceDim;
    step = 1;
  } else {
    from = amrex::SpaceDim-1;
    to = -1;
    step = -1;
  }
  
  for (int d = from; d != to; d += step)
    //for (int d = amrex::SpaceDim-1; d >= 0 ; d--)
    //for (int d = 0; d < amrex::SpaceDim ; d++)
    {
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
      
      for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);
	  
	  const auto& arr = Sborder.array(mfi);
	  const auto& fluxArr = fluxes[d].array(mfi);
	  
	  for(int k = lo.z; k <= hi.z+kOffset; k++)
	    {
	      for(int j = lo.y; j <= hi.y+jOffset; j++)
		{
		  for(int i = lo.x; i <= hi.x+iOffset; i++)
		    {
		      Vector<Real> flux = Maxwell_flux_Godunov(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);
		      //Vector<Real> flux = Maxwell_flux_HLLE(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);
		      //Vector<Real> flux = flux_HLLC(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, dim[d]);   
		      for(int n=NUM_STATE_FLUID; n<NUM_STATE; n++)
			//for(int n=0; n<NUM_STATE; n++)
			{
			  //amrex::Print() << flux[n] << " ";
			  fluxArr(i,j,k,n) = flux[n-NUM_STATE_FLUID];
			  //fluxArr(i,j,k,n) = flux[n];
			}	   
		    }
		}
	    }
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {	    
		      for(int n=NUM_STATE_FLUID; n<NUM_STATE; n++)
			//for(int n=0; n<NUM_STATE; n++)
			{
			  // Conservative update formula
			  arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
			}
		    }
		}
	    }	  	 	 
	}
      // We need to compute boundary conditions again after each update
      Sborder.FillBoundary(geom.periodicity());
      
      // Fill non-periodic physical boundaries  
      FillDomainBoundary(Sborder, geom, bc);  
      
    }
}
void CAMReXmp::implicitMaxwellSolverSetUp()
{
  // Set boundary conditions for MLABecLaplacian  
  setDomainBC(mlmg_lobc_X, mlmg_hibc_X, BX);
  setDomainBC(mlmg_lobc_Y, mlmg_hibc_Y, BY);
  setDomainBC(mlmg_lobc_Z, mlmg_hibc_Z, BZ);
  
  // The A coefficient (acoef) multiplies the T^(n+1) element of the
  // linear system.  For the heat equation, this is \rho(T) * c_v(T),
  // and is set by an call to a fortran function.
  computeAlpha(acoef);
  
  // The B coefficient (bcoeffs) multiply the gradient within the
  // linear solve - due to the nature of the solver, they are defined
  // at cell edges.  For the heat equation, these coefficients are the
  // thermal conductivity, kappa.
  computeBeta(bcoeffs);
  
}
//template<int size>
void CAMReXmp::implicitMaxwellSolver(MultiFab& Sborder, const Real* dx, Real dt) 
{
  MultiFab Sold(grids, dmap, NUM_STATE, NUM_GROW);
  MultiFab::Copy(Sold, Sborder, 0, 0, NUM_STATE, NUM_GROW);
  
  LPInfo info;
  info.setAgglomeration(1);
  info.setConsolidation(1);
  info.setMetricTerm(false);
  
  // Implicit solve using MLABecLaplacian class
  MLABecLaplacian mlabecX({geom}, {grids}, {dmap}, info);
  mlabecX.setMaxOrder(max_order);
  MLABecLaplacian mlabecY({geom}, {grids}, {dmap}, info);
  mlabecY.setMaxOrder(max_order);
  MLABecLaplacian mlabecZ({geom}, {grids}, {dmap}, info);
  mlabecZ.setMaxOrder(max_order);
  
  // Define the current multifab
  MultiFab J(grids,dmap,3,1);  
  
  // compute half time step update
  computeCurrentUpdate(J, Sborder, dx, 0.5*dt);
  
  // Set boundary conditions for MLABecLaplacian  
  mlabecX.setDomainBC(mlmg_lobc_X, mlmg_hibc_X);
  mlabecY.setDomainBC(mlmg_lobc_Y, mlmg_hibc_Y);
  mlabecZ.setDomainBC(mlmg_lobc_Z, mlmg_hibc_Z);
  
  MultiFab S_X(Sborder, amrex::make_alias, BX, 1);
  MultiFab S_Y(Sborder, amrex::make_alias, BY, 1);
  MultiFab S_Z(Sborder, amrex::make_alias, BZ, 1);  
  
  // Set boundary conditions for the current patch 
  mlabecX.setLevelBC(0,&S_X);
  mlabecY.setLevelBC(0,&S_Y);
  mlabecZ.setLevelBC(0,&S_Z);  
  
  // Coefficients a, b, are constant multipliers within the heat
  // equation, for all cases so far considered, these are a=1, b=dt
  // (any other coefficients can be placed within matrix multipliers).
  // Therefore we give these hard-coded values.
  Real a = 1.0;
  Real b = tau*tau*c*c*dt*dt;
  
  mlabecX.setScalars(a, b);
  mlabecX.setACoeffs(0, acoef);
  mlabecY.setScalars(a, b);
  mlabecY.setACoeffs(0, acoef);
  mlabecZ.setScalars(a, b);
  mlabecZ.setACoeffs(0, acoef);
  
  mlabecX.setBCoeffs(0,amrex::GetArrOfConstPtrs(bcoeffs));
  mlabecY.setBCoeffs(0,amrex::GetArrOfConstPtrs(bcoeffs));
  mlabecZ.setBCoeffs(0,amrex::GetArrOfConstPtrs(bcoeffs));  
  
  // Set the RHS to be multiplied appropriately (by rho * c_v, i.e.\ acoef)
  // for(MFIter RHSmfi(Rhs,true); RHSmfi.isValid(); ++RHSmfi)
  // {
  //   const Box& box = RHSmfi.tilebox();
  //   Rhs[RHSmfi].mult(acoef[RHSmfi],box,0,0,1);
  // }
  std::array<MultiFab,3> Rhs;
  Rhs[0].define(grids, dmap, 1, 0);
  Rhs[1].define(grids, dmap, 1, 0);
  Rhs[2].define(grids, dmap, 1, 0);
  MultiFab::Copy(Rhs[0], Sborder, 0, BX, 1, 0);
  MultiFab::Copy(Rhs[1], Sborder, 0, BY, 1, 0);
  MultiFab::Copy(Rhs[2], Sborder, 0, BZ, 1, 0);
  computeRhs(Rhs, J, Sborder, dx, dt);
  //computeRhsX(Rhs[0], J, Sborder, 0.0, dt, dx[0], dx[1]);
  //computeRhsY(Rhs[1], J, Sborder, 0.0, dt, dx[0], dx[1]);
  //computeRhsZ(Rhs[2], J, Sborder, 0.0, dt, dx[0], dx[1]);    

  MLMG mlmgX(mlabecX);
  MLMG mlmgY(mlabecY);
  MLMG mlmgZ(mlabecZ);
  
  mlmgX.setMaxFmgIter(max_fmg_iter);
  mlmgX.setVerbose(verbose);
  mlmgY.setMaxFmgIter(max_fmg_iter);
  mlmgY.setVerbose(verbose);
  mlmgZ.setMaxFmgIter(max_fmg_iter);
  mlmgZ.setVerbose(verbose);
  
  const Real S_X_abs = soln_tol*Rhs[0].norm0();
  const Real S_Y_abs = soln_tol*Rhs[1].norm0();
  const Real S_Z_abs = soln_tol*Rhs[2].norm0();  

  // Solve to get S^(n+1)
  mlmgX.solve({&S_X}, {&Rhs[0]}, soln_tol, S_X_abs);
  mlmgY.solve({&S_Y}, {&Rhs[1]}, soln_tol, S_Y_abs);
  mlmgZ.solve({&S_Z}, {&Rhs[2]}, soln_tol, S_Z_abs);

  // We need to compute boundary conditions again after each update
  Sborder.FillBoundary(geom.periodicity());
  // Fill non-periodic physical boundaries
  FillDomainBoundary(Sborder, geom, bc);

  // Update electric field
  implicitMaxwellSolverElecFieldUpdate(Sborder, Sold, J, dx, dt);
}
void CAMReXmp::implicitMaxwellSolverElecFieldUpdate(MultiFab& Sborder, const amrex::MultiFab& Sold, const amrex::MultiFab& J, const Real* dx, Real dt)
{
  MFIter::allowMultipleMFIters(true);
  
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

	  const auto& arr = Sborder.array(mfi);
	  const auto& arr_old = Sold.array(mfi);
	  const auto& J_nPlusHalf = J.array(mfi); 

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      Real dxBy = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,BX+(1+d)%3),
						    arr(i+iOffset,j+jOffset,k+kOffset,BX+(1+d)%3),dx[d]);
		      Real dxBz = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,BX+(2+d)%3),
						    arr(i+iOffset,j+jOffset,k+kOffset,BX+(2+d)%3),dx[d]);
		      // when using implicit source treatment do not include the current
		      if (sourceMethod!="IM" && d==0)
			{
			  arr(i,j,k,EX) -= dt*J_nPlusHalf(i,j,k,0)/(lambda_d*lambda_d*l_r);
			  arr(i,j,k,EY) -= dt*J_nPlusHalf(i,j,k,1)/(lambda_d*lambda_d*l_r);
			  arr(i,j,k,EZ) -= dt*J_nPlusHalf(i,j,k,2)/(lambda_d*lambda_d*l_r);
			}
		      /*if (geom.Coord()==1 && d==0)
			{
			  const Real y = geom.ProbLo()[1]+(double(j)+0.5)*dx[1];
			  arr(i,j,k,BX) -= dt*(arr(i,j,k,EZ)/y);      
			  arr(i,j,k,EX) += dt*(c*c*arr(i,j,k,BZ)/y);
			  }*/
		      arr(i,j,k,EX+d) += 0.0;
		      arr(i,j,k,EX+(1+d)%3) -= tau*dt*c*c*dxBz;
		      arr(i,j,k,EX+(2+d)%3) += tau*dt*c*c*dxBy;

		      if (RKOrder==2)
			{
			  Real dxByOld = computeDerivative(arr_old(i-iOffset,j-jOffset,k-kOffset,BX+(1+d)%3),
							   arr_old(i+iOffset,j+jOffset,k+kOffset,BX+(1+d)%3),dx[d]);
			  Real dxBzOld = computeDerivative(arr_old(i-iOffset,j-jOffset,k-kOffset,BX+(2+d)%3),
							   arr_old(i+iOffset,j+jOffset,k+kOffset,BX+(2+d)%3),dx[d]);
			  
			  arr(i,j,k,EX+d) += 0.0;
			  arr(i,j,k,EX+(1+d)%3) -= (1.0-tau)*dt*c*c*dxBzOld;
			  arr(i,j,k,EX+(2+d)%3) += (1.0-tau)*dt*c*c*dxByOld;
			}
		    }
		}
	    }
	}
    }
  // We need to compute boundary conditions again after each update 
  Sborder.FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries
  FillDomainBoundary(Sborder, geom, bc);  
}
void CAMReXmp::computeCurrentUpdate(MultiFab& J, const MultiFab& Sborder, const Real* dx, Real dt)
{
  MFIter::allowMultipleMFIters(true);
  
  const int iOffset = 1;
  const int jOffset = (amrex::SpaceDim!=1 ? 1 : 0);
  const int kOffset = (amrex::SpaceDim==3 ? 1 : 0);

  for(MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& J_nPlusHalf = J.array(mfi);

      for(int k = lo.z-kOffset; k <= hi.z+kOffset; k++)
	{
	  for(int j = lo.y-jOffset; j <= hi.y+jOffset; j++)
	    {
	      for(int i = lo.x-iOffset; i <= hi.x+iOffset; i++)
		{
		  Vector<Real> currentUpdated = currentUpdate(arr, i, j, k, dx, dt);
		  J_nPlusHalf(i,j,k,0) = currentUpdated[0];
		  J_nPlusHalf(i,j,k,1) = currentUpdated[1];
		  J_nPlusHalf(i,j,k,2) = currentUpdated[2];
		  /*J_nPlusHalf(i,j,k,0) = computeJx_nPlusHalf(arr, i, j, k, dx[0], dx[1], 2.0*dt);
		  J_nPlusHalf(i,j,k,1) = computeJy_nPlusHalf(arr, i, j, k, dx[0], dx[1], 2.0*dt);
		  J_nPlusHalf(i,j,k,2) = computeJz_nPlusHalf(arr, i, j, k, dx[0], dx[1], 2.0*dt);*/
		}
	    }
	}
    }
}

void CAMReXmp::setDomainBC (std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
			       std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc, int index)
{
  const BCRec& bc = get_desc_lst()[Phi_Type].getBC(index);
  Geometry const* gg = AMReX::top()->getDefaultGeometry();

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      if (gg->isPeriodic(idim))
	{
	  mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
	}
      else
	{
	  int pbc = bc.lo(idim);

	  if (pbc == EXT_DIR)
	    {
	      mlmg_lobc[idim] = LinOpBCType::Dirichlet;
	    }
	  else if (pbc == FOEXTRAP      ||
		   pbc == HOEXTRAP      || 
		   pbc == REFLECT_EVEN)
	    {
	      mlmg_lobc[idim] = LinOpBCType::Neumann;
	    }
	  else if (pbc == REFLECT_ODD)
	    {
	      mlmg_lobc[idim] = LinOpBCType::reflect_odd;
	    }
	  else
	    {
	      mlmg_lobc[idim] = LinOpBCType::bogus;
	    }
	  
	  pbc = bc.hi(idim);
	  if (pbc == EXT_DIR)
	    {
	      mlmg_hibc[idim] = LinOpBCType::Dirichlet;
	    }
	  else if (pbc == FOEXTRAP      ||
		   pbc == HOEXTRAP      || 
		   pbc == REFLECT_EVEN)
	    {
	      mlmg_hibc[idim] = LinOpBCType::Neumann;
	    }
	  else if (pbc == REFLECT_ODD)
	    {
	      mlmg_hibc[idim] = LinOpBCType::reflect_odd;
	    }
	  else
	    {
	      mlmg_hibc[idim] = LinOpBCType::bogus;
	    }
	}
    }
}

// Compute alpha using a fortran routine
void CAMReXmp::computeAlpha(MultiFab& alpha)
{
  MFIter::allowMultipleMFIters(true);
  
  alpha.define(grids, dmap, 1, 0);
  
#ifdef _OPENMP
#pragma omp parallel
#endif
  for(MFIter mfi(alpha, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
  
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
	  
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = alpha.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  arr(i,j,k) = 1.0;
		}
	    }
	}
    }  
}


// Compute beta using a fortran routine
void CAMReXmp::computeBeta(std::array<MultiFab,AMREX_SPACEDIM>& bcoeffs)
{
  MFIter::allowMultipleMFIters(true);
  
  MultiFab bTmp(grids, dmap, 1, 1);
  
  for(int n = 0; n < BL_SPACEDIM; n++)
    {
      const BoxArray& ba = convert(bTmp.boxArray(), IntVect::TheDimensionVector(n));
      bcoeffs[n].define(ba, bTmp.DistributionMap(), 1, 0);
    }
  
  for(int n = 0; n < BL_SPACEDIM; n++)
    {
      for(MFIter mfi(bcoeffs[n], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);
	  
	  // Indexable arrays for the data, and the directional flux
	  // Based on the vertex-centred definition of the flux array, the
	  // data array runs from e.g. [0,N] and the flux array from [0,N+1]
	  const auto& arr = bcoeffs[n].array(mfi);
	  
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      arr(i,j,k) = 1.0;
		    }
		}
	    }
	  
	}
    }
}
void CAMReXmp::computeRhs(std::array<MultiFab, 3>& Rhs, MultiFab& current, MultiFab& state,
			  const Real* dx, Real dt)
{
  MFIter::allowMultipleMFIters(true);
  /*
  Rhs[0].define(grids, dmap, 1, 0);
  Rhs[1].define(grids, dmap, 1, 0);
  Rhs[2].define(grids, dmap, 1, 0);
  */
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (int d = 0; d < amrex::SpaceDim ; d++)
    {
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      for(MFIter mfi(Rhs[0], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  // Indexable arrays for the data, and the directional flux
	  // Based on the vertex-centred definition of the flux array, the
	  // data array runs from e.g. [0,N] and the flux array from [0,N+1]
	  const auto& rhsX = Rhs[0].array(mfi);
	  const auto& rhsY = Rhs[1].array(mfi);
	  const auto& rhsZ = Rhs[2].array(mfi);
	  const auto& arr = state.array(mfi);
	  const auto& J_nPlusHalf = current.array(mfi);

	  // array used to update the rhs variables
	  std::array<Real, 3> rhsArray = {0.0, 0.0, 0.0};
	  
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      // contribution for Rhs-Y (or corresponding cylic relations, Z for d==1, X for d==2)
		      Real dxJz_nPlusHalf = computeDerivative(J_nPlusHalf(i-iOffset,j-jOffset,k-kOffset,(2+d)%3), J_nPlusHalf(i+iOffset,j+jOffset,k+kOffset,(2+d)%3), dx[d]);
		      Real dxEz = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,EX+(2+d)%3), arr(i+iOffset,j+jOffset,k+kOffset,EX+(2+d)%3), dx[d]);
		      // contribution for Rhs-Z (or corresponding cylic relations, X for d==1, Y for d==2) 
		      Real dxJy_nPlusHalf = computeDerivative(J_nPlusHalf(i-iOffset,j-jOffset,k-kOffset,(1+d)%3), J_nPlusHalf(i+iOffset,j+jOffset,k+kOffset,(1+d)%3), dx[d]);
		      Real dxEy = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,EX+(1+d)%3), arr(i+iOffset,j+jOffset,k+kOffset,EX+(1+d)%3), dx[d]);
		      
		      rhsArray[d] = 0.0;
		      rhsArray[(1+d)%3] = + dt*dxEz - tau*dt*dt*dxJz_nPlusHalf/(lambda_d*lambda_d*l_r);
		      rhsArray[(2+d)%3] = - dt*dxEy + tau*dt*dt*dxJy_nPlusHalf/(lambda_d*lambda_d*l_r);
		      
		      // update rhs variables
		      //if (d==0)
			{
			  //rhsX(i,j,k) = arr(i,j,k,BX);
			  //rhsY(i,j,k) = arr(i,j,k,BY);
			  //rhsZ(i,j,k) = arr(i,j,k,BZ);
			  /*
			  // add cylindrical source terms
			  if (geom.Coord()==1)
			    {
			      // rhsZ += -dt*Ez/x... in cartesian
			      // y <- x in cylindrical
			      const Real y = geom.ProbLo()[1]+(double(j)+0.5)*dx[1];
			      rhsX(i,j,k) += -dt*arr(i,j,k,EZ)/y + tau*dt*dt*J_nPlusHalf(i+1,j,k,2)/(lambda_d*lambda_d*l_r*y);
			      
			      // cylindrical source terms for the vector laplace operator
			      if (RKOrder==2)
				{
				  // rhsX -= ...Bx/x^2 in cartesian
				  // y <- x in cylindrical
				  //rhsY(i,j,k) -= tau*(1.0-tau)*c*c*dt*dt*arr(i,j,k,BY)/(y*y);
				  //rhsZ(i,j,k) -= tau*(1.0-tau)*c*c*dt*dt*arr(i,j,k,BZ)/(y*y);
				}
			    }
			  */
			}
		      
		      if (RKOrder==2)
			{
			  Real dxdxBx = computeSecondDerivative(arr(i-iOffset,j-jOffset,k-kOffset,BX+d), arr(i,j,k,BX+d),
								arr(i+iOffset,j+jOffset,k+kOffset,BX+d), dx[d]);
			  Real dxdxBy = computeSecondDerivative(arr(i-iOffset,j-jOffset,k-kOffset,BX+(1+d)%3), arr(i,j,k,BX+(1+d)%3),
								arr(i+iOffset,j+jOffset,k+kOffset,BX+(1+d)%3), dx[d]);
			  Real dxdxBz = computeSecondDerivative(arr(i-iOffset,j-jOffset,k-kOffset,BX+(2+d)%3), arr(i,j,k,BX+(2+d)%3),
								arr(i+iOffset,j+jOffset,k+kOffset,BX+(2+d)%3), dx[d]);

			  rhsArray[d] += tau*(1.0-tau)*c*c*dt*dt*dxdxBx;
			  rhsArray[(1+d)%3] += tau*(1.0-tau)*c*c*dt*dt*dxdxBy;
			  rhsArray[(2+d)%3] += tau*(1.0-tau)*c*c*dt*dt*dxdxBz;
			}

		      rhsX(i,j,k) += rhsArray[0];
		      rhsY(i,j,k) += rhsArray[1];
		      rhsZ(i,j,k) += rhsArray[2];		     

		    }
		}
	    }
	}
    }
}

// For debugging purposes - this will print all the data from a multifab
void CAMReXmp::printComponents(MultiFab& mfIn)
{
  MFIter::allowMultipleMFIters(true);
  
  for(MFIter mfi(mfIn, true); mfi.isValid(); ++mfi)
  {
    // Real* dat = mfIn[mfi].dataPtr();
    // const Box& abx = mfIn[mfi].box();
    // const int* alo = abx.loVect();
    // const int* ahi = abx.hiVect();

    // const Box& bx = mfi.tilebox();
    // const int* lo = bx.loVect();
    // const int* hi = bx.hiVect();

    // // amrex::Print() << "Printing components for box " << bx << "\n";
    // 
    // printComp(AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
	//       BL_TO_FORTRAN_3D(mfIn[mfi]));

    // Example: access MF elements without Fortran function
    // TODO: Remove the Fortran version 

    // Get the tilebox and its bounds
    Box bx = mfi.tilebox();
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    // This is the array structure that we can manipulate directly 
    const auto& arr = mfIn.array(mfi);

    // Somewhere near the middle of the y coordinate 
    int mid = (bx.bigEnd(1) + bx.smallEnd(1)) / 2;

    // Loop over x
    // Note that we could also loop over j from lo.y to hi.y and z from lo.z to hi.z
    // When doing multiple loops, we want to go k first, then j, then i (for memory layout)
    for(int i = lo.x; i <= hi.x; i++)
    {
        amrex::Print() << "At (" << i << "," << mid << "), " 
            << "component = " << arr(i,mid,0) << std::endl;
    }
  } 
}

void CAMReXmp::computeRhsX(MultiFab& Rhs, MultiFab& current, MultiFab& state,
			     Real time, Real dt, Real dx, Real dy)
{
  Rhs.define(grids, dmap, 1, 0);
  
#ifdef _OPENMP
#pragma omp parallel
#endif
  for(MFIter mfi(Rhs, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      /*Real* dat = temp[mfi].dataPtr();
      const Box& abx = Rhs[mfi].box();
      const int* alo = abx.loVect();
      const int* ahi = abx.hiVect();
      const int* lo  = bx.loVect();
      const int* hi  = bx.hiVect();
      const Real* powIdx = powerIdx.dataPtr();
      const Real* powVal = powerVal.dataPtr();
      const int  powSize = powerIdx.size();
      const Real* eIdx = emisIdx.dataPtr();
      const Real* eVal = emisVal.dataPtr();
      const int  eSize = emisIdx.size();*/
      //setRHS(AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
      //     BL_TO_FORTRAN_3D(Rhs[mfi]), BL_TO_FORTRAN_3D(alpha[mfi]), BL_TO_FORTRAN_3D(temp[mfi]),
      //     AMREX_ZFILL(geom.CellSize()),
      //     AMREX_ZFILL(geom.ProbLo()), AMREX_ZFILL(geom.ProbHi()), &time, &dt,
      //   powIdx, powVal, &powSize, eIdx, eVal, &eSize);

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs.array(mfi);
      const auto& arr = state.array(mfi);
      const auto& J_nPlusHalf = current.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  /////////////////////////////////////////////////////////////////////////
		  // 2D
		  //Real Jz_nPlusHalf_jMinus1 = computeJz_nPlusHalf(arr, i, j-1, k, dx, dy, dt);
		  //Real Jz_nPlusHalf_jPlus1 = computeJz_nPlusHalf(arr, i, j+1, k, dx, dy, dt);
		  Real dyJz_nPlusHalf = computeDerivative(J_nPlusHalf(i,j-1,k,2), J_nPlusHalf(i,j+1,k,2), dy);
		  Real dyEz = computeDerivative(arr(i,j-1,k,EZ), arr(i,j+1,k,EZ), dy);
		  //Real dxdxBx = computeSecondDerivative(arr(i-1,j,k,BX), arr(i,j,k,BX), arr(i+1,j,k,BX), dx);
		  //Real dydyBx = computeSecondDerivative(arr(i,j-1,k,BX), arr(i,j,k,BX), arr(i,j+1,k,BX), dx);		  
		  rhs(i,j,k) = arr(i,j,k,BX) - dt*dyEz + dt*dt*dyJz_nPlusHalf/(lambda_d*lambda_d*l_r);
		  /////////////////////////////////////////////////////////////////////////
		  // 1D
		  //rhs(i,j,k) = arr(i,j,k,BX);
		}
	    }
	}
    }
  
}
void CAMReXmp::computeRhsY(MultiFab& Rhs, MultiFab& current, MultiFab& state,
			     Real time, Real dt, Real dx, Real dy)
{
  Rhs.define(grids, dmap, 1, 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
  for(MFIter mfi(Rhs, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      /*Real* dat = Rhs[mfi].dataPtr();
      const Box& abx = Rhs[mfi].box();
      const int* alo = abx.loVect();
      const int* ahi = abx.hiVect();
      const int* lo  = bx.loVect();
      const int* hi  = bx.hiVect();
      const Real* powIdx = powerIdx.dataPtr();
      const Real* powVal = powerVal.dataPtr();
      const int  powSize = powerIdx.size();
      const Real* eIdx = emisIdx.dataPtr();
      const Real* eVal = emisVal.dataPtr();
      const int  eSize = emisIdx.size();*/
      //setRHS(AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
      //     BL_TO_FORTRAN_3D(Rhs[mfi]), BL_TO_FORTRAN_3D(alpha[mfi]), BL_TO_FORTRAN_3D(temp[mfi]),
      //     AMREX_ZFILL(geom.CellSize()),
      //     AMREX_ZFILL(geom.ProbLo()), AMREX_ZFILL(geom.ProbHi()), &time, &dt,
      //   powIdx, powVal, &powSize, eIdx, eVal, &eSize);

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
	  
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs.array(mfi);
      const auto& arr = state.array(mfi);
      const auto& J_nPlusHalf = current.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  //Real Jz_nPlusHalf_iMinus1 = computeJz_nPlusHalf(arr, i-1, j, k, dx, dy, dt);
		  //Real Jz_nPlusHalf_iPlus1 = computeJz_nPlusHalf(arr, i+1, j, k, dx, dy, dt);
		  Real dxJz_nPlusHalf = computeDerivative(J_nPlusHalf(i-1,j,k,2), J_nPlusHalf(i+1,j,k,2), dx);
		  Real dxEz = computeDerivative(arr(i-1,j,k,EZ), arr(i+1,j,k,EZ), dx);
		  //Real dxdxBy = computeSecondDerivative(arr(i-1,j,k,BY), arr(i,j,k,BY), arr(i+1,j,k,BY), dx);
		  //Real dydyBy = computeSecondDerivative(arr(i,j-1,k,BY), arr(i,j,k,BY), arr(i,j+1,k,BY), dx);
		  //std::cout << dxEz << " " << dxJz_nPlusHalf << std::endl;
		  rhs(i,j,k) = arr(i,j,k,BY) + dt*dxEz - dt*dt*dxJz_nPlusHalf/(lambda_d*lambda_d*l_r);
		}
	    }
	}
    }
}
void CAMReXmp::computeRhsZ(MultiFab& Rhs, MultiFab& current, MultiFab& state,
			     Real time, Real dt, Real dx, Real dy)
{
  Rhs.define(grids, dmap, 1, 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
  for(MFIter mfi(Rhs, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      /*Real* dat = Rhs[mfi].dataPtr();
      const Box& abx = Rhs[mfi].box();
      const int* alo = abx.loVect();
      const int* ahi = abx.hiVect();
      const int* lo  = bx.loVect();
      const int* hi  = bx.hiVect();
      const Real* powIdx = powerIdx.dataPtr();
      const Real* powVal = powerVal.dataPtr();
      const int  powSize = powerIdx.size();
      const Real* eIdx = emisIdx.dataPtr();
      const Real* eVal = emisVal.dataPtr();
      const int  eSize = emisIdx.size();*/
      //setRHS(AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
      //     BL_TO_FORTRAN_3D(Rhs[mfi]), BL_TO_FORTRAN_3D(alpha[mfi]), BL_TO_FORTRAN_3D(temp[mfi]),
      //     AMREX_ZFILL(geom.CellSize()),
      //     AMREX_ZFILL(geom.ProbLo()), AMREX_ZFILL(geom.ProbHi()), &time, &dt,
      //   powIdx, powVal, &powSize, eIdx, eVal, &eSize);

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
	  
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs.array(mfi);
      const auto& arr = state.array(mfi);
      const auto& J_nPlusHalf = current.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  //Real Jx_nPlusHalf_jMinus1 = computeJx_nPlusHalf(arr, i, j-1, k, dx, dy, dt);
		  //Real Jx_nPlusHalf_jPlus1 = computeJx_nPlusHalf(arr, i, j+1, k, dx, dy, dt);
		  //Real Jy_nPlusHalf_iMinus1 = computeJy_nPlusHalf(arr, i-1, j, k, dx, dy, dt);
                  //Real Jy_nPlusHalf_iPlus1 = computeJy_nPlusHalf(arr, i+1, j, k, dx, dy, dt);
		  /////////////////////////////////////////////////////////////////////////
		  // 2D
		  Real dyJx_nPlusHalf = computeDerivative(J_nPlusHalf(i,j-1,k,0), J_nPlusHalf(i,j+1,k,0), dy);
		  Real dyEx = computeDerivative(arr(i,j-1,k,EX), arr(i,j+1,k,EX), dy);
		  ///////////////////////////////////////////////////////////////////////// 
		  Real dxJy_nPlusHalf = computeDerivative(J_nPlusHalf(i-1,j,k,1), J_nPlusHalf(i+1,j,k,1), dx);
		  Real dxEy = computeDerivative(arr(i-1,j,k,EY), arr(i+1,j,k,EY), dx);
		  //Real dxdxBz = computeSecondDerivative(arr(i-1,j,k,BZ), arr(i,j,k,BZ), arr(i+1,j,k,BZ), dx);
		  //Real dydyBz = computeSecondDerivative(arr(i,j-1,k,BZ), arr(i,j,k,BZ), arr(i,j+1,k,BZ), dx);
		  /////////////////////////////////////////////////////////////////////////
		  // 2D
		  rhs(i,j,k) = arr(i,j,k,BZ) - dt*(dxEy-dyEx) + dt*dt*(dxJy_nPlusHalf-dyJx_nPlusHalf)/(lambda_d*lambda_d*l_r);
		  /////////////////////////////////////////////////////////////////////////
		  // 1D
		  //rhs(i,j,k) = arr(i,j,k,BZ) - dt*dxEy + dt*dt*dxJy_nPlusHalf/(lambda_d*lambda_d*l_r); 
		}
	    }
	}
    }
}
void CAMReXmp::hyperbolicMaxwellSolverDivFree(Array<MultiFab,AMREX_SPACEDIM>& S_EM, MultiFab& fluxesEM, MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{

  // Update cell-centred z-components of EM fields
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
	  //const auto& arrEM = S_EM[d].array(mfi);      
      
	  for(int k = lo.z; k <= hi.z+kOffset; k++)
	    {
	      for(int j = lo.y; j <= hi.y+jOffset; j++)
		{
		  for(int i = lo.x; i <= hi.x+iOffset; i++)
		    {	    
		      // Fluxes for the z-components of the EM fields because it is 2D code
		      Vector<Real> fluxEM = Maxwell_flux_Godunov(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);
		      fluxArr(i,j,k,BZ) = fluxEM[BZ_LOCAL];
		      fluxArr(i,j,k,EZ) = fluxEM[EZ_LOCAL];
	    
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
		      // Update cell-centred z-components becuause it is 2D code
		      arr(i,j,k,BZ) = arr(i,j,k,BZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, BZ) - fluxArr(i,j,k,BZ));
		      arr(i,j,k,EZ) = arr(i,j,k,EZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, EZ) - fluxArr(i,j,k,EZ));
	   
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
    
    }
  
  // Compute cell-centred Ez source terms
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());

      Array4<Real> arr = Sborder.array(mfi);
      Array4<Real> arrEMX = S_EM[0].array(mfi);
      Array4<Real> arrEMY = S_EM[1].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
		}
	    }
	}      
    }

  // We need to compute boundary conditions again after each update
  Sborder.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries
  FillDomainBoundary(Sborder, geom, bc);
  
    
  // Compute edge components of the EM fields    
  for (MFIter mfi(fluxesEM, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi); 
      const auto& fluxArrEM = fluxesEM.array(mfi);
      //const auto& fluxArrX = fluxes[0].array(mfi);
      //const auto& fluxArrY = fluxes[1].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{		 
		  // intermediate values using y-face in the x-direction
		  Vector<Real> cornerEM_Y = Maxwell_corner(arrEM_Y, i, j, k, 1, 0, 0, dx[0], dt, 0);

		  // intermediate values using x-face in the y-direction
		  Vector<Real> cornerEM_X = Maxwell_corner(arrEM_X, i, j, k, 0, 1, 0, dx[1], dt, 1);	       
		  
		  fluxArrEM(i,j,k,BX_LOCAL) = cornerEM_X[BX_LOCAL];
		  fluxArrEM(i,j,k,BY_LOCAL) = cornerEM_Y[BY_LOCAL];
		  fluxArrEM(i,j,k,BZ_LOCAL) = cornerEM_X[BZ_LOCAL] + cornerEM_Y[BZ_LOCAL];
		  fluxArrEM(i,j,k,EX_LOCAL) = cornerEM_X[EX_LOCAL];
		  fluxArrEM(i,j,k,EY_LOCAL) = cornerEM_Y[EY_LOCAL];
		  fluxArrEM(i,j,k,EZ_LOCAL) = cornerEM_X[EZ_LOCAL] + cornerEM_Y[EZ_LOCAL];
		  
		}
	    }
	}       
    }

  // Could be useful in a 3D code
  // // Compute edge components of the EM fields
  // for (int d = 0; d < amrex::SpaceDim ; d++)   
  // {

  //   const int iOffset = ( d == 0 ? 1 : 0);
  //   const int jOffset = ( d == 1 ? 1 : 0);
  //   const int kOffset = ( d == 2 ? 1 : 0);

  //   int d_EM = (d==0) ? 1 : 0;
    
  //   // Loop over all the patches at this level
  //   for (MFIter mfi(S_EM[d_EM], true); mfi.isValid(); ++mfi)
  //   {
  //     const Box& bx = mfi.tilebox();

  //     const Dim3 lo = lbound(bx);
  //     const Dim3 hi = ubound(bx);

  //     // Indexable arrays for the data, and the directional flux
  //     // Based on the corner-centred definition of the flux array, the
  //     // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
  //     const auto& arrEM = S_EM[d_EM].array(mfi);
  //     const auto& fluxArrEM = fluxesEM.array(mfi);
  //     const auto& fluxArr = fluxes[d_EM].array(mfi);
      
  //     //std::cout << fluxArrEM << " " << arrEM << " " << lo.x << " " << hi.x+iOffset << " " << lo.y << " " << hi.y+jOffset << std::endl;      
  //     //amrex::Abort();
      
  //     for(int k = lo.z; k <= hi.z+kOffset; k++)
  //     {
  // 	for(int j = lo.y; j <= hi.y+jOffset; j++)
  // 	{
  // 	  for(int i = lo.x; i <= hi.x+iOffset; i++)
  // 	  {
  // 	    // electromagnetic field values reconstructed at the corners
  // 	    Vector<Real> cornerEM = Maxwell_corner(arrEM, i, j, k, iOffset, jOffset, kOffset, d);
	    
  // 	    fluxArrEM(i,j,k,BX_LOCAL+(1+d)%3) += cornerEM[BX_LOCAL+(1+d)%3];
  // 	    fluxArrEM(i,j,k,BX_LOCAL+(2+d)%3) += cornerEM[BX_LOCAL+(2+d)%3];
  // 	    fluxArrEM(i,j,k,EX_LOCAL+(1+d)%3) += cornerEM[EX_LOCAL+(1+d)%3];
  // 	    fluxArrEM(i,j,k,EX_LOCAL+(2+d)%3) += cornerEM[EX_LOCAL+(2+d)%3];

  // 	    /*
  // 	    if ((i==256 && j>252))
  // 	    std::cout << d << " " << j << " " << arrEM(i,j,k,EX_LOCAL+d_EM) << " " << -std::pow(-1,d)*c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) - fluxArrEM(i,j,k,BZ_LOCAL)) << " " << fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) << " " << fluxArrEM(i,j,k,BZ_LOCAL) << " " << -dt*1.0/(lambda_d*lambda_d*l_r)*currentFace << " " << fluxArr(i,j,k,0) << " " << fluxArr(i,j,k,5) << std::endl;*/
  // 	  }
  // 	}
  //     }
  //   }
      
  //   // We need to compute boundary conditions again after each update
  //   //S_EM[0].FillBoundary(geom.periodicity());
  //   //S_EM[1].FillBoundary(geom.periodicity());
     
  //   // added by 2020D 
  //   // Fill non-periodic physical boundaries                          
  //   //FillDomainBoundary(S_EM[0], geom, bc_EM);
  //   //FillDomainBoundary(S_EM[1], geom, bc_EM);
    
  // }


  // Update face-centred EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
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
	  //const auto& arr = Sborder.array(mfi);
	  const auto& arrEM = S_EM[d_EM].array(mfi);
	  //const auto& arrEMother = S_EM[d].array(mfi);
	  const auto& fluxArrEM = fluxesEM.array(mfi);
	  const auto& fluxArr = fluxes[d_EM].array(mfi);

	  const Dim3 hiDomain = ubound(geom.Domain());
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      // 2D code; z-component updated using cell-centred scheme
		      arrEM(i,j,k,BX_LOCAL+d_EM) += std::pow(-1,d)*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EZ_LOCAL) - fluxArrEM(i,j,k,EZ_LOCAL));
		      arrEM(i,j,k,EX_LOCAL+d_EM) -= std::pow(-1,d)*c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) - fluxArrEM(i,j,k,BZ_LOCAL));

		      // 3D code
		      /*arrEM(i,j,k,BX_LOCAL+(1+d)%3) += (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(2+d)%3));
			arrEM(i,j,k,BX_LOCAL+(2+d)%3) -= (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(1+d)%3));
			arrEM(i,j,k,EX_LOCAL+(1+d)%3) -= c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(2+d)%3));
			arrEM(i,j,k,EX_LOCAL+(2+d)%3) += c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(1+d)%3));
		      */
	    
		      // source terms
		      Real currentFace = r_i*fluxArr(i,j,k,RHO_I) + r_e*fluxArr(i,j,k,RHO_E);
		      arrEM(i,j,k,EX_LOCAL+d_EM) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;

		      //if (d_EM==1 && i==256 && (j>252))
		      //std::cout << arrEM(i,j,k,EX_LOCAL+d_EM) << " " << currentFace << std::endl;
		      // 3D code
		      //arrEM(i,j,k,EX_LOCAL+(1+d)%3) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;
		      //arrEM(i,j,k,EX_LOCAL+(2+d)%3) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;

		      // NOTE
		      // In future, when debugging, check that the flux (current) is zero at the boundaries
		      // i==1 and not i==0 because we need value at i==1 to be updated first
		      /*
			if (i==1 && bc[MOMX_I+d_EM].lo(0) == BCType::reflect_odd)
			arrEM(0,j,k,EX_LOCAL+d_EM) = arrEM(1,j,k,EX_LOCAL+d_EM);	      
			if (i==hiDomain.x+1 && bc[MOMX_I+d_EM].hi(0) == BCType::reflect_odd)
			arrEM(hiDomain.x+1,j,k,EX_LOCAL+d_EM) = arrEM(hiDomain.x,j,k,EX_LOCAL+d_EM);	    
			if (j==1 && bc[MOMX_I+d_EM].lo(1) == BCType::reflect_odd)
			arrEM(i,0,k,EX_LOCAL+d_EM) = arrEM(i,1,k,EX_LOCAL+d_EM);
			if (j==hiDomain.y+1 && bc[MOMX_I+d_EM].hi(1) == BCType::reflect_odd)
			arrEM(i,hiDomain.y+1,k,EX_LOCAL+d_EM) = arrEM(i,hiDomain.y,k,EX_LOCAL+d_EM);
		      */
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
  /*
  // check if the current is zero at the domain boundaries
  if (bc[MOMX_I].lo(0) == BCType::reflect_odd)
    {
      for (MFIter mfi(S_EM[0], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Dim3 hiDomain = ubound(geom.Domain());

	  Array4<Real> arr = Sborder.array(mfi);
	  Array4<Real> arrEMX = S_EM[0].array(mfi);
	  Array4<Real> arrEMY = S_EM[1].array(mfi);

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      if (i==0)
			//arrEMX(0,j,k,EX_LOCAL) = arrEMX(1,j,k,EX_LOCAL) -
			//dx[0]*(arr(i,j,k,DIVE)-(arrEMY(i,j+1,k,EY_LOCAL)-arrEMY(i,j,k,EY_LOCAL))/dx[1]);
			//arrEMX(0,j,k,EX_LOCAL) = arrEMX(1,j,k,EX_LOCAL);
			arrEMX(0,j,k,EX_LOCAL) = 0.0;
		    }
		}
	    }      
	}      
    }
  if (bc[MOMX_I].hi(0) == BCType::reflect_odd)
    {
      for (MFIter mfi(S_EM[0], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Dim3 hiDomain = ubound(geom.Domain());

	  Array4<Real> arr = Sborder.array(mfi);
	  Array4<Real> arrEMX = S_EM[0].array(mfi);
	  Array4<Real> arrEMY = S_EM[1].array(mfi);

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      if (i==hiDomain.x)
			//arrEMX(hiDomain.x+1,j,k,EX_LOCAL) = arrEMX(hiDomain.x,j,k,EX_LOCAL) +
			//dx[0]*(arr(i,j,k,DIVE)-(arrEMY(i,j+1,k,EY_LOCAL)-arrEMY(i,j,k,EY_LOCAL))/dx[1]);
			//arrEMX(hiDomain.x+1,j,k,EX_LOCAL) = arrEMX(hiDomain.x,j,k,EX_LOCAL);
			arrEMX(hiDomain.x+1,j,k,EX_LOCAL) = 0.0;
		    }
		}
	    }      
	}      
    }
  if (bc[MOMY_I].lo(1) == BCType::reflect_odd)
    {
      for (MFIter mfi(S_EM[1], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Dim3 hiDomain = ubound(geom.Domain());

	  Array4<Real> arr = Sborder.array(mfi);
	  Array4<Real> arrEMX = S_EM[0].array(mfi);
	  Array4<Real> arrEMY = S_EM[1].array(mfi);

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      if (j==0)
			//arrEMY(i,0,k,EY_LOCAL) = arrEMY(i,1,k,EY_LOCAL) -
			//dx[1]*(arr(i,j,k,DIVE)-(arrEMX(i+1,j,k,EX_LOCAL)-arrEMX(i,j,k,EX_LOCAL))/dx[0]);
			//arrEMY(i,0,k,EY_LOCAL) = arrEMY(i,1,k,EY_LOCAL);
			arrEMY(i,0,k,EY_LOCAL) = 0.0;
		    }
		}
	    }      
	}      
    }
  if (bc[MOMY_I].hi(1) == BCType::reflect_odd)
    {
      for (MFIter mfi(S_EM[1], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Dim3 hiDomain = ubound(geom.Domain());

	  Array4<Real> arr = Sborder.array(mfi);
	  Array4<Real> arrEMX = S_EM[0].array(mfi);
	  Array4<Real> arrEMY = S_EM[1].array(mfi);

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {		      		      
		      if (j==hiDomain.y)
			//arrEMY(i,hiDomain.y+1,k,EY_LOCAL) = arrEMY(i,hiDomain.y,k,EY_LOCAL) +
			//dx[1]*(arr(i,j,k,DIVE)-(arrEMX(i+1,j,k,EX_LOCAL)-arrEMX(i,j,k,EX_LOCAL))/dx[0]);
			//arrEMY(i,hiDomain.y+1,k,EY_LOCAL) = arrEMY(i,hiDomain.y,k,EY_LOCAL);
			arrEMY(i,hiDomain.y+1,k,EY_LOCAL) = 0.0;
		    }
		}
	    }      
	}      
    }
  */
  // We need to compute boundary conditions again after each update
  S_EM[0].FillBoundary(geom.periodicity());
  S_EM[1].FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM[0], geom, bc_EM);    
  FillDomainBoundary(S_EM[1], geom, bc_EM);  

}
void CAMReXmp::hyperbolicMaxwellSolverDivFreeBalsara(Array<MultiFab,AMREX_SPACEDIM>& S_EM, MultiFab& fluxesEM, MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], MultiFab& S0, const Real* dx, Real dt)
{
  // Future slopes definition suggestion
  // Useful when dealing with high-order reconstruction
  // To get first-order slope in y-direction for Bx: slopes[BX_LOCALL].array(i,j,k,Y)
  // Useful to define indices int X,Y,Z in second-order linear reconstruction
  // For third-order quadratic reconstruction use int X,Y,Z,XX,YY,ZZ,XY,YZ,XZ
  // Although for x-components fields, only exist Y,Z,YY,ZZ,YZ (see Balsara et al. 2017)
  /*
  Array<MultiFab,6> slopes;
  // number of slopes
  // 2 slopes for second-order linear reconstruction
  // Ex = Ex^0 + Ex^y*(y/Delta_y) + Ex^z*(z/Delta_z), where Ex^y and Ex^z are the slopes
  // For 2D code, only 1 slope for x- and y-compoenents
  // For clarity, will use 3 slopes for second-order, first one represents slope in x-direction
  // 9 slopes for third-order, and so on
  int nSlopes = 3;
  slopes[BX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);
  slopes[BY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
  slopes[EX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);

  // For 2D code, 2 slopes for z-components, although it is cell-centred
  slopes[BZ_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[BZ_LOCAL].define(grids, dmap, nSlopes, 1);
  */

  Array<MultiFab,AMREX_SPACEDIM> slopesFC;
  // six components, one for each field direction
  slopesFC[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, 1);
  slopesFC[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, 1);

  // only for 2D
  Array<MultiFab,AMREX_SPACEDIM> slopesCC;
  slopesCC[0].define(grids, dmap, 6, 1);
  slopesCC[1].define(grids, dmap, 6, 1);
  
  // Three components, one for each direction
  MultiFab slopesCharge(grids, dmap, AMREX_SPACEDIM, 1);

  // Compute cell-centred density slopes
  for (MFIter mfi(slopesCharge, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // uses old charge
      const Array4<Real> arr = S0.array(mfi);
      Array4<Real> slopes = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  Real u_iMinus1 = (r_i*arr(i-1,j,k,RHO_I) + r_e*arr(i-1,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		  Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		  Real u_iPlus1 = (r_i*arr(i+1,j,k,RHO_I) + r_e*arr(i+1,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		  Real u_jMinus1 = (r_i*arr(i,j-1,k,RHO_I) + r_e*arr(i,j-1,k,RHO_E))/(lambda_d*lambda_d*l_r);
		  Real u_jPlus1 = (r_i*arr(i,j+1,k,RHO_I) + r_e*arr(i,j+1,k,RHO_E))/(lambda_d*lambda_d*l_r);
		  // slope measure
		  Vector<Real> delta_i = get_delta_i({u_iMinus1,u_jMinus1},{u_i,u_i},{u_iPlus1,u_jPlus1});
		  // slope ratio
		  Real ri = get_r(u_iMinus1,u_i,u_iPlus1);
		  Real rj = get_r(u_jMinus1,u_i,u_jPlus1);
		  // slope limiter
		  Real epsilon_i = get_epsilon(ri);
		  Real epsilon_j = get_epsilon(rj);

		  // slopes for the current
		  slopes(i,j,k,0) = epsilon_i*delta_i[0];
		  slopes(i,j,k,1) = epsilon_j*delta_i[1];

		}
	    }
	}      
    }
  
  // Compute cell-centred z-components slopes
  for (MFIter mfi(slopesCC[0], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Array4<Real> arr = Sborder.array(mfi);
      Array4<Real> slopesX = slopesCC[0].array(mfi);
      Array4<Real> slopesY = slopesCC[1].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  for (int n : {BZ,EZ}){
		    Real u_iMinus1 = arr(i-1,j,k,n);
		    Real u_i = arr(i,j,k,n);
		    Real u_iPlus1 = arr(i+1,j,k,n);
		    Real u_jMinus1 = arr(i,j-1,k,n);
		    Real u_jPlus1 = arr(i,j+1,k,n);
		    // slope measure
		    Vector<Real> delta_i = get_delta_i({u_iMinus1,u_jMinus1},{u_i,u_i},{u_iPlus1,u_jPlus1});
		    // slope ratio
		    Real ri = get_r(u_iMinus1,u_i,u_iPlus1);
		    Real rj = get_r(u_jMinus1,u_i,u_jPlus1);
		    // slope limiter
		    Real epsilon_i = get_epsilon(ri);
		    Real epsilon_j = get_epsilon(rj);

		    // slopes for the z-components in x and y directions
		    slopesX(i,j,k,n-NUM_STATE_FLUID) = epsilon_i*delta_i[0];
		    slopesY(i,j,k,n-NUM_STATE_FLUID) = epsilon_j*delta_i[1];
		    
		  }
		}
	    }
	}      
    }

  // Compute y-components slopes in x-direction
  for (MFIter mfi(slopesFC[1], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Array4<Real> arr = S_EM[1].array(mfi);
      Array4<Real> slopes = slopesFC[1].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  for (int n : {BY_LOCAL,EY_LOCAL}){
		    Real u_iMinus1 = arr(i-1,j,k,n);
		    Real u_i = arr(i,j,k,n);		    
		    Real u_iPlus1 = arr(i+1,j,k,n);
		    /*// slope measure
		    Real delta_i = get_delta_i(u_iMinus1,u_i,u_iPlus1);
		    // slope ratio
		    Real ri = get_r(u_iMinus1,u_i,u_iPlus1);
		    // slope limiter
		    Real epsilon_i = get_epsilon(ri);
		    */
		    Real slope = TVD_slope(u_iMinus1,u_i,u_iPlus1);
		    // slopes
		    //slopes(i,j,k,n) = epsilon_i*delta_i;
		    slopes(i,j,k,n) = TVD_slope(u_iMinus1,u_i,u_iPlus1); 
		  }
		}
	    }
	}      
    }
  // Compute x-components slopes in y-direction
  for (MFIter mfi(slopesFC[0], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Array4<Real> arr = S_EM[0].array(mfi);
      Array4<Real> slopes = slopesFC[0].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  for (int n : {BX_LOCAL,EX_LOCAL}){
		    Real u_i = arr(i,j,k,n);
		    Real u_jMinus1 = arr(i,j-1,k,n);
		    Real u_jPlus1 = arr(i,j+1,k,n);
		    // slope measure
		    Real delta_i = get_delta_i(u_jMinus1,u_i,u_jPlus1);
		    // slope ratio
		    Real rj = get_r(u_jMinus1,u_i,u_jPlus1);
		    // slope limiter
		    Real epsilon_j = get_epsilon(rj);

		    // slopes
		    slopes(i,j,k,n) = epsilon_j*delta_i;
		  }
		}
	    }
	}      
    }  

  // Compute edge components of the EM fields    
  for (MFIter mfi(fluxesEM, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      // use old array, Sborder should also work
      //const auto& arr = Sborder.array(mfi);
      const auto& arr = S0.array(mfi);
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi); 
      const auto& fluxArrEM = fluxesEM.array(mfi);

      const auto& slopesFC_X = slopesFC[0].array(mfi);
      const auto& slopesFC_Y = slopesFC[1].array(mfi);

      // only for 2D
      const auto& slopesCC_X = slopesCC[0].array(mfi);
      const auto& slopesCC_Y = slopesCC[1].array(mfi);
  
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{		 
		  
		  Real a0,ax,ay,az,axx,axy,axz,b0,bx,by,bz,bxy,byy,byz,c0,cx,cy,cz,cxz,cyz,czz;
		  Real x,y;

		  // Magnetic field
		  // LD state
		  x = dx[0]/2.0, y = dx[1]/2.0;
		  ay = (slopesFC_X(i,j-1,k,BX_LOCAL)+slopesFC_X(i-1,j-1,k,BX_LOCAL))/2.0;
		  axy = (slopesFC_X(i,j-1,k,BX_LOCAL)-slopesFC_X(i-1,j-1,k,BX_LOCAL));
		  az = 0.0;
		  axz = 0.0;
		  bx = (slopesFC_Y(i-1,j,k,BY_LOCAL)+slopesFC_Y(i-1,j-1,k,BY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i-1,j,k,BY_LOCAL)-slopesFC_Y(i-1,j-1,k,BY_LOCAL));
		  bz = 0.0;
		  byz = 0.0;
		  cx = slopesCC_X(i-1,j-1,k,BZ_LOCAL);
		  cxz = 0.0;
		  cy = slopesCC_Y(i-1,j-1,k,BZ_LOCAL);
		  cyz = 0.0;
		  axx = -dx[0]/2.0*(bxy/dx[1]); // +cxz/dx[2]
		  byy = -dx[1]/2.0*(axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i,j-1,k,BX_LOCAL)+arrEM_X(i-1,j-1,k,BX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i,j-1,k,BX_LOCAL)-arrEM_X(i-1,j-1,k,BX_LOCAL));
		  b0 = (arrEM_Y(i-1,j,k,BY_LOCAL)+arrEM_Y(i-1,j-1,k,BY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i-1,j,k,BY_LOCAL)-arrEM_Y(i-1,j-1,k,BY_LOCAL));
		  c0 = arr(i-1,j-1,k,BZ) - czz/6.0; //czz*dx[2]*dx[2]/6.0
		  cz = 0.0;
		  Real BxLD = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real ByLD = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real BzLD = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  
		  
		  // RD state
		  x = -dx[0]/2.0, y = dx[1]/2.0;
		  ay = (slopesFC_X(i+1,j-1,k,BX_LOCAL)+slopesFC_X(i,j-1,k,BX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j-1,k,BX_LOCAL)-slopesFC_X(i,j-1,k,BX_LOCAL));
		  az = 0.0;
		  axz = 0.0;
		  bx = (slopesFC_Y(i,j,k,BY_LOCAL)+slopesFC_Y(i,j-1,k,BY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i,j,k,BY_LOCAL)-slopesFC_Y(i,j-1,k,BY_LOCAL));
		  bz = 0.0;
		  byz = 0.0;
		  cx = slopesCC_X(i,j-1,k,BZ_LOCAL);
		  cxz = 0.0;
		  cy = slopesCC_Y(i,j-1,k,BZ_LOCAL);
		  cyz = 0.0;
		  axx = -dx[0]/2.0*(bxy/dx[1]); // +cxz/dx[2]
		  byy = -dx[1]/2.0*(axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j-1,k,BX_LOCAL)+arrEM_X(i,j-1,k,BX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i+1,j-1,k,BX_LOCAL)-arrEM_X(i,j-1,k,BX_LOCAL));
		  b0 = (arrEM_Y(i,j,k,BY_LOCAL)+arrEM_Y(i,j-1,k,BY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j,k,BY_LOCAL)-arrEM_Y(i,j-1,k,BY_LOCAL));
		  c0 = arr(i,j-1,k,BZ) - czz/6.0; //czz*dx[2]*dx[2]/6.0
		  cz = 0.0;
		  Real BxRD = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real ByRD = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real BzRD = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

		  // LU state
		  x = dx[0]/2.0, y = -dx[1]/2.0;
		  ay = (slopesFC_X(i,j,k,BX_LOCAL)+slopesFC_X(i-1,j,k,BX_LOCAL))/2.0;
		  axy = (slopesFC_X(i,j,k,BX_LOCAL)-slopesFC_X(i-1,j,k,BX_LOCAL));
		  az = 0.0;
		  axz = 0.0;
		  bx = (slopesFC_Y(i-1,j+1,k,BY_LOCAL)+slopesFC_Y(i-1,j,k,BY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i-1,j+1,k,BY_LOCAL)-slopesFC_Y(i-1,j,k,BY_LOCAL));
		  bz = 0.0;
		  byz = 0.0;
		  cx = slopesCC_X(i-1,j,k,BZ_LOCAL);
		  cxz = 0.0;
		  cy = slopesCC_Y(i-1,j,k,BZ_LOCAL);
		  cyz = 0.0;
		  axx = -dx[0]/2.0*(bxy/dx[1]); // +cxz/dx[2]
		  byy = -dx[1]/2.0*(axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i,j,k,BX_LOCAL)+arrEM_X(i-1,j,k,BX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i,j,k,BX_LOCAL)-arrEM_X(i-1,j,k,BX_LOCAL));
		  b0 = (arrEM_Y(i-1,j+1,k,BY_LOCAL)+arrEM_Y(i-1,j,k,BY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i-1,j+1,k,BY_LOCAL)-arrEM_Y(i-1,j,k,BY_LOCAL));
		  c0 = arr(i-1,j,k,BZ) - czz/6.0; //czz*dx[2]*dx[2]/6.0
		  cz = 0.0;
		  Real BxLU = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real ByLU = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real BzLU = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

		  // RU state
		  x = -dx[0]/2.0, y = -dx[1]/2.0;
		  ay = (slopesFC_X(i+1,j,k,BX_LOCAL)+slopesFC_X(i,j,k,BX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j,k,BX_LOCAL)-slopesFC_X(i,j,k,BX_LOCAL));
		  az = 0.0;
		  axz = 0.0;
		  bx = (slopesFC_Y(i,j+1,k,BY_LOCAL)+slopesFC_Y(i,j,k,BY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i,j+1,k,BY_LOCAL)-slopesFC_Y(i,j,k,BY_LOCAL));
		  bz = 0.0;
		  byz = 0.0;
		  cx = slopesCC_X(i,j,k,BZ_LOCAL);
		  cxz = 0.0;
		  cy = slopesCC_Y(i,j,k,BZ_LOCAL);
		  cyz = 0.0;
		  axx = -dx[0]/2.0*(bxy/dx[1]); // +cxz/dx[2]
		  byy = -dx[1]/2.0*(axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i+1,j,k,BX_LOCAL)-arrEM_X(i,j,k,BX_LOCAL));
		  b0 = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j+1,k,BY_LOCAL)-arrEM_Y(i,j,k,BY_LOCAL));
		  c0 = arr(i,j,k,BZ) - czz/6.0; //czz*dx[2]*dx[2]/6.0
		  cz = 0.0;
		  Real BxRU = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real ByRU = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real BzRU = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

		  // Electric field
		  // LD state
		  x = dx[0]/2.0, y = dx[1]/2.0;
		  ay = (slopesFC_X(i,j-1,k,EX_LOCAL)+slopesFC_X(i-1,j-1,k,EX_LOCAL))/2.0;
		  axy = (slopesFC_X(i,j-1,k,EX_LOCAL)-slopesFC_X(i-1,j-1,k,EX_LOCAL));
		  az = 0.0;
		  axz = 0.0;
		  bx = (slopesFC_Y(i-1,j,k,EY_LOCAL)+slopesFC_Y(i-1,j-1,k,EY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i-1,j,k,EY_LOCAL)-slopesFC_Y(i-1,j-1,k,EY_LOCAL));
		  bz = 0.0;
		  byz = 0.0;
		  cx = slopesCC_X(i-1,j-1,k,EZ_LOCAL);
		  cxz = 0.0;
		  cy = slopesCC_Y(i-1,j-1,k,EZ_LOCAL);
		  cyz = 0.0;
		  axx = dx[0]/2.0*(slopesQ(i-1,j-1,k,0)-bxy/dx[1]); // +cxz/dx[2]
		  byy = dx[1]/2.0*(slopesQ(i-1,j-1,k,1)-axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i,j-1,k,EX_LOCAL)+arrEM_X(i-1,j-1,k,EX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i,j-1,k,EX_LOCAL)-arrEM_X(i-1,j-1,k,EX_LOCAL));
		  b0 = (arrEM_Y(i-1,j,k,EY_LOCAL)+arrEM_Y(i-1,j-1,k,EY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i-1,j,k,EY_LOCAL)-arrEM_Y(i-1,j-1,k,EY_LOCAL));
		  c0 = arr(i-1,j-1,k,EZ) - czz/6.0; //czz*dx[2]*dx[2]/6.0
		  cz = 0.0;
		  Real ExLD = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real EyLD = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real EzLD = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

		  // RD state
		  x = -dx[0]/2.0, y = dx[1]/2.0;
		  ay = (slopesFC_X(i+1,j-1,k,EX_LOCAL)+slopesFC_X(i,j-1,k,EX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j-1,k,EX_LOCAL)-slopesFC_X(i,j-1,k,EX_LOCAL));
		  az = 0.0;
		  axz = 0.0;
		  bx = (slopesFC_Y(i,j,k,EY_LOCAL)+slopesFC_Y(i,j-1,k,EY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i,j,k,EY_LOCAL)-slopesFC_Y(i,j-1,k,EY_LOCAL));
		  bz = 0.0;
		  byz = 0.0;
		  cx = slopesCC_X(i,j-1,k,EZ_LOCAL);
		  cxz = 0.0;
		  cy = slopesCC_Y(i,j-1,k,EZ_LOCAL);
		  cyz = 0.0;
		  axx = dx[0]/2.0*(slopesQ(i,j-1,k,0)-bxy/dx[1]); // +cxz/dx[2]
		  byy = dx[1]/2.0*(slopesQ(i,j-1,k,1)-axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j-1,k,EX_LOCAL)+arrEM_X(i,j-1,k,EX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i+1,j-1,k,EX_LOCAL)-arrEM_X(i,j-1,k,EX_LOCAL));
		  b0 = (arrEM_Y(i,j,k,EY_LOCAL)+arrEM_Y(i,j-1,k,EY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j,k,EY_LOCAL)-arrEM_Y(i,j-1,k,EY_LOCAL));
		  c0 = arr(i,j-1,k,EZ) - czz/6.0;
		  cz = 0.0;
		  Real ExRD = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real EyRD = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real EzRD = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

		  // LU state
		  x = dx[0]/2.0, y = -dx[1]/2.0;
		  ay = (slopesFC_X(i,j,k,EX_LOCAL)+slopesFC_X(i-1,j,k,EX_LOCAL))/2.0;
		  axy = (slopesFC_X(i,j,k,EX_LOCAL)-slopesFC_X(i-1,j,k,EX_LOCAL));
		  az = 0.0;
		  axz = 0.0;
		  bx = (slopesFC_Y(i-1,j+1,k,EY_LOCAL)+slopesFC_Y(i-1,j,k,EY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i-1,j+1,k,EY_LOCAL)-slopesFC_Y(i-1,j,k,EY_LOCAL));
		  bz = 0.0;
		  byz = 0.0;
		  cx = slopesCC_X(i-1,j,k,EZ_LOCAL);
		  cxz = 0.0;
		  cy = slopesCC_Y(i-1,j,k,EZ_LOCAL);
		  cyz = 0.0;
		  axx = dx[0]/2.0*(slopesQ(i-1,j,k,0)-bxy/dx[1]); // +cxz/dx[2]
		  byy = dx[1]/2.0*(slopesQ(i-1,j,k,1)-axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i,j,k,EX_LOCAL)+arrEM_X(i-1,j,k,EX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i,j,k,EX_LOCAL)-arrEM_X(i-1,j,k,EX_LOCAL));
		  b0 = (arrEM_Y(i-1,j+1,k,EY_LOCAL)+arrEM_Y(i-1,j,k,EY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i-1,j+1,k,EY_LOCAL)-arrEM_Y(i-1,j,k,EY_LOCAL));
		  c0 = arr(i-1,j,k,EZ) - czz/6.0; //czz*dx[2]*dx[2]/6.0
		  cz = 0.0;
		  Real ExLU = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real EyLU = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real EzLU = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

		  // RU state
		  x = -dx[0]/2.0, y = -dx[1]/2.0;
		  ay = (slopesFC_X(i+1,j,k,EX_LOCAL)+slopesFC_X(i,j,k,EX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j,k,EX_LOCAL)-slopesFC_X(i,j,k,EX_LOCAL));
		  az = 0.0;
		  axz = 0.0;
		  bx = (slopesFC_Y(i,j+1,k,EY_LOCAL)+slopesFC_Y(i,j,k,EY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i,j+1,k,EY_LOCAL)-slopesFC_Y(i,j,k,EY_LOCAL));
		  bz = 0.0;
		  byz = 0.0;
		  cx = slopesCC_X(i,j,k,EZ_LOCAL);
		  cxz = 0.0;
		  cy = slopesCC_Y(i,j,k,EZ_LOCAL);
		  cyz = 0.0;
		  axx = dx[0]/2.0*(slopesQ(i,j,k,0)-bxy/dx[1]); // +cxz/dx[2]
		  byy = dx[1]/2.0*(slopesQ(i,j,k,1)-axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - axx/6.0;	
		  ax = (arrEM_X(i+1,j,k,EX_LOCAL)-arrEM_X(i,j,k,EX_LOCAL));
		  b0 = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j+1,k,EY_LOCAL)-arrEM_Y(i,j,k,EY_LOCAL));
		  c0 = arr(i,j,k,EZ) - czz/6.0; //czz*dx[2]*dx[2]/6.0
		  cz = 0.0;
		  Real ExRU = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real EyRU = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real EzRU = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  
		  
		  // If the reconstruction is consistent, then EyRU=EyRD
		  Real EyR = (EyRU+EyRD)/2.0;
		  Real EyL = (EyLU+EyLD)/2.0;
		  Real ExU = (ExRU+ExLU)/2.0;
		  Real ExD = (ExRD+ExLD)/2.0;
		  fluxArrEM(i,j,k,BZ_LOCAL) = 0.25*(BzLD+BzRD+BzLU+BzRU) - 0.5/c*(EyR-EyL) + 0.5/c*(ExU-ExD);

		  Real ByR = (ByRU+ByRD)/2.0;
		  Real ByL = (ByLU+ByLD)/2.0;
		  Real BxU = (BxRU+BxLU)/2.0;
		  Real BxD = (BxRD+BxLD)/2.0;
		  fluxArrEM(i,j,k,EZ_LOCAL) = 0.25*(EzLD+EzRD+EzLU+EzRU) + 0.5*c*(ByR-ByL) - 0.5*c*(BxU-BxD);

		  if (i<3 && j==5)
		    std::cout << "Old " << i << " "  << fluxArrEM(i,j,k,BZ_LOCAL) << " " << fluxArrEM(i,j,k,EZ_LOCAL)  << std::endl;

		  
		}
	    }
	}       
    }

  // Update face-centred EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
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
	  //const auto& arr = Sborder.array(mfi);
	  const auto& arrEM = S_EM[d_EM].array(mfi);
	  //const auto& arrEMother = S_EM[d].array(mfi);
	  const auto& fluxArrEM = fluxesEM.array(mfi);
	  const auto& fluxArr = fluxes[d_EM].array(mfi);

	  const Dim3 hiDomain = ubound(geom.Domain());
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      // 2D code; z-component updated using cell-centred scheme
		      arrEM(i,j,k,BX_LOCAL+d_EM) += std::pow(-1,d)*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EZ_LOCAL) - fluxArrEM(i,j,k,EZ_LOCAL));
		      arrEM(i,j,k,EX_LOCAL+d_EM) -= std::pow(-1,d)*c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) - fluxArrEM(i,j,k,BZ_LOCAL));

		      // 3D code
		      /*arrEM(i,j,k,BX_LOCAL+(1+d)%3) += (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(2+d)%3));
			arrEM(i,j,k,BX_LOCAL+(2+d)%3) -= (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(1+d)%3));
			arrEM(i,j,k,EX_LOCAL+(1+d)%3) -= c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(2+d)%3));
			arrEM(i,j,k,EX_LOCAL+(2+d)%3) += c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(1+d)%3));
		      */
		      
		      // source terms
		      Real currentFace = r_i*fluxArr(i,j,k,RHO_I) + r_e*fluxArr(i,j,k,RHO_E);
		      arrEM(i,j,k,EX_LOCAL+d_EM) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;

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

  // Compute cell-centred EM fields from face-centred
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi); 
      //const auto& fluxArrEM = fluxesEM.array(mfi);
      //const auto& fluxArrX = fluxes[0].array(mfi);
      //const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& slopesFC_X = slopesFC[0].array(mfi);
      const auto& slopesFC_Y = slopesFC[1].array(mfi);

      // only for 2D
      const auto& slopesCC_X = slopesCC[0].array(mfi);
      const auto& slopesCC_Y = slopesCC[1].array(mfi);
  
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 
		  
  		  Real a0,axx,axy,axz,b0,bxy,byy,byz,c0,cxz,cyz,czz;
  		  Real x,y;

  		  // Magnetic field       
		  axy = (slopesFC_X(i+1,j,k,BX_LOCAL)-slopesFC_X(i,j,k,BX_LOCAL));
  		  axz = 0.0;
		  bxy = (slopesFC_Y(i,j+1,k,BY_LOCAL)-slopesFC_Y(i,j,k,BY_LOCAL));
  		  byz = 0.0;
  		  cxz = 0.0;
  		  cyz = 0.0;
		  axx = -dx[0]/2.0*(bxy/dx[1]); // +cxz/dx[2]
		  byy = -dx[1]/2.0*(axy/dx[0]); // ++cyz/dx[2]
  		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - axx/6.0;
		  b0 = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - byy/6.0;
  		  c0 = arr(i,j,k,BZ) - czz/6.0;

  		  arr(i,j,k,BX) = a0;
  		  arr(i,j,k,BY) = b0;
		  
  		  // Electric field
  		  axy = (slopesFC_X(i+1,j,k,EX_LOCAL)-slopesFC_X(i,j,k,EX_LOCAL));
  		  axz = 0.0;
  		  bxy = (slopesFC_Y(i,j+1,k,EY_LOCAL)-slopesFC_Y(i,j,k,EY_LOCAL));
  		  byz = 0.0;
  		  cxz = 0.0;
  		  cyz = 0.0;
		  axx = dx[0]/2.0*(slopesQ(i,j,k,0)-bxy/dx[1]); // +cxz/dx[2]
		  byy = dx[1]/2.0*(slopesQ(i,j,k,1)-axy/dx[0]); // ++cyz/dx[2]
  		  czz = 0.0;
  		  a0 = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - axx/6.0;
  		  b0 = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - byy/6.0;
  		  c0 = arr(i,j,k,EZ) - czz/6.0;

  		  arr(i,j,k,EX) = a0;
  		  arr(i,j,k,EY) = b0;
		  
  		}
  	    }
  	}       
    }

  // Only 2D
  // Compute edge components for y-components at x-faces
  for (MFIter mfi(fluxes[0], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // use old data
      //const auto& arr = Sborder.array(mfi);
      const auto& arr = S0.array(mfi);
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi); 
      //const auto& fluxArrEM = fluxesEM.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      //const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& slopesFC_X = slopesFC[0].array(mfi);
      const auto& slopesFC_Y = slopesFC[1].array(mfi);

      // only for 2D
      const auto& slopesCC_X = slopesCC[0].array(mfi);
      const auto& slopesCC_Y = slopesCC[1].array(mfi);
  
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 
		  
  		  Real a0,ax,ay,az,axx,axy,axz,b0,bx,by,bz,bxy,byy,byz,c0,cx,cy,cz,cxz,cyz,czz;
  		  Real x,y;
		  
  		  // Magnetic field
  		  // L state
  		  x = dx[0]/2.0, y = 0.0;
  		  ay = (slopesFC_X(i,j,k,BX_LOCAL)+slopesFC_X(i-1,j,k,BX_LOCAL))/2.0;
		  axy = (slopesFC_X(i,j,k,BX_LOCAL)-slopesFC_X(i-1,j,k,BX_LOCAL));
  		  az = 0.0;
  		  axz = 0.0;
  		  bx = (slopesFC_Y(i-1,j+1,k,BY_LOCAL)+slopesFC_Y(i-1,j,k,BY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i-1,j+1,k,BY_LOCAL)-slopesFC_Y(i-1,j,k,BY_LOCAL));
  		  bz = 0.0;
  		  byz = 0.0;
  		  cx = slopesCC_X(i-1,j,k,BZ_LOCAL);
  		  cxz = 0.0;
  		  cy = slopesCC_Y(i-1,j,k,BZ_LOCAL);
  		  cyz = 0.0;
		  axx = -dx[0]/2.0*(bxy/dx[1]); // +cxz/dx[2]
		  byy = -dx[1]/2.0*(axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i,j,k,BX_LOCAL)+arrEM_X(i-1,j,k,BX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i,j,k,BX_LOCAL)-arrEM_X(i-1,j,k,BX_LOCAL));
		  b0 = (arrEM_Y(i-1,j+1,k,BY_LOCAL)+arrEM_Y(i-1,j,k,BY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i-1,j+1,k,BY_LOCAL)-arrEM_Y(i-1,j,k,BY_LOCAL));
  		  c0 = arr(i-1,j,k,BZ) - czz/6.0;
  		  cz = 0.0;
		  Real BxL = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real ByL = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real BzL = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

  		  // R state
  		  x = -dx[0]/2.0, y = 0.0;
  		  ay = (slopesFC_X(i+1,j,k,BX_LOCAL)+slopesFC_X(i,j,k,BX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j,k,BX_LOCAL)-slopesFC_X(i,j,k,BX_LOCAL));
  		  az = 0.0;
  		  axz = 0.0;
  		  bx = (slopesFC_Y(i,j+1,k,BY_LOCAL)+slopesFC_Y(i,j,k,BY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i,j+1,k,BY_LOCAL)-slopesFC_Y(i,j,k,BY_LOCAL));
  		  bz = 0.0;
  		  byz = 0.0;
  		  cx = slopesCC_X(i,j,k,BZ_LOCAL);
  		  cxz = 0.0;
  		  cy = slopesCC_Y(i,j,k,BZ_LOCAL);
  		  cyz = 0.0;
		  axx = -dx[0]/2.0*(bxy/dx[1]); // +cxz/dx[2]
		  byy = -dx[1]/2.0*(axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i+1,j,k,BX_LOCAL)-arrEM_X(i,j,k,BX_LOCAL));
		  b0 = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j+1,k,BY_LOCAL)-arrEM_Y(i,j,k,BY_LOCAL));
  		  c0 = arr(i,j,k,BZ) - czz/6.0;
  		  cz = 0.0;
		  Real BxR = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real ByR = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real BzR = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

  		  // Electric field
  		  // L state
  		  x = dx[0]/2.0, y = 0.0;
  		  ay = (slopesFC_X(i,j,k,EX_LOCAL)+slopesFC_X(i-1,j,k,EX_LOCAL))/2.0;
		  axy = (slopesFC_X(i,j,k,EX_LOCAL)-slopesFC_X(i-1,j,k,EX_LOCAL));
  		  az = 0.0;
  		  axz = 0.0;
  		  bx = (slopesFC_Y(i-1,j+1,k,EY_LOCAL)+slopesFC_Y(i-1,j,k,EY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i-1,j+1,k,EY_LOCAL)-slopesFC_Y(i-1,j,k,EY_LOCAL));
  		  bz = 0.0;
  		  byz = 0.0;
  		  cx = slopesCC_X(i-1,j,k,EZ_LOCAL);
  		  cxz = 0.0;
  		  cy = slopesCC_Y(i-1,j,k,EZ_LOCAL);
  		  cyz = 0.0;
		  axx = dx[0]/2.0*(slopesQ(i-1,j,k,0)-bxy/dx[1]); // +cxz/dx[2]  
		  byy = dx[1]/2.0*(slopesQ(i-1,j,k,1)-axy/dx[0]); // ++cyz/dx[2]
  		  czz = 0.0;
		  a0 = (arrEM_X(i,j,k,EX_LOCAL)+arrEM_X(i-1,j,k,EX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i,j,k,EX_LOCAL)-arrEM_X(i-1,j,k,EX_LOCAL));
		  b0 = (arrEM_Y(i-1,j+1,k,EY_LOCAL)+arrEM_Y(i-1,j,k,EY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i-1,j+1,k,EY_LOCAL)-arrEM_Y(i-1,j,k,EY_LOCAL));
  		  c0 = arr(i-1,j,k,EZ) - czz/6.0;
  		  cz = 0.0;
		  Real ExL = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real EyL = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real EzL = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

  		  // R state
  		  x = -dx[0]/2.0, y = 0.0;
  		  ay = (slopesFC_X(i+1,j,k,EX_LOCAL)+slopesFC_X(i,j,k,EX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j,k,EX_LOCAL)-slopesFC_X(i,j,k,EX_LOCAL));
  		  az = 0.0;
  		  axz = 0.0;
  		  bx = (slopesFC_Y(i,j+1,k,EY_LOCAL)+slopesFC_Y(i,j,k,EY_LOCAL))/2.0;  	
		  bxy = (slopesFC_Y(i,j+1,k,EY_LOCAL)-slopesFC_Y(i,j,k,EY_LOCAL));
  		  bz = 0.0;
  		  byz = 0.0;
  		  cx = slopesCC_X(i,j,k,EZ_LOCAL);
  		  cxz = 0.0;
  		  cy = slopesCC_Y(i,j,k,EZ_LOCAL);
  		  cyz = 0.0;
		  axx = dx[0]/2.0*(slopesQ(i,j,k,0)-bxy/dx[1]); // +cxz/dx[2]
		  byy = dx[1]/2.0*(slopesQ(i,j,k,1)-axy/dx[0]); // ++cyz/dx[2]
  		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i+1,j,k,EX_LOCAL)-arrEM_X(i,j,k,EX_LOCAL));
		  b0 = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j+1,k,EY_LOCAL)-arrEM_Y(i,j,k,EY_LOCAL));
  		  c0 = arr(i,j,k,EZ) - czz/6.0;
  		  cz = 0.0;
		  Real ExR = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real EyR = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real EzR = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

  		  fluxArrX(i,j,k,BY) = 0.5*(ByL+ByR) + 0.5/c*(EzR-EzL);
  		  fluxArrX(i,j,k,EY) = 0.5*(EyL+EyR) - 0.5*c*(BzR-BzL);
  		}
  	    }
  	}
    }
  
  // Only 2D
  // Compute edge components for x-components at y-faces
  for (MFIter mfi(fluxes[1], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      // use old data
      //const auto& arr = Sborder.array(mfi);
      const auto& arr = S0.array(mfi);
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi); 
      //const auto& fluxArrEM = fluxesEM.array(mfi);
      //const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& slopesFC_X = slopesFC[0].array(mfi);
      const auto& slopesFC_Y = slopesFC[1].array(mfi);

      // only for 2D
      const auto& slopesCC_X = slopesCC[0].array(mfi);
      const auto& slopesCC_Y = slopesCC[1].array(mfi);
  
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 
		  
  		  Real a0,ax,ay,az,axx,axy,axz,b0,bx,by,bz,bxy,byy,byz,c0,cx,cy,cz,cxz,cyz,czz;
  		  Real x,y;
		  	
  		  // D state
  		  x = 0.0, y = dx[1]/2.0;
  		  ay = (slopesFC_X(i+1,j-1,k,BX_LOCAL)+slopesFC_X(i,j-1,k,BX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j-1,k,BX_LOCAL)-slopesFC_X(i,j-1,k,BX_LOCAL));
  		  az = 0.0;
  		  axz = 0.0;
  		  bx = (slopesFC_Y(i,j,k,BY_LOCAL)+slopesFC_Y(i,j-1,k,BY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i,j,k,BY_LOCAL)-slopesFC_Y(i,j-1,k,BY_LOCAL));
  		  bz = 0.0;
  		  byz = 0.0;
  		  cx = slopesCC_X(i,j-1,k,BZ_LOCAL);
  		  cxz = 0.0;
  		  cy = slopesCC_Y(i,j-1,k,BZ_LOCAL);
  		  cyz = 0.0;
		  axx = -dx[0]/2.0*(bxy/dx[1]); // +cxz/dx[2]
		  byy = -dx[1]/2.0*(axy/dx[0]); // ++cyz/dx[2]
  		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j-1,k,BX_LOCAL)+arrEM_X(i,j-1,k,BX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i+1,j-1,k,BX_LOCAL)-arrEM_X(i,j-1,k,BX_LOCAL));
		  b0 = (arrEM_Y(i,j,k,BY_LOCAL)+arrEM_Y(i,j-1,k,BY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j,k,BY_LOCAL)-arrEM_Y(i,j-1,k,BY_LOCAL));
  		  c0 = arr(i,j-1,k,BZ) - czz/6.0;
  		  cz = 0.0;
		  Real BxD = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real ByD = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real BzD = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  
		  
  		  // U state
  		  x = 0.0, y = -dx[1]/2.0;
  		  ay = (slopesFC_X(i+1,j,k,BX_LOCAL)+slopesFC_X(i,j,k,BX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j,k,BX_LOCAL)-slopesFC_X(i,j,k,BX_LOCAL));
  		  az = 0.0;
  		  axz = 0.0;
  		  bx = (slopesFC_Y(i,j+1,k,BY_LOCAL)+slopesFC_Y(i,j,k,BY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i,j+1,k,BY_LOCAL)-slopesFC_Y(i,j,k,BY_LOCAL));
  		  bz = 0.0;
  		  byz = 0.0;
  		  cx = slopesCC_X(i,j,k,BZ_LOCAL);
  		  cxz = 0.0;
  		  cy = slopesCC_Y(i,j,k,BZ_LOCAL);
  		  cyz = 0.0;
		  axx = -dx[0]/2.0*(bxy/dx[1]); // +cxz/dx[2]
		  byy = -dx[1]/2.0*(axy/dx[0]); // ++cyz/dx[2]
  		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i+1,j,k,BX_LOCAL)-arrEM_X(i,j,k,BX_LOCAL));
		  b0 = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j+1,k,BY_LOCAL)-arrEM_Y(i,j,k,BY_LOCAL));
  		  c0 = arr(i,j,k,BZ) - czz/6.0;
  		  cz = 0.0;
		  Real BxU = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real ByU = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real BzU = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  
		  
  		  // Electric field
  		  // D state
  		  x = 0.0, y = dx[1]/2.0;
  		  ay = (slopesFC_X(i+1,j-1,k,EX_LOCAL)+slopesFC_X(i,j-1,k,EX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j-1,k,EX_LOCAL)-slopesFC_X(i,j-1,k,EX_LOCAL));
  		  az = 0.0;
  		  axz = 0.0;
  		  bx = (slopesFC_Y(i,j,k,EY_LOCAL)+slopesFC_Y(i,j-1,k,EY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i,j,k,EY_LOCAL)-slopesFC_Y(i,j-1,k,EY_LOCAL));
  		  bz = 0.0;
  		  byz = 0.0;
  		  cx = slopesCC_X(i,j-1,k,EZ_LOCAL);
  		  cxz = 0.0;
  		  cy = slopesCC_Y(i,j-1,k,EZ_LOCAL);
  		  cyz = 0.0;
		  axx = dx[0]/2.0*(slopesQ(i,j-1,k,0)-bxy/dx[1]); // +cxz/dx[2]
		  byy = dx[1]/2.0*(slopesQ(i,j-1,k,1)-axy/dx[0]); // ++cyz/dx[2]
		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j-1,k,EX_LOCAL)+arrEM_X(i,j-1,k,EX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i+1,j-1,k,EX_LOCAL)-arrEM_X(i,j-1,k,EX_LOCAL));
		  b0 = (arrEM_Y(i,j,k,EY_LOCAL)+arrEM_Y(i,j-1,k,EY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j,k,EY_LOCAL)-arrEM_Y(i,j-1,k,EY_LOCAL));
  		  c0 = arr(i,j-1,k,EZ) - czz/6.0;
  		  cz = 0.0;
		  Real ExD = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real EyD = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real EzD = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

  		  // U state
  		  x = 0.0, y = -dx[1]/2.0;
  		  ay = (slopesFC_X(i+1,j,k,EX_LOCAL)+slopesFC_X(i,j,k,EX_LOCAL))/2.0;
		  axy = (slopesFC_X(i+1,j,k,EX_LOCAL)-slopesFC_X(i,j,k,EX_LOCAL));
  		  az = 0.0;
  		  axz = 0.0;
  		  bx = (slopesFC_Y(i,j+1,k,EY_LOCAL)+slopesFC_Y(i,j,k,EY_LOCAL))/2.0;
		  bxy = (slopesFC_Y(i,j+1,k,EY_LOCAL)-slopesFC_Y(i,j,k,EY_LOCAL));
  		  bz = 0.0;
  		  byz = 0.0;
  		  cx = slopesCC_X(i,j,k,EZ_LOCAL);
  		  cxz = 0.0;
  		  cy = slopesCC_Y(i,j,k,EZ_LOCAL);
  		  cyz = 0.0;
		  axx = dx[0]/2.0*(slopesQ(i,j,k,0)-bxy/dx[1]); // +cxz/dx[2]
		  byy = dx[1]/2.0*(slopesQ(i,j,k,1)-axy/dx[0]); // ++cyz/dx[2]		  
  		  czz = 0.0;
		  a0 = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - axx/6.0;
		  ax = (arrEM_X(i+1,j,k,EX_LOCAL)-arrEM_X(i,j,k,EX_LOCAL));
		  b0 = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - byy/6.0;
		  by = (arrEM_Y(i,j+1,k,EY_LOCAL)-arrEM_Y(i,j,k,EY_LOCAL));
  		  c0 = arr(i,j,k,EZ) - czz/6.0;
  		  cz = 0.0;
		  Real ExU = a0 + ax*(x/dx[0]) + ay*(y/dx[1]) // + az*
		    + axx*((x/dx[0])*(x/dx[0]) - 1.0/12.0) + axy*(x/dx[0])*(y/dx[1]); // + axz*
		  Real EyU = b0 + bx*(x/dx[0]) + by*(y/dx[1]) // + bz*
		    + byy*((y/dx[1])*(y/dx[1]) - 1.0/12.0) + bxy*(x/dx[0])*(y/dx[1]); // + bxz*
		  Real EzU = c0 + cx*(x/dx[0]) + cy*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

  		  fluxArrY(i,j,k,BX) = 0.5*(BxD+BxU) - 0.5/c*(EzU-EzD);
  		  fluxArrY(i,j,k,EX) = 0.5*(ExD+ExU) + 0.5*c*(BzU-BzD); 
  		}
  	    }
  	}
    }

  // Update cell-centred z-components of EM fields
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);
	        
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  // Update cell-centred z-components becuause it is 2D code
  		  arr(i,j,k,BZ) = arr(i,j,k,BZ) - (dt / dx[0]) * (fluxArrX(i+1,j,k,EY) - fluxArrX(i,j,k,EY)) + (dt / dx[1]) * (fluxArrY(i,j+1,k,EX) - fluxArrY(i,j,k,EX));
  		  arr(i,j,k,EZ) = arr(i,j,k,EZ) + c*c*(dt / dx[0]) * (fluxArrX(i+1,j,k,BY) - fluxArrX(i,j,k,BY)) - c*c*(dt / dx[1]) * (fluxArrY(i,j+1,k,BX) - fluxArrY(i,j,k,BX));		  

  		  // source terms
  		  arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
  		}
  	    }
  	}
    }
    
}
void CAMReXmp::MaxwellSolverDivFree(Array<MultiFab,AMREX_SPACEDIM>& S_EM, MultiFab& fluxesEM, MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], MultiFab& S0, const Real* dx, Real dt)
{

  // Note that fluxesEM and fluxes are different definitions
  // fluxesEM are EM star states, fluxes are fluxes
  
  // Useful when dealing with high-order reconstruction
  // To get first-order slope in y-direction for Bx: slopes[BX_LOCALL].array(i,j,k,Y)
  // Useful to define indices int X,Y,Z in second-order linear reconstruction
  // For third-order quadratic reconstruction use int X,Y,Z,XX,YY,ZZ,XY,YZ,XZ
  // Although for x-components fields, only exist Y,Z,YY,ZZ,YZ (see Balsara et al. 2017)
  
  Array<MultiFab,6> slopes;
  // number of slopes
  // 2 slopes for second-order linear reconstruction
  // Ex = Ex^0 + Ex^y*(y/Delta_y) + Ex^z*(z/Delta_z), where Ex^y and Ex^z are the slopes
  // For 2D code, only 1 slope for x- and y-compoenents
  // For clarity, will use 3 slopes for second-order, first one represents slope in x-direction
  // 9 slopes for third-order, and so on
  int nSlopes = 3;
  // indices for the slopes
  // for 2D, do not need Z index
  // for third order, use also XX,YY,ZZ,XY,YZ,XZ
  int X=0,Y=1,Z=2;
  slopes[BX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);
  slopes[BY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
  slopes[EX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);

  // For 2D code, 2 slopes for z-components, although it is cell-centred
  slopes[BZ_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EZ_LOCAL].define(grids, dmap, nSlopes, 1);
    
  // Three components, one for each direction
  MultiFab slopesCharge(grids, dmap, AMREX_SPACEDIM, 1);

  // Compute slopes of the charge densities
  for (int d = 0; d < AMREX_SPACEDIM ; d++)   
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      // Compute cell-centred density slopes
      for (MFIter mfi(slopesCharge, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  // uses old charge
	  const Array4<Real> arr = S0.array(mfi);
	  Array4<Real> slopes = slopesCharge.array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y-1; j <= hi.y+1; j++)
		{
		  for(int i = lo.x-1; i <= hi.x+1; i++)
		    {
		      Real u_iMinus1 = (r_i*arr(i-iOffset,j-jOffset,k-kOffset,RHO_I)
					+ r_e*arr(i-iOffset,j-jOffset,k-kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_iPlus1 = (r_i*arr(i+iOffset,j+jOffset,k+kOffset,RHO_I)
				       + r_e*arr(i+iOffset,j+jOffset,k+kOffset,RHO_E))/(lambda_d*lambda_d*l_r);

		      // slopes for the current
		      slopes(i,j,k,d) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		    }
		}
	    }      
	}
    }
  
  // Compute cell-centred z-components slopes in x- and y-direction
  for (int d = 0; d < AMREX_SPACEDIM ; d++)
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
      
      for (MFIter mfi(slopes[BZ_LOCAL], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Array4<Real> arr = Sborder.array(mfi);
	  Array4<Real> slopesBZ = slopes[BZ_LOCAL].array(mfi);
	  Array4<Real> slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y-1; j <= hi.y+1; j++)
		{
		  for(int i = lo.x-1; i <= hi.x+1; i++)
		    {

		      Real u_iMinus1 = arr(i-iOffset,j-jOffset,k-kOffset,BZ);
		      Real u_i = arr(i,j,k,BZ);
		      Real u_iPlus1 = arr(i+iOffset,j+jOffset,k+kOffset,BZ);
		      slopesBZ(i,j,k,d) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		      u_iMinus1 = arr(i-iOffset,j-jOffset,k-kOffset,EZ);
		      u_i = arr(i,j,k,EZ);
		      u_iPlus1 = arr(i+iOffset,j+jOffset,k+kOffset,EZ);
		      slopesEZ(i,j,k,d) = TVD_slope(u_iMinus1,u_i,u_iPlus1);		    

		    }
		}
	    }      
	}
    }
  
  // Compute y-components slopes in x-direction
  for (MFIter mfi(slopes[BY_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Array4<Real> arr = S_EM[1].array(mfi);
      Array4<Real> slopesBY = slopes[BY_LOCAL].array(mfi);
      Array4<Real> slopesEY = slopes[EY_LOCAL].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{

		  Real u_iMinus1 = arr(i-1,j,k,BY_LOCAL);
		  Real u_i = arr(i,j,k,BY_LOCAL);		    
		  Real u_iPlus1 = arr(i+1,j,k,BY_LOCAL);
		  slopesBY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		  u_iMinus1 = arr(i-1,j,k,EY_LOCAL);
		  u_i = arr(i,j,k,EY_LOCAL);		    
		  u_iPlus1 = arr(i+1,j,k,EY_LOCAL);
		  slopesEY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		}
	    }
	}      
    }
  // Compute x-components slopes in y-direction
  for (MFIter mfi(slopes[BX_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Array4<Real> arr = S_EM[0].array(mfi);
      Array4<Real> slopesBX = slopes[BX_LOCAL].array(mfi);
      Array4<Real> slopesEX = slopes[EX_LOCAL].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{

		  Real u_iMinus1 = arr(i,j-1,k,BX_LOCAL);
		  Real u_i = arr(i,j,k,BX_LOCAL);		    
		  Real u_iPlus1 = arr(i,j+1,k,BX_LOCAL);
		  slopesBX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		  u_iMinus1 = arr(i,j-1,k,EX_LOCAL);
		  u_i = arr(i,j,k,EX_LOCAL);		    
		  u_iPlus1 = arr(i,j+1,k,EX_LOCAL);
		  slopesEX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		}
	    }
	}      
    }


  // Coefficients for the magnetic and electric fields
  // For second order, there are 21 coeff.
  // i.e. a0,ax,ay,az,axx,...,b0,bx,...,c0,cx,...,czz
  // use 1 ghost cell
  int a0=0,ax=1,ay=2,az=3,axx=4,axy=5,axz=6;
  int b0=7,bx=8,by=9,bz=10,byy=11,bxy=12,byz=13;
  int c0=14,cx=15,cy=16,cz=17,czz=18,cxz=19,cyz=20;
  
  MultiFab Bcoeff(grids, dmap, 21, 1);
  MultiFab Ecoeff(grids, dmap, 21, 1);

  // Compute coefficients for the field function
  for (MFIter mfi(Bcoeff, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      // coeff
      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      // data
      const auto& arr = S0.array(mfi);
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi);       

      // slopes
      const auto& slopesBX = slopes[BX_LOCAL].array(mfi);
      const auto& slopesBY = slopes[BY_LOCAL].array(mfi);
      const auto& slopesBZ = slopes[BZ_LOCAL].array(mfi);
      const auto& slopesEX = slopes[EX_LOCAL].array(mfi);
      const auto& slopesEY = slopes[EY_LOCAL].array(mfi);
      const auto& slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{		 
		  
		  Bc(i,j,k,ay) = (slopesBX(i+1,j,k,Y)+slopesBX(i,j,k,Y))/2.0;
		  Bc(i,j,k,axy) = (slopesBX(i+1,j,k,Y)-slopesBX(i,j,k,Y));
		  Bc(i,j,k,az) = 0.0;
		  Bc(i,j,k,axz) = 0.0;
		  Bc(i,j,k,bx) = (slopesBY(i,j+1,k,X)+slopesBY(i,j,k,X))/2.0;
		  Bc(i,j,k,bxy) = (slopesBY(i,j+1,k,X)-slopesBY(i,j,k,X));
		  Bc(i,j,k,bz) = 0.0;
		  Bc(i,j,k,byz) = 0.0;
		  Bc(i,j,k,cx) = slopesBZ(i,j,k,X);
		  Bc(i,j,k,cxz) = 0.0;
		  Bc(i,j,k,cy) = slopesBZ(i,j,k,Y);
		  Bc(i,j,k,cyz) = 0.0;
		  Bc(i,j,k,axx) = -dx[0]/2.0*(Bc(i,j,k,bxy)/dx[1]); // +cxz/dx[2]
		  Bc(i,j,k,byy) = -dx[1]/2.0*(Bc(i,j,k,axy)/dx[0]); // ++cyz/dx[2]
		  Bc(i,j,k,czz) = 0.0;
		  Bc(i,j,k,a0) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - Bc(i,j,k,axx)/6.0;
		  Bc(i,j,k,ax) = (arrEM_X(i+1,j,k,BX_LOCAL)-arrEM_X(i,j,k,BX_LOCAL));
		  Bc(i,j,k,b0) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - Bc(i,j,k,byy)/6.0;
		  Bc(i,j,k,by) = (arrEM_Y(i,j+1,k,BY_LOCAL)-arrEM_Y(i,j,k,BY_LOCAL));
		  Bc(i,j,k,c0) = arr(i,j,k,BZ) - Bc(i,j,k,czz)/6.0;
		  Bc(i,j,k,cz) = 0.0;

		  Ec(i,j,k,ay) = (slopesEX(i+1,j,k,Y)+slopesEX(i,j,k,Y))/2.0;
		  Ec(i,j,k,axy) = (slopesEX(i+1,j,k,Y)-slopesEX(i,j,k,Y));
		  Ec(i,j,k,az) = 0.0;
		  Ec(i,j,k,axz) = 0.0;
		  Ec(i,j,k,bx) = (slopesEY(i,j+1,k,X)+slopesEY(i,j,k,X))/2.0;
		  Ec(i,j,k,bxy) = (slopesEY(i,j+1,k,X)-slopesEY(i,j,k,X));
		  Ec(i,j,k,bz) = 0.0;
		  Ec(i,j,k,byz) = 0.0;
		  Ec(i,j,k,cx) = slopesEZ(i,j,k,X);
		  Ec(i,j,k,cxz) = 0.0;
		  Ec(i,j,k,cy) = slopesEZ(i,j,k,Y);
		  Ec(i,j,k,cyz) = 0.0;
		  Ec(i,j,k,axx) = -dx[0]/2.0*(Ec(i,j,k,bxy)/dx[1]-slopesQ(i,j,k,0)); // +cxz/dx[2]
		  Ec(i,j,k,byy) = -dx[1]/2.0*(Ec(i,j,k,axy)/dx[0]-slopesQ(i,j,k,1)); // ++cyz/dx[2]
		  Ec(i,j,k,czz) = 0.0;
		  Ec(i,j,k,a0) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - Ec(i,j,k,axx)/6.0;
		  Ec(i,j,k,ax) = (arrEM_X(i+1,j,k,EX_LOCAL)-arrEM_X(i,j,k,EX_LOCAL));
		  Ec(i,j,k,b0) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - Ec(i,j,k,byy)/6.0;
		  Ec(i,j,k,by) = (arrEM_Y(i,j+1,k,EY_LOCAL)-arrEM_Y(i,j,k,EY_LOCAL));
		  Ec(i,j,k,c0) = arr(i,j,k,EZ) - Ec(i,j,k,czz)/6.0; //czz*dx[2]*dx[2]/6.0
		  Ec(i,j,k,cz) = 0.0;		  

		}
	    }
	}
    }

  // Compute edge components of the EM fields    
  for (MFIter mfi(fluxesEM, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      // use old array, Sborder should also work
      //const auto& arr = Sborder.array(mfi);
      const auto& arr = S0.array(mfi);
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi); 
      const auto& fluxArrEM = fluxesEM.array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{		 
		  
		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
		  // LD state
		  x = dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = -1, jOffset = -1;
		  Vector<Real> EM_LD = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  		  
		  // RD state
		  x = -dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_RD = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // LU state
		  x = dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_LU = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // RU state
		  x = -dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_RU = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // If the reconstruction is consistent, then EyRU=EyRD
		  Real EyR = (EM_RU[EY_LOCAL]+EM_RD[EY_LOCAL])/2.0;
		  Real EyL = (EM_LU[EY_LOCAL]+EM_LD[EY_LOCAL])/2.0;
		  Real ExU = (EM_RU[EX_LOCAL]+EM_LU[EX_LOCAL])/2.0;
		  Real ExD = (EM_RD[EX_LOCAL]+EM_LD[EX_LOCAL])/2.0;
		  fluxArrEM(i,j,k,BZ_LOCAL) = 0.25*(EM_RU[BZ_LOCAL]+EM_RD[BZ_LOCAL]+EM_LU[BZ_LOCAL]+EM_LD[BZ_LOCAL])
		    - 0.5/c*(EyR-EyL) + 0.5/c*(ExU-ExD);

		  Real ByR = (EM_RU[BY_LOCAL]+EM_RD[BY_LOCAL])/2.0;
		  Real ByL = (EM_LU[BY_LOCAL]+EM_LD[BY_LOCAL])/2.0;
		  Real BxU = (EM_RU[BX_LOCAL]+EM_LU[BX_LOCAL])/2.0;
		  Real BxD = (EM_RD[BX_LOCAL]+EM_LD[BX_LOCAL])/2.0;
		  fluxArrEM(i,j,k,EZ_LOCAL) = 0.25*(EM_RU[EZ_LOCAL]+EM_RD[EZ_LOCAL]+EM_LU[EZ_LOCAL]+EM_LD[EZ_LOCAL])
		    + 0.5*c*(ByR-ByL) - 0.5*c*(BxU-BxD);
		}
	    }
	}       
    }

  // Update face-centred EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
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
	  //const auto& arr = Sborder.array(mfi);
	  const auto& arrEM = S_EM[d_EM].array(mfi);
	  //const auto& arrEMother = S_EM[d].array(mfi);
	  const auto& fluxArrEM = fluxesEM.array(mfi);
	  const auto& fluxArr = fluxes[d_EM].array(mfi);

	  const Dim3 hiDomain = ubound(geom.Domain());
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      // 2D code; z-component updated using cell-centred scheme
		      arrEM(i,j,k,BX_LOCAL+d_EM) += std::pow(-1,d)*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EZ_LOCAL) - fluxArrEM(i,j,k,EZ_LOCAL));
		      arrEM(i,j,k,EX_LOCAL+d_EM) -= std::pow(-1,d)*c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) - fluxArrEM(i,j,k,BZ_LOCAL));

		      // 3D code
		      /*arrEM(i,j,k,BX_LOCAL+(1+d)%3) += (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(2+d)%3));
			arrEM(i,j,k,BX_LOCAL+(2+d)%3) -= (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(1+d)%3));
			arrEM(i,j,k,EX_LOCAL+(1+d)%3) -= c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(2+d)%3));
			arrEM(i,j,k,EX_LOCAL+(2+d)%3) += c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(1+d)%3));
		      */
		      
		      // source terms
		      Real currentFace = r_i*fluxArr(i,j,k,RHO_I) + r_e*fluxArr(i,j,k,RHO_E);
		      arrEM(i,j,k,EX_LOCAL+d_EM) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;
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

  // Compute cell-centred EM fields from face-centred
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi); 

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

  		  arr(i,j,k,BX) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - Bc(i,j,k,axx)/6.0;
  		  arr(i,j,k,BY) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - Bc(i,j,k,byy)/6.0;
  		  arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - Ec(i,j,k,axx)/6.0;
  		  arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - Ec(i,j,k,byy)/6.0;
		  
  		}
  	    }
  	}       
    }

  // Only 2D
  // Compute edge components for y-components at x-faces
  for (MFIter mfi(fluxes[0], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // use old data
      //const auto& arr = Sborder.array(mfi);
      const auto& arr = S0.array(mfi);
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi); 
      //const auto& fluxArrEM = fluxesEM.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      //const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // L state
  		  x = dx[0]/2.0, y = 0.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_L = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		   
  		  // R state
  		  x = -dx[0]/2.0, y = 0.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_R = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // Note that this is not the flux, but the star states
  		  fluxArrX(i,j,k,BY) = 0.5*(EM_R[BY_LOCAL]+EM_L[BY_LOCAL])
		    + 0.5/c*(EM_R[EZ_LOCAL]-EM_L[EZ_LOCAL]);
  		  fluxArrX(i,j,k,EY) = 0.5*(EM_R[EY_LOCAL]+EM_L[EY_LOCAL])
		    - 0.5*c*(EM_R[BZ_LOCAL]-EM_L[BZ_LOCAL]);
  		}
  	    }
  	}
    }
  
  // Only 2D
  // Compute edge components for x-components at y-faces
  for (MFIter mfi(fluxes[1], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      // use old data
      //const auto& arr = Sborder.array(mfi);
      const auto& arr = S0.array(mfi);
      const auto& arrEM_X = S_EM[0].array(mfi);
      const auto& arrEM_Y = S_EM[1].array(mfi); 
      //const auto& fluxArrEM = fluxesEM.array(mfi);
      //const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 		  

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;
		  	
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;

  		  // D state
  		  x = 0.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_D = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
  		  // U state
  		  x = 0.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_U = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // Note that this is not the flux, but the star states		  
  		  fluxArrY(i,j,k,BX) = 0.5*(EM_U[BX_LOCAL]+EM_D[BX_LOCAL])
		    - 0.5/c*(EM_U[EZ_LOCAL]-EM_D[EZ_LOCAL]);
  		  fluxArrY(i,j,k,EX) = 0.5*(EM_U[EX_LOCAL]+EM_D[EX_LOCAL])
		    + 0.5*c*(EM_U[BZ_LOCAL]-EM_D[BZ_LOCAL]); 
  		}
  	    }
  	}
    }

  // Update cell-centred z-components of EM fields
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);
	        
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  // Update cell-centred z-components becuause it is 2D code
  		  arr(i,j,k,BZ) = arr(i,j,k,BZ) - (dt / dx[0]) * (fluxArrX(i+1,j,k,EY) - fluxArrX(i,j,k,EY)) + (dt / dx[1]) * (fluxArrY(i,j+1,k,EX) - fluxArrY(i,j,k,EX));
  		  arr(i,j,k,EZ) = arr(i,j,k,EZ) + c*c*(dt / dx[0]) * (fluxArrX(i+1,j,k,BY) - fluxArrX(i,j,k,BY)) - c*c*(dt / dx[1]) * (fluxArrY(i,j+1,k,BX) - fluxArrY(i,j,k,BX));		  

  		  // source terms
  		  arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
  		}
  	    }
  	}
    }    
}
void CAMReXmp::MaxwellSolverDivFreeNothing(Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  return;
}
void CAMReXmp::MaxwellSolverDivFreeWENO(Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{

  // Note that fluxesEM and fluxes are different definitions
  // fluxesEM are EM star states, fluxes are fluxes
  
  // Useful when dealing with high-order reconstruction
  // To get first-order slope in y-direction for Bx: slopes[BX_LOCALL].array(i,j,k,Y)
  // Useful to define indices int X,Y,Z in second-order linear reconstruction
  // For third-order quadratic reconstruction use int X,Y,Z,XX,YY,ZZ,XY,YZ,XZ
  // Although for x-components fields, only exist Y,Z,YY,ZZ,YZ (see Balsara et al. 2017)
  
  Array<MultiFab,6> slopes;
  // number of slopes
  // 9 slopes for third-order parabolic reconstruction
  // Ex = Ex^0 + Ex^y*(y/Delta_y) + Ex^z*(z/Delta_z) + Ex^yy*((y/Delta_y)^2-1/12) + ...
  // where Ex^y, Ex^z, Ex^yy are the slopes
  // For 2D code, only 2 slopes for x- and y-compoenents
  // For clarity, will use 9 slopes for third-order, first one represents slope in x-direction
  int nSlopes = 9;
  // indices for the slopes
  // for 2D, do not need Z indices
  int X=0,Y=1,Z=2,XX=3,YY=4,ZZ=5,XY=6,YZ=7,XZ=8;
  slopes[BX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);  
  slopes[EX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);
#if (AMREX_SPACEDIM >= 2)
  slopes[BY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
#else  
  // For 1D code it is cell-centred
  slopes[BY_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
  
#if (AMREX_SPACEDIM != 3)  
  // For 2D code it is cell-centred
  slopes[BZ_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EZ_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
  
  // Three components, one for each direction
  MultiFab slopesCharge(grids, dmap, nSlopes, 1);

  // Compute slopes of the charge densities  
  for (MFIter mfi(slopesCharge, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());
      
      // uses old charge
      const Array4<Real> arr = S_source.array(mfi);
      Array4<Real> slopes = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  Vector<Real> data_i = get_data_stencil(arr, i, j, k, 2, 0, 0, RHO_I);
		  Vector<Real> data_e = get_data_stencil(arr, i, j, k, 2, 0, 0, RHO_E);
		  Vector<Real> data_charge = get_charge_scaled(data_i,data_e);
		  std::array<Real, 2> slopesX = WENO3_slope(data_charge);
		  slopes(i,j,k,X) = slopesX[0];
		  slopes(i,j,k,XX) = slopesX[1];

		  data_i = get_data_stencil(arr, i, j, k, 0, 2, 0, RHO_I);
		  data_e = get_data_stencil(arr, i, j, k, 0, 2, 0, RHO_E);
		  data_charge = get_charge_scaled(data_i,data_e);
		  std::array<Real, 2> slopesY = WENO3_slope(data_charge);
		  slopes(i,j,k,Y) = slopesY[0];
		  slopes(i,j,k,YY) = slopesY[1];
		  
		  data_i = get_data_stencil(arr, i, j, k, 1, 1, 0, RHO_I);
		  data_e = get_data_stencil(arr, i, j, k, 1, 1, 0, RHO_E);
		  data_charge = get_charge_scaled(data_i,data_e);
		  Real slopesCross = WENO3_slopeCross(data_charge,
						      {slopes(i,j,k,X),slopes(i,j,k,XX),slopes(i,j,k,Y),slopes(i,j,k,YY)});
		  slopes(i,j,k,XY) = slopesCross;
		  
		  if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
		    {
		      int iOffset = 1, jOffset = 0, kOffset = 0;
		      Real u_iMinus1 = (r_i*arr(i-iOffset,j-jOffset,k-kOffset,RHO_I)
					+ r_e*arr(i-iOffset,j-jOffset,k-kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_iPlus1 = (r_i*arr(i+iOffset,j+jOffset,k+kOffset,RHO_I)
				       + r_e*arr(i+iOffset,j+jOffset,k+kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      slopes(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
		      slopes(i,j,k,XX) = 0.0, slopes(i,j,k,XY) = 0.0; 
		    }
		  if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
		    {
		      int iOffset = 0, jOffset = 1, kOffset = 0;
		      Real u_iMinus1 = (r_i*arr(i-iOffset,j-jOffset,k-kOffset,RHO_I)
					+ r_e*arr(i-iOffset,j-jOffset,k-kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_iPlus1 = (r_i*arr(i+iOffset,j+jOffset,k+kOffset,RHO_I)
				       + r_e*arr(i+iOffset,j+jOffset,k+kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      slopes(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
		      slopes(i,j,k,YY) = 0.0, slopes(i,j,k,XY) = 0.0; 			  
		    }
		}
	    }
	}      
    }
  
  // Compute cell-centred z-components slopes in x- and y-direction
  //for (int d = 0; d < AMREX_SPACEDIM ; d++)
    {

      //const int iOffset = ( d == 0 ? 1 : 0);
      //const int jOffset = ( d == 1 ? 1 : 0);
      //const int kOffset = ( d == 2 ? 1 : 0);
      
      for (MFIter mfi(slopes[BZ_LOCAL], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Dim3 hiDomain = ubound(geom.Domain());
	  
	  const Array4<Real> arr = S_source.array(mfi);
	  Array4<Real> slopesBZ = slopes[BZ_LOCAL].array(mfi);
	  Array4<Real> slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y-1; j <= hi.y+1; j++)
		{
		  for(int i = lo.x-1; i <= hi.x+1; i++)
		    {
		      // x slopes
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, BZ);
		      std::array<Real, 2> slopesX = WENO3_slope(dataX);
		      slopesBZ(i,j,k,X) = slopesX[0];
		      slopesBZ(i,j,k,XX) = slopesX[1];

		      dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, EZ);
		      slopesX = WENO3_slope(dataX);
		      slopesEZ(i,j,k,X) = slopesX[0];
		      slopesEZ(i,j,k,XX) = slopesX[1];
		      
		      // y slopes
		      dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, BZ);
		      slopesX = WENO3_slope(dataX);
		      slopesBZ(i,j,k,Y) = slopesX[0];
		      slopesBZ(i,j,k,YY) = slopesX[1];
		     		      
		      dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, EZ);
		      slopesX = WENO3_slope(dataX);
		      slopesEZ(i,j,k,Y) = slopesX[0];
		      slopesEZ(i,j,k,YY) = slopesX[1];
		      
		      // cross slopes
		      Vector<Real> dataXY = get_data_stencil(arr, i, j, k, 1, 1, 0, BZ);
		      Real slopesCross = WENO3_slopeCross(dataXY,
							  {slopesBZ(i,j,k,X),slopesBZ(i,j,k,XX),slopesBZ(i,j,k,Y),slopesBZ(i,j,k,YY)});
		      slopesBZ(i,j,k,XY) = slopesCross;

		      dataXY = get_data_stencil(arr, i, j, k, 1, 1, 0, EZ);
		      slopesCross = WENO3_slopeCross(dataXY,
						     {slopesEZ(i,j,k,X),slopesEZ(i,j,k,XX),slopesEZ(i,j,k,Y),slopesEZ(i,j,k,YY)});
		      slopesEZ(i,j,k,XY) = slopesCross;		      		      

		      if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
			{			      
			  {
			    Real u_iMinus1 = arr(i-1,j,k,BZ);
			    Real u_i = arr(i,j,k,BZ);
			    Real u_iPlus1 = arr(i+1,j,k,BZ);
			    slopesBZ(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesBZ(i,j,k,XX) = 0.0, slopesBZ(i,j,k,XY) = 0.0;
			  }
			  {				  
			    Real u_iMinus1 = arr(i-1,j,k,EZ);
			    Real u_i = arr(i,j,k,EZ);
			    Real u_iPlus1 = arr(i+1,j,k,EZ);
			    slopesEZ(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesEZ(i,j,k,XX) = 0.0, slopesEZ(i,j,k,XY) = 0.0;
			  }			      
			}
		      if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
			{
			  {			
			    Real u_iMinus1 = arr(i,j-1,k,BZ);
			    Real u_i = arr(i,j,k,BZ);
			    Real u_iPlus1 = arr(i,j+1,k,BZ);
			    slopesBZ(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesBZ(i,j,k,YY) = 0.0, slopesBZ(i,j,k,XY) = 0.0;
			    
			  }			  
			  {
			    Real u_iMinus1 = arr(i,j-1,k,EZ);
			    Real u_i = arr(i,j,k,EZ);
			    Real u_iPlus1 = arr(i,j+1,k,EZ);
			    slopesEZ(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesEZ(i,j,k,YY) = 0.0, slopesEZ(i,j,k,XY) = 0.0;
			  }			      
			}
		    }
		}
	    }      
	}
    }
  
  // Compute y-components slopes in x-direction
  for (MFIter mfi(slopes[BY_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());

      const Array4<Real> arr = S_EM_source[1].array(mfi);
      Array4<Real> slopesBY = slopes[BY_LOCAL].array(mfi);
      Array4<Real> slopesEY = slopes[EY_LOCAL].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  Vector<Real> dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, BY_LOCAL);
		  std::array<Real, 2> slopesX = WENO3_slope(dataX);
		  slopesBY(i,j,k,X) = slopesX[0];
		  slopesBY(i,j,k,XX) = slopesX[1];

		  dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, EY_LOCAL);
		  slopesX = WENO3_slope(dataX);
		  slopesEY(i,j,k,X) = slopesX[0];
		  slopesEY(i,j,k,XX) = slopesX[1];
		  
		  if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
		    {			      
		      {
			Real u_iMinus1 = arr(i-1,j,k,BY_LOCAL);
			Real u_i = arr(i,j,k,BY_LOCAL);		    
			Real u_iPlus1 = arr(i+1,j,k,BY_LOCAL);
			slopesBY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesBY(i,j,k,XX) = 0.0;
		      }
		      {
			Real u_iMinus1 = arr(i-1,j,k,EY_LOCAL);
			Real u_i = arr(i,j,k,EY_LOCAL);		    
			Real u_iPlus1 = arr(i+1,j,k,EY_LOCAL);
			slopesEY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesEY(i,j,k,XX) = 0.0;
		      }			      
		    }
		}
	    }
	}      
    }
  // Compute x-components slopes in y-direction
  for (MFIter mfi(slopes[BX_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());

      const Array4<Real> arr = S_EM_source[0].array(mfi);
      Array4<Real> slopesBX = slopes[BX_LOCAL].array(mfi);
      Array4<Real> slopesEX = slopes[EX_LOCAL].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  Vector<Real> dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, BX_LOCAL);
		  std::array<Real, 2> slopesX = WENO3_slope(dataX);
		  slopesBX(i,j,k,Y) = slopesX[0];
		  slopesBX(i,j,k,YY) = slopesX[1];

		  dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, EX_LOCAL);
		  slopesX = WENO3_slope(dataX);
		  slopesEX(i,j,k,Y) = slopesX[0];
		  slopesEX(i,j,k,YY) = slopesX[1];
		  
		  if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
		    {
		      {
			Real u_iMinus1 = arr(i,j-1,k,BX_LOCAL);
			Real u_i = arr(i,j,k,BX_LOCAL);		    
			Real u_iPlus1 = arr(i,j+1,k,BX_LOCAL);
			slopesBX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesBX(i,j,k,YY) = 0.0;			      
		      }
		      {
			Real u_iMinus1 = arr(i,j-1,k,EX_LOCAL);
			Real u_i = arr(i,j,k,EX_LOCAL);		    
			Real u_iPlus1 = arr(i,j+1,k,EX_LOCAL);
			slopesEX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesEX(i,j,k,YY) = 0.0;			      
		      }			      
		    }
		}
	    }
	}      
    }


  // Coefficients for the magnetic and electric fields
  // For third order, there are 48 coeff.
  // i.e. a0,ax,ay,az,axx,...,b0,bx,...,c0,cx,...,czz
  // use 1 ghost cell
  int a0=0,ax=1,ay=2,az=3,axx=4,ayy=5,azz=6,axy=7,ayz=8,axz=9,axxx=10,axxy=11,axxz=12,axyy=13,axzz=14,axyz=15;
  int b0=16,bx=17,by=18,bz=19,bxx=20,byy=21,bzz=22,bxy=23,byz=24,bxz=25,byyy=26,bxyy=27,byyz=28,bxxy=29,byzz=30,bxyz=31;
  int c0=32,cx=33,cy=34,cz=35,cxx=36,cyy=37,czz=38,cxy=39,cyz=40,cxz=41,czzz=42,cxzz=43,cyzz=44,cxxz=45,cyyz=46,cxyz=47;
  MultiFab Bcoeff(grids, dmap, 48, 1);
  MultiFab Ecoeff(grids, dmap, 48, 1);

  // Compute coefficients for the field function
  for (MFIter mfi(Bcoeff, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      // coeff
      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      // data
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi);       

      // slopes
      const auto& slopesBX = slopes[BX_LOCAL].array(mfi);
      const auto& slopesBY = slopes[BY_LOCAL].array(mfi);
      const auto& slopesBZ = slopes[BZ_LOCAL].array(mfi);
      const auto& slopesEX = slopes[EX_LOCAL].array(mfi);
      const auto& slopesEY = slopes[EY_LOCAL].array(mfi);
      const auto& slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  // third order terms
		  Bc(i,j,k,ayy) = (slopesBX(i+1,j,k,YY)+slopesBX(i,j,k,YY))/2.0;
		  Bc(i,j,k,axyy) = slopesBX(i+1,j,k,YY)-slopesBX(i,j,k,YY);
		  Bc(i,j,k,azz) = 0.0;
		  Bc(i,j,k,axzz) = 0.0;
		  Bc(i,j,k,ayz) = 0.0;
		  Bc(i,j,k,axyz) = 0.0;
		  Bc(i,j,k,bxx) = (slopesBY(i,j+1,k,XX)+slopesBY(i,j,k,XX))/2.0;
		  Bc(i,j,k,bxxy) = slopesBY(i,j+1,k,XX)-slopesBY(i,j,k,XX);
		  Bc(i,j,k,bzz) = 0.0;
		  Bc(i,j,k,byzz) = 0.0;
		  Bc(i,j,k,bxz) = 0.0;
		  Bc(i,j,k,bxyz) = 0.0;
		  Bc(i,j,k,cxx) = slopesBZ(i,j,k,XX);
		  Bc(i,j,k,cxxz) = 0.0;
		  Bc(i,j,k,cyy) = slopesBZ(i,j,k,YY);
		  Bc(i,j,k,cyyz) = 0.0;
		  Bc(i,j,k,cxy) = slopesBZ(i,j,k,XY);
		  Bc(i,j,k,cxyz) = 0.0;
		  Bc(i,j,k,axxx) = -dx[0]/3.0*(Bc(i,j,k,bxxy)/dx[1]);
		  Bc(i,j,k,byyy) = -dx[1]/3.0*(Bc(i,j,k,axyy)/dx[0]);
		  Bc(i,j,k,czzz) = 0.0;
		  Bc(i,j,k,axxy) = 0.0;
		  Bc(i,j,k,bxyy) = 0.0;
		  Bc(i,j,k,byyz) = 0.0;
		  Bc(i,j,k,cyzz) = 0.0;
		  Bc(i,j,k,cxzz) = 0.0;
		  Bc(i,j,k,axxz) = 0.0;
		  // second order terms
		  Bc(i,j,k,ay) = (slopesBX(i+1,j,k,Y)+slopesBX(i,j,k,Y))/2.0 - Bc(i,j,k,axxy)/6.0;
		  Bc(i,j,k,axy) = (slopesBX(i+1,j,k,Y)-slopesBX(i,j,k,Y));
		  Bc(i,j,k,az) = 0.0;
		  Bc(i,j,k,axz) = 0.0;
		  Bc(i,j,k,bx) = (slopesBY(i,j+1,k,X)+slopesBY(i,j,k,X))/2.0 - Bc(i,j,k,bxyy)/6.0;
		  Bc(i,j,k,bxy) = (slopesBY(i,j+1,k,X)-slopesBY(i,j,k,X));
		  Bc(i,j,k,bz) = 0.0;
		  Bc(i,j,k,byz) = 0.0;
		  Bc(i,j,k,cx) = slopesBZ(i,j,k,X) - Bc(i,j,k,cxzz)/6.0;
		  Bc(i,j,k,cxz) = 0.0;
		  Bc(i,j,k,cy) = slopesBZ(i,j,k,Y) - Bc(i,j,k,cyzz)/6.0;
		  Bc(i,j,k,cyz) = 0.0;
		  Bc(i,j,k,axx) = -dx[0]/2.0*(Bc(i,j,k,bxy)/dx[1]); // +cxz/dx[2]
		  Bc(i,j,k,byy) = -dx[1]/2.0*(Bc(i,j,k,axy)/dx[0]); // ++cyz/dx[2]
		  Bc(i,j,k,czz) = 0.0;
		  Bc(i,j,k,a0) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - Bc(i,j,k,axx)/6.0;
		  Bc(i,j,k,ax) = (arrEM_X(i+1,j,k,BX_LOCAL)-arrEM_X(i,j,k,BX_LOCAL)) - Bc(i,j,k,axxx)/10.0;
		  Bc(i,j,k,b0) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - Bc(i,j,k,byy)/6.0;
		  Bc(i,j,k,by) = (arrEM_Y(i,j+1,k,BY_LOCAL)-arrEM_Y(i,j,k,BY_LOCAL)) - Bc(i,j,k,byyy)/10.0;
		  Bc(i,j,k,c0) = arr(i,j,k,BZ) - Bc(i,j,k,czz)/6.0;
		  Bc(i,j,k,cz) = 0.0;

		  // third order terms
		  Ec(i,j,k,ayy) = (slopesEX(i+1,j,k,YY)+slopesEX(i,j,k,YY))/2.0;
		  Ec(i,j,k,axyy) = slopesEX(i+1,j,k,YY)-slopesEX(i,j,k,YY);
		  Ec(i,j,k,azz) = 0.0;
		  Ec(i,j,k,axzz) = 0.0;
		  Ec(i,j,k,ayz) = 0.0;
		  Ec(i,j,k,axyz) = 0.0;
		  Ec(i,j,k,bxx) = (slopesEY(i,j+1,k,XX)+slopesEY(i,j,k,XX))/2.0;
		  Ec(i,j,k,bxxy) = slopesEY(i,j+1,k,XX)-slopesEY(i,j,k,XX);
		  Ec(i,j,k,bzz) = 0.0;
		  Ec(i,j,k,byzz) = 0.0;
		  Ec(i,j,k,bxz) = 0.0;
		  Ec(i,j,k,bxyz) = 0.0;
		  Ec(i,j,k,cxx) = slopesEZ(i,j,k,XX);
		  Ec(i,j,k,cxxz) = 0.0;
		  Ec(i,j,k,cyy) = slopesEZ(i,j,k,YY);
		  Ec(i,j,k,cyyz) = 0.0;
		  Ec(i,j,k,cxy) = slopesEZ(i,j,k,XY);
		  Ec(i,j,k,cxyz) = 0.0;
		  Ec(i,j,k,axxx) = -dx[0]/3.0*(Ec(i,j,k,bxxy)/dx[1]-slopesQ(i,j,k,XX));
		  Ec(i,j,k,byyy) = -dx[1]/3.0*(Ec(i,j,k,axyy)/dx[0]-slopesQ(i,j,k,YY));
		  Ec(i,j,k,czzz) = 0.0;
		  Ec(i,j,k,axxy) = dx[0]*dx[1]/(2.0*(dx[0]+dx[1]))*slopesQ(i,j,k,XY);
		  Ec(i,j,k,bxyy) = Ec(i,j,k,axxy);
		  Ec(i,j,k,byyz) = 0.0;
		  Ec(i,j,k,cyzz) = 0.0;
		  Ec(i,j,k,cxzz) = 0.0;
		  Ec(i,j,k,axxz) = 0.0;
		  // second order terms
		  Ec(i,j,k,ay) = (slopesEX(i+1,j,k,Y)+slopesEX(i,j,k,Y))/2.0 - Ec(i,j,k,axxy)/6.0;
		  Ec(i,j,k,axy) = (slopesEX(i+1,j,k,Y)-slopesEX(i,j,k,Y));
		  Ec(i,j,k,az) = 0.0;
		  Ec(i,j,k,axz) = 0.0;
		  Ec(i,j,k,bx) = (slopesEY(i,j+1,k,X)+slopesEY(i,j,k,X))/2.0 - Ec(i,j,k,bxyy)/6.0;
		  Ec(i,j,k,bxy) = (slopesEY(i,j+1,k,X)-slopesEY(i,j,k,X));
		  Ec(i,j,k,bz) = 0.0;
		  Ec(i,j,k,byz) = 0.0;
		  Ec(i,j,k,cx) = slopesEZ(i,j,k,X) - Ec(i,j,k,cxzz)/6.0;
		  Ec(i,j,k,cxz) = 0.0;
		  Ec(i,j,k,cy) = slopesEZ(i,j,k,Y) - Ec(i,j,k,cyzz)/6.0;
		  Ec(i,j,k,cyz) = 0.0;
		  Ec(i,j,k,axx) = -dx[0]/2.0*(Ec(i,j,k,bxy)/dx[1]-slopesQ(i,j,k,X)); // +cxz/dx[2]
		  Ec(i,j,k,byy) = -dx[1]/2.0*(Ec(i,j,k,axy)/dx[0]-slopesQ(i,j,k,Y)); // ++cyz/dx[2]
		  Ec(i,j,k,czz) = 0.0;
		  Ec(i,j,k,a0) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - Ec(i,j,k,axx)/6.0;
		  Ec(i,j,k,ax) = (arrEM_X(i+1,j,k,EX_LOCAL)-arrEM_X(i,j,k,EX_LOCAL)) - Ec(i,j,k,axxx)/10.0;
		  Ec(i,j,k,b0) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - Ec(i,j,k,byy)/6.0;
		  Ec(i,j,k,by) = (arrEM_Y(i,j+1,k,EY_LOCAL)-arrEM_Y(i,j,k,EY_LOCAL)) - Ec(i,j,k,byyy)/10.0;
		  Ec(i,j,k,c0) = arr(i,j,k,EZ) - Ec(i,j,k,czz)/6.0; //czz*dx[2]*dx[2]/6.0
		  Ec(i,j,k,cz) = 0.0;		  

		}
	    }
	}
    }

  // Compute edge components of the EM fields    
  for (MFIter mfi(fluxesEM, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      const Dim3 hiDomain = ubound(geom.Domain());
      
      // use old array, S_source should also work
      //const auto& arr = S_source.array(mfi);
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi); 
      const auto& fluxArrEM = fluxesEM.array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{		 
		  
		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
		  // LD state
		  x = dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = -1, jOffset = -1;
		  Vector<Real> EM_LD = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  		  
		  // RD state
		  x = -dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_RD = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // LU state
		  x = dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_LU = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // RU state
		  x = -dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_RU = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // If the reconstruction is consistent, then EyRU=EyRD
		  Real EyR = (EM_RU[EY_LOCAL]+EM_RD[EY_LOCAL])/2.0;
		  Real EyL = (EM_LU[EY_LOCAL]+EM_LD[EY_LOCAL])/2.0;
		  Real ExU = (EM_RU[EX_LOCAL]+EM_LU[EX_LOCAL])/2.0;
		  Real ExD = (EM_RD[EX_LOCAL]+EM_LD[EX_LOCAL])/2.0;
		  fluxArrEM(i,j,k,BZ_LOCAL) = 0.25*(EM_RU[BZ_LOCAL]+EM_RD[BZ_LOCAL]+EM_LU[BZ_LOCAL]+EM_LD[BZ_LOCAL])
		    - 0.5/c*(EyR-EyL) + 0.5/c*(ExU-ExD);

		  Real ByR = (EM_RU[BY_LOCAL]+EM_RD[BY_LOCAL])/2.0;
		  Real ByL = (EM_LU[BY_LOCAL]+EM_LD[BY_LOCAL])/2.0;
		  Real BxU = (EM_RU[BX_LOCAL]+EM_LU[BX_LOCAL])/2.0;
		  Real BxD = (EM_RD[BX_LOCAL]+EM_LD[BX_LOCAL])/2.0;

		  fluxArrEM(i,j,k,EZ_LOCAL) = 0.25*(EM_RU[EZ_LOCAL]+EM_RD[EZ_LOCAL]+EM_LU[EZ_LOCAL]+EM_LD[EZ_LOCAL])
		    + 0.5*c*(ByR-ByL) - 0.5*c*(BxU-BxD);
		  		  
		}
	    }
	}       
    }

  // Update face-centred EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      int d_EM = (d==0) ? 1 : 0;
    
      // Loop over all the patches at this level
      for (MFIter mfi(S_EM_dest[d_EM], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();

	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  // Indexable arrays for the data, and the directional flux
	  // Based on the corner-centred definition of the flux array, the
	  // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	  //const auto& arr = S_source.array(mfi);
	  const auto& arrEM = S_EM_dest[d_EM].array(mfi);
	  const auto& arrEMOld = S_EM_source[d_EM].array(mfi);	  
	  const auto& arr = S_source.array(mfi);
	  const auto& fluxArrEM = fluxesEM.array(mfi);
	  const auto& fluxArr = fluxes[d_EM].array(mfi);

	  const Dim3 hiDomain = ubound(geom.Domain());
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      // 2D code; z-component updated using cell-centred scheme
		      arrEM(i,j,k,BX_LOCAL+d_EM) = arrEMOld(i,j,k,BX_LOCAL+d_EM) + std::pow(-1,d)*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EZ_LOCAL) - fluxArrEM(i,j,k,EZ_LOCAL));
		      arrEM(i,j,k,EX_LOCAL+d_EM) = arrEMOld(i,j,k,EX_LOCAL+d_EM) - std::pow(-1,d)*c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) - fluxArrEM(i,j,k,BZ_LOCAL));
		      
		      // source terms
		      Real currentFace = r_i*fluxArr(i,j,k,RHO_I) + r_e*fluxArr(i,j,k,RHO_E);
		      arrEM(i,j,k,EX_LOCAL+d_EM) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;
		      
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

  // Re-compute some slopes to obtain 3rd order cell centred averages
  // However, in 1D problem, Bx should be constant, but the bxy component (from axx) is not constant
  // Ignore also to save some computational effort
  // // Compute y-components slopes in x-direction
  // for (MFIter mfi(slopes[BY_LOCAL], true); mfi.isValid(); ++mfi)
  //   {
  //     const Box& bx = mfi.tilebox();
      
  //     const Dim3 lo = lbound(bx);
  //     const Dim3 hi = ubound(bx);

  //     const Dim3 hiDomain = ubound(geom.Domain());

  //     const Array4<Real> arr = S_EM_dest[1].array(mfi);
  //     Array4<Real> slopesBY = slopes[BY_LOCAL].array(mfi);
  //     Array4<Real> slopesEY = slopes[EY_LOCAL].array(mfi);
      
  //     for(int k = lo.z; k <= hi.z; k++)
  // 	{
  // 	  for(int j = lo.y-1; j <= hi.y+1; j++)
  // 	    {
  // 	      for(int i = lo.x-1; i <= hi.x+1; i++)
  // 		{
  // 		  Vector<Real> dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, BY_LOCAL);
  // 		  std::array<Real, 2> slopesX = WENO3_slope(dataX);
  // 		  slopesBY(i,j,k,X) = slopesX[0];
		  
  // 		  dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, EY_LOCAL);
  // 		  slopesX = WENO3_slope(dataX);
  // 		  slopesEY(i,j,k,X) = slopesX[0];
		  
  // 		    }
  // 		}
  // 	    }
  // 	}      
  //   }
  // // Compute x-components slopes in y-direction
  // for (MFIter mfi(slopes[BX_LOCAL], true); mfi.isValid(); ++mfi)
  //   {
  //     const Box& bx = mfi.tilebox();
      
  //     const Dim3 lo = lbound(bx);
  //     const Dim3 hi = ubound(bx);

  //     const Dim3 hiDomain = ubound(geom.Domain());

  //     const Array4<Real> arr = S_EM_dest[0].array(mfi);
  //     Array4<Real> slopesBX = slopes[BX_LOCAL].array(mfi);
  //     Array4<Real> slopesEX = slopes[EX_LOCAL].array(mfi);
      
  //     for(int k = lo.z; k <= hi.z; k++)
  // 	{
  // 	  for(int j = lo.y-1; j <= hi.y+1; j++)
  // 	    {
  // 	      for(int i = lo.x-1; i <= hi.x+1; i++)
  // 		{
  // 		  Vector<Real> dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, BX_LOCAL);
  // 		  std::array<Real, 2> slopesX = WENO3_slope(dataX);
  // 		  slopesBX(i,j,k,Y) = slopesX[0];

  // 		  dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, EX_LOCAL);
  // 		  slopesX = WENO3_slope(dataX);
  // 		  slopesEX(i,j,k,Y) = slopesX[0];

  // 		    }
  // 		}
  // 	    }
  // 	}      
  //   }

  // // Compute coefficients for the field function
  // for (MFIter mfi(Bcoeff, true); mfi.isValid(); ++mfi)
  //   {
  //     const Box& box = mfi.tilebox();
      
  //     const Dim3 lo = lbound(box);
  //     const Dim3 hi = ubound(box);

  //     // coeff
  //     const auto& Bc = Bcoeff.array(mfi);
  //     const auto& Ec = Ecoeff.array(mfi);
      
  //     // slopes
  //     const auto& slopesBX = slopes[BX_LOCAL].array(mfi);
  //     const auto& slopesBY = slopes[BY_LOCAL].array(mfi);
  //     const auto& slopesEX = slopes[EX_LOCAL].array(mfi);
  //     const auto& slopesEY = slopes[EY_LOCAL].array(mfi);
      
  //     const auto& slopesQ = slopesCharge.array(mfi);
      
  //     for(int k = lo.z; k <= hi.z; k++)
  // 	{
  // 	  for(int j = lo.y-1; j <= hi.y+1; j++)
  // 	    {
  // 	      for(int i = lo.x-1; i <= hi.x+1; i++)
  // 		{

  // 		  Bc(i,j,k,axy) = (slopesBX(i+1,j,k,Y)-slopesBX(i,j,k,Y));
  // 		  Bc(i,j,k,bxy) = (slopesBY(i,j+1,k,X)-slopesBY(i,j,k,X));
  // 		  Bc(i,j,k,axx) = -dx[0]/2.0*(Bc(i,j,k,bxy)/dx[1]); // +cxz/dx[2]
  // 		  Bc(i,j,k,byy) = -dx[1]/2.0*(Bc(i,j,k,axy)/dx[0]); // ++cyz/dx[2]

  // 		  Ec(i,j,k,axy) = (slopesEX(i+1,j,k,Y)-slopesEX(i,j,k,Y));
  // 		  Ec(i,j,k,bxy) = (slopesEY(i,j+1,k,X)-slopesEY(i,j,k,X));
  // 		  Ec(i,j,k,axx) = -dx[0]/2.0*(Ec(i,j,k,bxy)/dx[1]-slopesQ(i,j,k,X)); // +cxz/dx[2]
  // 		  Ec(i,j,k,byy) = -dx[1]/2.0*(Ec(i,j,k,axy)/dx[0]-slopesQ(i,j,k,Y)); // ++cyz/dx[2]

  // 		}
  // 	    }
  // 	}
  //   }

  // Compute cell-centred EM fields from face-centred
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrEM_X = S_EM_dest[0].array(mfi);
      const auto& arrEM_Y = S_EM_dest[1].array(mfi); 

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 
  		  arr(i,j,k,BX) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0;// - Bc(i,j,k,axx)/6.0;
  		  arr(i,j,k,BY) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0;// - Bc(i,j,k,byy)/6.0;
  		  arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0;// - Ec(i,j,k,axx)/6.0;
  		  arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0;// - Ec(i,j,k,byy)/6.0;
		  
  		}
  	    }
  	}       
    }

  // Only 2D
  // Compute edge components for y-components at x-faces
  for (MFIter mfi(fluxes[0], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      const auto& fluxArrX = fluxes[0].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;
		  /*
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // L state
  		  x = dx[0]/2.0, y = 0.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_L = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		   
  		  // R state
  		  x = -dx[0]/2.0, y = 0.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_R = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  */
		  // 2 gaussian points
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // L state
		  iOffset = -1, jOffset = 0;
  		  x = dx[0]/2.0, y = dx[1]/(2.0*std::sqrt(3.0));		  
		  Vector<Real> EM_L1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = dx[0]/2.0, y = -dx[1]/(2.0*std::sqrt(3.0));		  
		  Vector<Real> EM_L2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);		  
		   
  		  // R state
		  iOffset = 0, jOffset = 0;
  		  x = -dx[0]/2.0, y = dx[1]/(2.0*std::sqrt(3.0));		 
		  Vector<Real> EM_R1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = -dx[0]/2.0, y = -dx[1]/(2.0*std::sqrt(3.0));		 
		  Vector<Real> EM_R2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  Vector<Real> EM_L, EM_R;
		  for (int n=0; n<6; n++)
		    {
		      EM_L.push_back(0.5*(EM_L1[n]+EM_L2[n]));
		      EM_R.push_back(0.5*(EM_R1[n]+EM_R2[n]));
		    }
		  
		  // Note that this is not the flux, but the star states
  		  fluxArrX(i,j,k,BY) = 0.5*(EM_R[BY_LOCAL]+EM_L[BY_LOCAL])
		    + 0.5/c*(EM_R[EZ_LOCAL]-EM_L[EZ_LOCAL]);
  		  fluxArrX(i,j,k,EY) = 0.5*(EM_R[EY_LOCAL]+EM_L[EY_LOCAL])
		    - 0.5*c*(EM_R[BZ_LOCAL]-EM_L[BZ_LOCAL]);
  		}
  	    }
  	}
    }
  
  // Only 2D
  // Compute edge components for x-components at y-faces
  for (MFIter mfi(fluxes[1], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 		  

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;
		  /*		  	
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;

  		  // D state
  		  x = 0.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_D = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
  		  // U state
  		  x = 0.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_U = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  */
		  // 2 gaussian points
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // D state
  		  iOffset = 0, jOffset = -1;
		  x = dx[0]/(2.0*std::sqrt(3.0)), y = dx[1]/2.0;		  
		  Vector<Real> EM_D1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = -dx[0]/(2.0*std::sqrt(3.0)), y = dx[1]/2.0;		  
		  Vector<Real> EM_D2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
  		  // U state
  		  iOffset = 0, jOffset = 0;
		  x = dx[0]/(2.0*std::sqrt(3.0)), y = -dx[1]/2.0;		  
		  Vector<Real> EM_U1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = -dx[0]/(2.0*std::sqrt(3.0)), y = -dx[1]/2.0;		  
		  Vector<Real> EM_U2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  Vector<Real> EM_D, EM_U;
		  for (int n=0; n<6; n++)
		    {
		      EM_D.push_back(0.5*(EM_D1[n]+EM_D2[n]));
		      EM_U.push_back(0.5*(EM_U1[n]+EM_U2[n]));
		    }
		  
		  // Note that this is not the flux, but the star states 	  
  		  fluxArrY(i,j,k,BX) = 0.5*(EM_U[BX_LOCAL]+EM_D[BX_LOCAL])
		    - 0.5/c*(EM_U[EZ_LOCAL]-EM_D[EZ_LOCAL]);
  		  fluxArrY(i,j,k,EX) = 0.5*(EM_U[EX_LOCAL]+EM_D[EX_LOCAL])
		    + 0.5*c*(EM_U[BZ_LOCAL]-EM_D[BZ_LOCAL]); 
  		}
  	    }
  	}
    }

  // Update cell-centred z-components of EM fields
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrOld = S_source.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);

      Array4<Real> slopes = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  // Update cell-centred z-components becuause it is 2D code
  		  arr(i,j,k,BZ) = arrOld(i,j,k,BZ) - (dt / dx[0]) * (fluxArrX(i+1,j,k,EY) - fluxArrX(i,j,k,EY)) + (dt / dx[1]) * (fluxArrY(i,j+1,k,EX) - fluxArrY(i,j,k,EX));
  		  arr(i,j,k,EZ) = arrOld(i,j,k,EZ) + c*c*(dt / dx[0]) * (fluxArrX(i+1,j,k,BY) - fluxArrX(i,j,k,BY)) - c*c*(dt / dx[1]) * (fluxArrY(i,j+1,k,BX) - fluxArrY(i,j,k,BX));		  

  		  // source terms
  		  arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
		  // 4 gaussian points
		  /*Real data_charge = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
		  Real slopesX = slopes(i,j,k,X);
		  Real slopesXX = slopes(i,j,k,XX);
		  Real slopesY = slopes(i,j,k,Y);
		  Real slopesYY = slopes(i,j,k,YY);
		  Real slopesXY = slopes(i,j,k,XY);
		  // 1st quadrature
		  Real x = -0.5, y = 0.0;
		  Real Lx = x, Lxx = x*x - 1.0/12.0;
		  Real Ly = y, Lyy = y*y - 1.0/12.0;
		  Real data_charge1 = data_charge + slopesX*Lx + slopesXX*Lxx +
		    slopesY*Ly + slopesYY*Lyy + slopesXY*Lx*Ly;
		  // 2nd quadrature
		  x = 0.5, y = 0.0;
		  Lx = x, Lxx = x*x - 1.0/12.0;
		  Ly = y, Lyy = y*y - 1.0/12.0;
		  Real data_charge2 = data_charge + slopesX*Lx + slopesXX*Lxx +
		    slopesY*Ly + slopesYY*Lyy + slopesXY*Lx*Ly;
		  // 3nd quadrature
		  x = 0.0, y = 0.5;
		  Lx = x, Lxx = x*x - 1.0/12.0;
		  Ly = y, Lyy = y*y - 1.0/12.0;
		  Real data_charge3 = data_charge + slopesX*Lx + slopesXX*Lxx +
		    slopesY*Ly + slopesYY*Lyy + slopesXY*Lx*Ly;
		  // 2nd quadrature
		  x = 0.0, y = -0.5;
		  Lx = x, Lxx = x*x - 1.0/12.0;
		  Ly = y, Lyy = y*y - 1.0/12.0;
		  Real data_charge4 = data_charge + slopesX*Lx + slopesXX*Lxx +
		    slopesY*Ly + slopesYY*Lyy + slopesXY*Lx*Ly;
		  arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*0.25*(data_charge1+data_charge2+
							   data_charge3+data_charge4);		  
		  */
  		}
  	    }
  	}
    }
  // We need to compute boundary conditions again after each update
  S_dest.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_dest, geom, bc);  

}
void CAMReXmp::MaxwellSolverDivFreeTVD(Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{

  // Note that fluxesEM and fluxes are different definitions
  // fluxesEM are EM star states, fluxes are fluxes
  
  // Useful when dealing with high-order reconstruction
  // To get first-order slope in y-direction for Bx: slopes[BX_LOCALL].array(i,j,k,Y)
  // Useful to define indices int X,Y,Z in second-order linear reconstruction
  // For third-order quadratic reconstruction use int X,Y,Z,XX,YY,ZZ,XY,YZ,XZ
  // Although for x-components fields, only exist Y,Z,YY,ZZ,YZ (see Balsara et al. 2017)
  
  Array<MultiFab,6> slopes;
  // number of slopes
  // 2 slopes for second-order linear reconstruction
  // Ex = Ex^0 + Ex^y*(y/Delta_y) + Ex^z*(z/Delta_z), where Ex^y and Ex^z are the slopes
  // For 2D code, only 1 slope for x- and y-compoenents
  // For clarity, will use 3 slopes for second-order, first one represents slope in x-direction
  // 9 slopes for third-order, and so on
  int nSlopes = 3;
  // indices for the slopes
  // for 2D, do not need Z index
  // for third order, use also XX,YY,ZZ,XY,YZ,XZ
  int X=0,Y=1,Z=2;
  slopes[BX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);  
  slopes[EX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);
#if (AMREX_SPACEDIM >= 2)
  slopes[BY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
#else  
  // For 1D code it is cell-centred
  slopes[BY_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
  
#if (AMREX_SPACEDIM != 3)  
  // For 2D code it is cell-centred
  slopes[BZ_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EZ_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
    
  // Three components, one for each direction
  MultiFab slopesCharge(grids, dmap, AMREX_SPACEDIM, 1);

  // Compute slopes of the charge densities
  for (int d = 0; d < AMREX_SPACEDIM ; d++)   
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      // Compute cell-centred density slopes
      for (MFIter mfi(slopesCharge, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
 	  const Dim3 hi = ubound(bx);

	  // uses old charge
	  const Array4<Real> arr = S_source.array(mfi);
	  Array4<Real> slopes = slopesCharge.array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y-1; j <= hi.y+1; j++)
		{
		  for(int i = lo.x-1; i <= hi.x+1; i++)
		    {
		      Real u_iMinus1 = (r_i*arr(i-iOffset,j-jOffset,k-kOffset,RHO_I)
					+ r_e*arr(i-iOffset,j-jOffset,k-kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_iPlus1 = (r_i*arr(i+iOffset,j+jOffset,k+kOffset,RHO_I)
				       + r_e*arr(i+iOffset,j+jOffset,k+kOffset,RHO_E))/(lambda_d*lambda_d*l_r);

		      // slopes for the current
		      slopes(i,j,k,d) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		    }
		}
	    }      
	}
    }
  
  // Compute cell-centred z-components slopes in x- and y-direction
  for (int d = 0; d < AMREX_SPACEDIM ; d++)
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
      
      for (MFIter mfi(slopes[BZ_LOCAL], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Dim3 hiDomain = ubound(geom.Domain());

	  const Array4<Real> arr = S_source.array(mfi);
	  Array4<Real> slopesBZ = slopes[BZ_LOCAL].array(mfi);
	  Array4<Real> slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y-1; j <= hi.y+1; j++)
		{
		  for(int i = lo.x-1; i <= hi.x+1; i++)
		    {

		      Real u_iMinus1 = arr(i-iOffset,j-jOffset,k-kOffset,BZ);
		      Real u_i = arr(i,j,k,BZ);
		      Real u_iPlus1 = arr(i+iOffset,j+jOffset,k+kOffset,BZ);
		      slopesBZ(i,j,k,d) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		      u_iMinus1 = arr(i-iOffset,j-jOffset,k-kOffset,EZ);
		      u_i = arr(i,j,k,EZ);
		      u_iPlus1 = arr(i+iOffset,j+jOffset,k+kOffset,EZ);
		      slopesEZ(i,j,k,d) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
		    }
		}
	    }      
	}
    }
  
  // Compute y-components slopes in x-direction
  for (MFIter mfi(slopes[BY_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());
      
      const Array4<Real> arr = S_EM_source[1].array(mfi);
      Array4<Real> slopesBY = slopes[BY_LOCAL].array(mfi);
      Array4<Real> slopesEY = slopes[EY_LOCAL].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{

		  Real u_iMinus1 = arr(i-1,j,k,BY_LOCAL);
		  Real u_i = arr(i,j,k,BY_LOCAL);		    
		  Real u_iPlus1 = arr(i+1,j,k,BY_LOCAL);
		  slopesBY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		  u_iMinus1 = arr(i-1,j,k,EY_LOCAL);
		  u_i = arr(i,j,k,EY_LOCAL);		    
		  u_iPlus1 = arr(i+1,j,k,EY_LOCAL);
		  slopesEY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		}
	    }
	}      
    }
  // Compute x-components slopes in y-direction
  for (MFIter mfi(slopes[BX_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());

      const Array4<Real> arr = S_EM_source[0].array(mfi);
      Array4<Real> slopesBX = slopes[BX_LOCAL].array(mfi);
      Array4<Real> slopesEX = slopes[EX_LOCAL].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{

		  Real u_iMinus1 = arr(i,j-1,k,BX_LOCAL);
		  Real u_i = arr(i,j,k,BX_LOCAL);		    
		  Real u_iPlus1 = arr(i,j+1,k,BX_LOCAL);
		  slopesBX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		  u_iMinus1 = arr(i,j-1,k,EX_LOCAL);
		  u_i = arr(i,j,k,EX_LOCAL);		    
		  u_iPlus1 = arr(i,j+1,k,EX_LOCAL);
		  slopesEX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		}
	    }
	}      
    }


  // Coefficients for the magnetic and electric fields
  // For second order, there are 21 coeff.
  // i.e. a0,ax,ay,az,axx,...,b0,bx,...,c0,cx,...,czz
  // use 1 ghost cell
  int a0=0,ax=1,ay=2,az=3,axx=4,axy=5,axz=6;
  int b0=7,bx=8,by=9,bz=10,byy=11,bxy=12,byz=13;
  int c0=14,cx=15,cy=16,cz=17,czz=18,cxz=19,cyz=20;
  
  MultiFab Bcoeff(grids, dmap, 21, 1);
  MultiFab Ecoeff(grids, dmap, 21, 1);

  // Compute coefficients for the field function
  for (MFIter mfi(Bcoeff, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      // coeff
      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      // data
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi);       

      // slopes
      const auto& slopesBX = slopes[BX_LOCAL].array(mfi);
      const auto& slopesBY = slopes[BY_LOCAL].array(mfi);
      const auto& slopesBZ = slopes[BZ_LOCAL].array(mfi);
      const auto& slopesEX = slopes[EX_LOCAL].array(mfi);
      const auto& slopesEY = slopes[EY_LOCAL].array(mfi);
      const auto& slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{		 
		  
		  Bc(i,j,k,ay) = (slopesBX(i+1,j,k,Y)+slopesBX(i,j,k,Y))/2.0;
		  Bc(i,j,k,axy) = (slopesBX(i+1,j,k,Y)-slopesBX(i,j,k,Y));
		  Bc(i,j,k,az) = 0.0;
		  Bc(i,j,k,axz) = 0.0;
		  Bc(i,j,k,bx) = (slopesBY(i,j+1,k,X)+slopesBY(i,j,k,X))/2.0;
		  Bc(i,j,k,bxy) = (slopesBY(i,j+1,k,X)-slopesBY(i,j,k,X));
		  Bc(i,j,k,bz) = 0.0;
		  Bc(i,j,k,byz) = 0.0;
		  Bc(i,j,k,cx) = slopesBZ(i,j,k,X);
		  Bc(i,j,k,cxz) = 0.0;
		  Bc(i,j,k,cy) = slopesBZ(i,j,k,Y);
		  Bc(i,j,k,cyz) = 0.0;
		  Bc(i,j,k,axx) = -dx[0]/2.0*(Bc(i,j,k,bxy)/dx[1]); // +cxz/dx[2]
		  Bc(i,j,k,byy) = -dx[1]/2.0*(Bc(i,j,k,axy)/dx[0]); // ++cyz/dx[2]
		  Bc(i,j,k,czz) = 0.0;
		  Bc(i,j,k,a0) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - Bc(i,j,k,axx)/6.0;
		  Bc(i,j,k,ax) = (arrEM_X(i+1,j,k,BX_LOCAL)-arrEM_X(i,j,k,BX_LOCAL));
		  Bc(i,j,k,b0) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - Bc(i,j,k,byy)/6.0;
		  Bc(i,j,k,by) = (arrEM_Y(i,j+1,k,BY_LOCAL)-arrEM_Y(i,j,k,BY_LOCAL));
		  Bc(i,j,k,c0) = arr(i,j,k,BZ) - Bc(i,j,k,czz)/6.0;
		  Bc(i,j,k,cz) = 0.0;

		  Ec(i,j,k,ay) = (slopesEX(i+1,j,k,Y)+slopesEX(i,j,k,Y))/2.0;
		  Ec(i,j,k,axy) = (slopesEX(i+1,j,k,Y)-slopesEX(i,j,k,Y));
		  Ec(i,j,k,az) = 0.0;
		  Ec(i,j,k,axz) = 0.0;
		  Ec(i,j,k,bx) = (slopesEY(i,j+1,k,X)+slopesEY(i,j,k,X))/2.0;
		  Ec(i,j,k,bxy) = (slopesEY(i,j+1,k,X)-slopesEY(i,j,k,X));
		  Ec(i,j,k,bz) = 0.0;
		  Ec(i,j,k,byz) = 0.0;
		  Ec(i,j,k,cx) = slopesEZ(i,j,k,X);
		  Ec(i,j,k,cxz) = 0.0;
		  Ec(i,j,k,cy) = slopesEZ(i,j,k,Y);
		  Ec(i,j,k,cyz) = 0.0;
		  Ec(i,j,k,axx) = -dx[0]/2.0*(Ec(i,j,k,bxy)/dx[1]-slopesQ(i,j,k,0)); // +cxz/dx[2]
		  Ec(i,j,k,byy) = -dx[1]/2.0*(Ec(i,j,k,axy)/dx[0]-slopesQ(i,j,k,1)); // ++cyz/dx[2]
		  Ec(i,j,k,czz) = 0.0;
		  Ec(i,j,k,a0) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - Ec(i,j,k,axx)/6.0;
		  Ec(i,j,k,ax) = (arrEM_X(i+1,j,k,EX_LOCAL)-arrEM_X(i,j,k,EX_LOCAL));
		  Ec(i,j,k,b0) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - Ec(i,j,k,byy)/6.0;
		  Ec(i,j,k,by) = (arrEM_Y(i,j+1,k,EY_LOCAL)-arrEM_Y(i,j,k,EY_LOCAL));
		  Ec(i,j,k,c0) = arr(i,j,k,EZ) - Ec(i,j,k,czz)/6.0; //czz*dx[2]*dx[2]/6.0
		  Ec(i,j,k,cz) = 0.0;		  

		}
	    }
	}
    }

  // Compute edge components of the EM fields    
  for (MFIter mfi(fluxesEM, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      // use old array, S_dest should also work
      //const auto& arr = S_dest.array(mfi);
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi); 
      const auto& fluxArrEM = fluxesEM.array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);


      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{		 
		  
		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
		  // LD state
		  x = dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = -1, jOffset = -1;
		  Vector<Real> EM_LD = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  		  
		  // RD state
		  x = -dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_RD = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // LU state
		  x = dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_LU = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // RU state
		  x = -dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_RU = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // If the reconstruction is consistent, then EyRU=EyRD
		  Real EyR = (EM_RU[EY_LOCAL]+EM_RD[EY_LOCAL])/2.0;
		  Real EyL = (EM_LU[EY_LOCAL]+EM_LD[EY_LOCAL])/2.0;
		  Real ExU = (EM_RU[EX_LOCAL]+EM_LU[EX_LOCAL])/2.0;
		  Real ExD = (EM_RD[EX_LOCAL]+EM_LD[EX_LOCAL])/2.0;
		  fluxArrEM(i,j,k,BZ_LOCAL) = 0.25*(EM_RU[BZ_LOCAL]+EM_RD[BZ_LOCAL]+EM_LU[BZ_LOCAL]+EM_LD[BZ_LOCAL])
		    - 0.5/c*(EyR-EyL) + 0.5/c*(ExU-ExD);

		  Real ByR = (EM_RU[BY_LOCAL]+EM_RD[BY_LOCAL])/2.0;
		  Real ByL = (EM_LU[BY_LOCAL]+EM_LD[BY_LOCAL])/2.0;
		  Real BxU = (EM_RU[BX_LOCAL]+EM_LU[BX_LOCAL])/2.0;
		  Real BxD = (EM_RD[BX_LOCAL]+EM_LD[BX_LOCAL])/2.0;
		  fluxArrEM(i,j,k,EZ_LOCAL) = 0.25*(EM_RU[EZ_LOCAL]+EM_RD[EZ_LOCAL]+EM_LU[EZ_LOCAL]+EM_LD[EZ_LOCAL])
		    + 0.5*c*(ByR-ByL) - 0.5*c*(BxU-BxD);
		}
	    }
	}       
    }

  // Update face-centred EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      int d_EM = (d==0) ? 1 : 0;
    
      // Loop over all the patches at this level
      for (MFIter mfi(S_EM_source[d_EM], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();

	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  // Indexable arrays for the data, and the directional flux
	  // Based on the corner-centred definition of the flux array, the
	  // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	  //const auto& arr = S_dest.array(mfi);
	  const auto& arrEM = S_EM_dest[d_EM].array(mfi);
	  const auto& arrEMOld = S_EM_source[d_EM].array(mfi);
	  const auto& fluxArrEM = fluxesEM.array(mfi);
	  const auto& fluxArr = fluxes[d_EM].array(mfi);

	  const Dim3 hiDomain = ubound(geom.Domain());
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      // 2D code; z-component updated using cell-centred scheme
		      arrEM(i,j,k,BX_LOCAL+d_EM) = arrEMOld(i,j,k,BX_LOCAL+d_EM) + std::pow(-1,d)*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EZ_LOCAL) - fluxArrEM(i,j,k,EZ_LOCAL));
		      arrEM(i,j,k,EX_LOCAL+d_EM) = arrEMOld(i,j,k,EX_LOCAL+d_EM) - std::pow(-1,d)*c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) - fluxArrEM(i,j,k,BZ_LOCAL));

		      // 3D code
		      /*arrEM(i,j,k,BX_LOCAL+(1+d)%3) += (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(2+d)%3));
			arrEM(i,j,k,BX_LOCAL+(2+d)%3) -= (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(1+d)%3));
			arrEM(i,j,k,EX_LOCAL+(1+d)%3) -= c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(2+d)%3));
			arrEM(i,j,k,EX_LOCAL+(2+d)%3) += c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(1+d)%3));
		      */
		      
		      // source terms
		      if (sourceMethod!="IM")
			{
			  Real currentFace = r_i*fluxArr(i,j,k,RHO_I) + r_e*fluxArr(i,j,k,RHO_E);
			  arrEM(i,j,k,EX_LOCAL+d_EM) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;
			}
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

  // Compute cell-centred EM fields from face-centred
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrEM_X = S_EM_dest[0].array(mfi);
      const auto& arrEM_Y = S_EM_dest[1].array(mfi); 

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

  		  arr(i,j,k,BX) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0;// - Bc(i,j,k,axx)/6.0;
  		  arr(i,j,k,BY) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0;// - Bc(i,j,k,byy)/6.0;
  		  arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0;// - Ec(i,j,k,axx)/6.0;
  		  arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0;// - Ec(i,j,k,byy)/6.0;
		  
  		}
  	    }
  	}       
    }

  // Only 2D
  // Compute edge components for y-components at x-faces
  for (MFIter mfi(fluxes[0], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      const auto& fluxArrX = fluxes[0].array(mfi);
      //const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // L state
  		  x = dx[0]/2.0, y = 0.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_L = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		   
  		  // R state
  		  x = -dx[0]/2.0, y = 0.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_R = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // Note that this is not the flux, but the star states
  		  fluxArrX(i,j,k,BY) = 0.5*(EM_R[BY_LOCAL]+EM_L[BY_LOCAL])
		    + 0.5/c*(EM_R[EZ_LOCAL]-EM_L[EZ_LOCAL]);
  		  fluxArrX(i,j,k,EY) = 0.5*(EM_R[EY_LOCAL]+EM_L[EY_LOCAL])
		    - 0.5*c*(EM_R[BZ_LOCAL]-EM_L[BZ_LOCAL]);
  		}
  	    }
  	}
    }
  
  // Only 2D
  // Compute edge components for x-components at y-faces
  for (MFIter mfi(fluxes[1], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 		  

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;
		  	
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;

  		  // D state
  		  x = 0.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_D = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
  		  // U state
  		  x = 0.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_U = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // Note that this is not the flux, but the star states		  
  		  fluxArrY(i,j,k,BX) = 0.5*(EM_U[BX_LOCAL]+EM_D[BX_LOCAL])
		    - 0.5/c*(EM_U[EZ_LOCAL]-EM_D[EZ_LOCAL]);
  		  fluxArrY(i,j,k,EX) = 0.5*(EM_U[EX_LOCAL]+EM_D[EX_LOCAL])
		    + 0.5*c*(EM_U[BZ_LOCAL]-EM_D[BZ_LOCAL]); 
  		}
  	    }
  	}
    }

  // Update cell-centred z-components of EM fields
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrOld = S_source.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);
	        
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  // Update cell-centred z-components becuause it is 2D code
  		  arr(i,j,k,BZ) = arrOld(i,j,k,BZ) - (dt / dx[0]) * (fluxArrX(i+1,j,k,EY) - fluxArrX(i,j,k,EY)) + (dt / dx[1]) * (fluxArrY(i,j+1,k,EX) - fluxArrY(i,j,k,EX));
  		  arr(i,j,k,EZ) = arrOld(i,j,k,EZ) + c*c*(dt / dx[0]) * (fluxArrX(i+1,j,k,BY) - fluxArrX(i,j,k,BY)) - c*c*(dt / dx[1]) * (fluxArrY(i,j+1,k,BX) - fluxArrY(i,j,k,BX));		  
		  
  		  // source terms
		  if (sourceMethod!="IM")
		    {
		      arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));

		      //Vector<Real> u_i = get_data_zone(arrOld,i,j,k,0,NUM_STATE_FLUID/2);
		      //Vector<Real> u_e = get_data_zone(arrOld,i,j,k,RHO_E,NUM_STATE_FLUID/2);
		      //Vector<Real> u_i_HLLC = HLLC(u_i,u_i,2);
		      //Vector<Real> u_e_HLLC = HLLC(u_e,u_e,2);
		      //Real currentFace = r_i*u_i_HLLC[0] + r_e*u_e_HLLC[0];
		      //arr(i,j,k,EZ) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;
		    }
  		}
  	    }
  	}
    }
  // We need to compute boundary conditions again after each update
  S_dest.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_dest, geom, bc);  
}
void CAMReXmp::MaxwellSolverDivFreeWENOcharacteristic(Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& fluxesEM, MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{

  // Note that fluxesEM and fluxes are different definitions
  // fluxesEM are EM star states, fluxes are fluxes
  
  // Useful when dealing with high-order reconstruction
  // To get first-order slope in y-direction for Bx: slopes[BX_LOCALL].array(i,j,k,Y)
  // Useful to define indices int X,Y,Z in second-order linear reconstruction
  // For third-order quadratic reconstruction use int X,Y,Z,XX,YY,ZZ,XY,YZ,XZ
  // Although for x-components fields, only exist Y,Z,YY,ZZ,YZ (see Balsara et al. 2017)
  
  Array<MultiFab,6> slopes;
  // number of slopes
  // 9 slopes for third-order parabolic reconstruction
  // Ex = Ex^0 + Ex^y*(y/Delta_y) + Ex^z*(z/Delta_z) + Ex^yy*((y/Delta_y)^2-1/12) + ...
  // where Ex^y, Ex^z, Ex^yy are the slopes
  // For 2D code, only 2 slopes for x- and y-compoenents
  // For clarity, will use 9 slopes for third-order, first one represents slope in x-direction
  int nSlopes = 9;
  // indices for the slopes
  // for 2D, do not need Z indices
  int X=0,Y=1,Z=2,XX=3,YY=4,ZZ=5,XY=6,YZ=7,XZ=8;
  slopes[BX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);  
  slopes[EX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);
#if (AMREX_SPACEDIM >= 2)
  slopes[BY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
#else  
  // For 1D code it is cell-centred
  slopes[BY_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
  
#if (AMREX_SPACEDIM != 3)  
  // For 2D code it is cell-centred
  slopes[BZ_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EZ_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
  
  // Three components, one for each direction
  MultiFab slopesCharge(grids, dmap, nSlopes, 1);

  // Compute slopes of the charge densities  
  for (MFIter mfi(slopesCharge, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());
      
      // uses old charge
      const Array4<Real> arr = S_source.array(mfi);
      Array4<Real> slopes = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  Vector<Real> data_i = get_data_stencil(arr, i, j, k, 2, 0, 0, RHO_I);
		  Vector<Real> data_e = get_data_stencil(arr, i, j, k, 2, 0, 0, RHO_E);
		  Vector<Real> data_charge = get_charge_scaled(data_i,data_e);
		  std::array<Real, 2> slopesX = WENO3_slope(data_charge);
		  slopes(i,j,k,X) = slopesX[0];
		  slopes(i,j,k,XX) = slopesX[1];

		  data_i = get_data_stencil(arr, i, j, k, 0, 2, 0, RHO_I);
		  data_e = get_data_stencil(arr, i, j, k, 0, 2, 0, RHO_E);
		  data_charge = get_charge_scaled(data_i,data_e);
		  std::array<Real, 2> slopesY = WENO3_slope(data_charge);
		  slopes(i,j,k,Y) = slopesY[0];
		  slopes(i,j,k,YY) = slopesY[1];
		  
		  data_i = get_data_stencil(arr, i, j, k, 1, 1, 0, RHO_I);
		  data_e = get_data_stencil(arr, i, j, k, 1, 1, 0, RHO_E);
		  data_charge = get_charge_scaled(data_i,data_e);
		  Real slopesCross = WENO3_slopeCross(data_charge,
						      {slopes(i,j,k,X),slopes(i,j,k,XX),slopes(i,j,k,Y),slopes(i,j,k,YY)});
		  slopes(i,j,k,XY) = slopesCross;
		  if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
		    {
		      int iOffset = 1, jOffset = 0, kOffset = 0;
		      Real u_iMinus1 = (r_i*arr(i-iOffset,j-jOffset,k-kOffset,RHO_I)
					+ r_e*arr(i-iOffset,j-jOffset,k-kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_iPlus1 = (r_i*arr(i+iOffset,j+jOffset,k+kOffset,RHO_I)
				       + r_e*arr(i+iOffset,j+jOffset,k+kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      slopes(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
		      slopes(i,j,k,XX) = 0.0, slopes(i,j,k,XY) = 0.0; 
		    }
		  if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
		    {
		      int iOffset = 0, jOffset = 1, kOffset = 0;
		      Real u_iMinus1 = (r_i*arr(i-iOffset,j-jOffset,k-kOffset,RHO_I)
					+ r_e*arr(i-iOffset,j-jOffset,k-kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_iPlus1 = (r_i*arr(i+iOffset,j+jOffset,k+kOffset,RHO_I)
				       + r_e*arr(i+iOffset,j+jOffset,k+kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      slopes(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
		      slopes(i,j,k,YY) = 0.0, slopes(i,j,k,XY) = 0.0; 			  
		    }
		}
	    }
	}      
    }
  
  // Compute cell-centred z-components slopes in x- and y-direction
  //for (int d = 0; d < AMREX_SPACEDIM ; d++)
    {

      //const int iOffset = ( d == 0 ? 1 : 0);
      //const int jOffset = ( d == 1 ? 1 : 0);
      //const int kOffset = ( d == 2 ? 1 : 0);
      
      for (MFIter mfi(slopes[BZ_LOCAL], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Dim3 hiDomain = ubound(geom.Domain());
	  
	  const Array4<Real> arr = S_source.array(mfi);
	  Array4<Real> slopesBZ = slopes[BZ_LOCAL].array(mfi);
	  Array4<Real> slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y-1; j <= hi.y+1; j++)
		{
		  for(int i = lo.x-1; i <= hi.x+1; i++)
		    {
		      // x slopes
		      Vector<Vector<Real>> leftx, rightx;
		      leftx.push_back({0.5*c,0.0,0.0,0.5});
		      leftx.push_back({0.0,-0.5*c,0.5,0.0});
		      leftx.push_back({-0.5*c,0.0,0.0,0.5});
		      leftx.push_back({0.0,0.5*c,0.5,0.0});

		      rightx.push_back({1.0/c,0.0,-1.0/c,0.0});
		      rightx.push_back({0.0,-1.0/c,0.0,1.0/c});
		      rightx.push_back({0.0,1.0,0.0,1.0});
		      rightx.push_back({1.0,0.0,1.0,0.0});

		      // Conserved variables covering the stencils
		      Vector<Real> u_iMinus2 = {arr(i-2,j,k,BY),arr(i-2,j,k,BZ),arr(i-2,j,k,EY),arr(i-2,j,k,EZ)};
		      Vector<Real> u_iMinus1 = {arr(i-1,j,k,BY),arr(i-1,j,k,BZ),arr(i-1,j,k,EY),arr(i-1,j,k,EZ)};
		      Vector<Real> u_i = {arr(i,j,k,BY),arr(i,j,k,BZ),arr(i,j,k,EY),arr(i,j,k,EZ)};
		      Vector<Real> u_iPlus1 = {arr(i+1,j,k,BY),arr(i+1,j,k,BZ),arr(i+1,j,k,EY),arr(i+1,j,k,EZ)};
		      Vector<Real> u_iPlus2 = {arr(i+2,j,k,BY),arr(i+2,j,k,BZ),arr(i+2,j,k,EY),arr(i+2,j,k,EZ)};
		      // Characteristic variables covering the stencils
		      Vector<Real> w_iMinus2,w_iMinus1,w_i,w_iPlus1,w_iPlus2;
		      for (int n=0; n<leftx.size(); n++)
			{
			  w_iMinus2.push_back(dotProduct(leftx[n],u_iMinus2));
			  w_iMinus1.push_back(dotProduct(leftx[n],u_iMinus1));
			  w_i.push_back(dotProduct(leftx[n],u_i));
			  w_iPlus1.push_back(dotProduct(leftx[n],u_iPlus1));
			  w_iPlus2.push_back(dotProduct(leftx[n],u_iPlus2));
			}

		      // Stencils of characteristic variables
		      Vector<Real> w1 = {w_iMinus2[0],w_iMinus1[0],w_i[0],
					 w_iPlus1[0],w_iPlus2[0]};
		      Vector<Real> w2 = {w_iMinus2[1],w_iMinus1[1],w_i[1],
					 w_iPlus1[1],w_iPlus2[1]};
		      Vector<Real> w3 = {w_iMinus2[2],w_iMinus1[2],w_i[2],
					 w_iPlus1[2],w_iPlus2[2]};
		      Vector<Real> w4 = {w_iMinus2[3],w_iMinus1[3],w_i[3],
					 w_iPlus1[3],w_iPlus2[3]};		  
		  
		      std::array<Real, 2> slopesw1 = WENO3_slope(w1);		      
		      std::array<Real, 2> slopesw2 = WENO3_slope(w2);
		      std::array<Real, 2> slopesw3 = WENO3_slope(w3);
		      std::array<Real, 2> slopesw4 = WENO3_slope(w4);

		      slopesBZ(i,j,k,X) = dotProduct(rightx[1],{slopesw1[0],slopesw2[0],slopesw3[0],slopesw4[0]});
		      slopesBZ(i,j,k,XX) = dotProduct(rightx[1],{slopesw1[1],slopesw2[1],slopesw3[1],slopesw4[1]});

		      slopesEZ(i,j,k,X) = dotProduct(rightx[3],{slopesw1[0],slopesw2[0],slopesw3[0],slopesw4[0]});
		      slopesEZ(i,j,k,XX) = dotProduct(rightx[3],{slopesw1[1],slopesw2[1],slopesw3[1],slopesw4[1]});

		      // y slopes
		      Vector<Vector<Real>> lefty, righty;
		      lefty.push_back({-0.5*c,0.0,0.0,0.5});
		      lefty.push_back({0.0,0.5*c,0.5,0.0});
		      lefty.push_back({0.5*c,0.0,0.0,0.5});
		      lefty.push_back({0.0,-0.5*c,0.5,0.0});

		      righty.push_back({-1.0/c,0.0,1.0/c,0.0});
		      righty.push_back({0.0,1.0/c,0.0,-1.0/c});
		      righty.push_back({0.0,1.0,0.0,1.0});
		      righty.push_back({1.0,0.0,1.0,0.0});

		      // Conserved variables covering the stencils
		      u_iMinus2 = {arr(i,j-2,k,BX),arr(i,j-2,k,BZ),arr(i,j-2,k,EX),arr(i,j-2,k,EZ)};
		      u_iMinus1 = {arr(i,j-1,k,BX),arr(i,j-1,k,BZ),arr(i,j-1,k,EX),arr(i,j-1,k,EZ)};
		      u_i = {arr(i,j,k,BX),arr(i,j,k,BZ),arr(i,j,k,EX),arr(i,j,k,EZ)};
		      u_iPlus1 = {arr(i,j+1,k,BX),arr(i,j+1,k,BZ),arr(i,j+1,k,EX),arr(i,j+1,k,EZ)};
		      u_iPlus2 = {arr(i,j+2,k,BX),arr(i,j+2,k,BZ),arr(i,j+2,k,EX),arr(i,j+2,k,EZ)};
		      // Characteristic variables covering the stencils
		      for (int n=0; n<lefty.size(); n++)
			{
			  w_iMinus2[n] = dotProduct(lefty[n],u_iMinus2);
			  w_iMinus1[n] = dotProduct(lefty[n],u_iMinus1);
			  w_i[n] = dotProduct(lefty[n],u_i);
			  w_iPlus1[n] = dotProduct(lefty[n],u_iPlus1);
			  w_iPlus2[n] = dotProduct(lefty[n],u_iPlus2);
			}

		      // Stencils of characteristic variables
		      w1 = {w_iMinus2[0],w_iMinus1[0],w_i[0],
			    w_iPlus1[0],w_iPlus2[0]};
		      w2 = {w_iMinus2[1],w_iMinus1[1],w_i[1],
			    w_iPlus1[1],w_iPlus2[1]};
		      w3 = {w_iMinus2[2],w_iMinus1[2],w_i[2],
			    w_iPlus1[2],w_iPlus2[2]};
		      w4 = {w_iMinus2[3],w_iMinus1[3],w_i[3],
			    w_iPlus1[3],w_iPlus2[3]};		  
		  
		      slopesw1 = WENO3_slope(w1);		      
		      slopesw2 = WENO3_slope(w2);
		      slopesw3 = WENO3_slope(w3);
		      slopesw4 = WENO3_slope(w4);

		      slopesBZ(i,j,k,Y) = dotProduct(righty[1],{slopesw1[0],slopesw2[0],slopesw3[0],slopesw4[0]});
		      slopesBZ(i,j,k,YY) = dotProduct(righty[1],{slopesw1[1],slopesw2[1],slopesw3[1],slopesw4[1]});
		     		      
		      slopesEZ(i,j,k,Y) = dotProduct(righty[3],{slopesw1[0],slopesw2[0],slopesw3[0],slopesw4[0]});
		      slopesEZ(i,j,k,YY) = dotProduct(righty[3],{slopesw1[1],slopesw2[1],slopesw3[1],slopesw4[1]});

		      // cross slopes
		      Vector<Real> dataXY = get_data_stencil(arr, i, j, k, 1, 1, 0, BZ);
		      Real slopesCross = WENO3_slopeCross(dataXY,
							  {slopesBZ(i,j,k,X),slopesBZ(i,j,k,XX),slopesBZ(i,j,k,Y),slopesBZ(i,j,k,YY)});
		      slopesBZ(i,j,k,XY) = slopesCross;
		      
		      dataXY = get_data_stencil(arr, i, j, k, 1, 1, 0, EZ);
		      slopesCross = WENO3_slopeCross(dataXY,
						     {slopesEZ(i,j,k,X),slopesEZ(i,j,k,XX),slopesEZ(i,j,k,Y),slopesEZ(i,j,k,YY)});
		      slopesEZ(i,j,k,XY) = slopesCross;
		      
		      if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
			{			      
			  {
			    Real u_iMinus1 = arr(i-1,j,k,BZ);
			    Real u_i = arr(i,j,k,BZ);
			    Real u_iPlus1 = arr(i+1,j,k,BZ);
			    slopesBZ(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesBZ(i,j,k,XX) = 0.0, slopesBZ(i,j,k,XY) = 0.0;
			  }
			  {				  
			    Real u_iMinus1 = arr(i-1,j,k,EZ);
			    Real u_i = arr(i,j,k,EZ);
			    Real u_iPlus1 = arr(i+1,j,k,EZ);
			    slopesEZ(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesEZ(i,j,k,XX) = 0.0, slopesEZ(i,j,k,XY) = 0.0;
			  }			      
			}
		      if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
			{
			  {			
			    Real u_iMinus1 = arr(i,j-1,k,BZ);
			    Real u_i = arr(i,j,k,BZ);
			    Real u_iPlus1 = arr(i,j+1,k,BZ);
			    slopesBZ(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesBZ(i,j,k,YY) = 0.0, slopesBZ(i,j,k,XY) = 0.0;
			    
			  }			  
			  {
			    Real u_iMinus1 = arr(i,j-1,k,EZ);
			    Real u_i = arr(i,j,k,EZ);
			    Real u_iPlus1 = arr(i,j+1,k,EZ);
			    slopesEZ(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesEZ(i,j,k,YY) = 0.0, slopesEZ(i,j,k,XY) = 0.0;
			  }			      
			}
		    }
		}
	    }      
	}
    }    

  // Compute y-components slopes in x-direction
  for (MFIter mfi(slopes[BY_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());

      const Array4<Real> arr = S_EM_source[1].array(mfi);
      Array4<Real> slopesBY = slopes[BY_LOCAL].array(mfi);
      Array4<Real> slopesEY = slopes[EY_LOCAL].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  // x slopes
		  Vector<Vector<Real>> leftx, rightx;
		  leftx.push_back({0.5*c,0.0,0.0,0.5});
		  leftx.push_back({0.0,-0.5*c,0.5,0.0});
		  leftx.push_back({-0.5*c,0.0,0.0,0.5});
		  leftx.push_back({0.0,0.5*c,0.5,0.0});

		  rightx.push_back({1.0/c,0.0,-1.0/c,0.0});
		  rightx.push_back({0.0,-1.0/c,0.0,1.0/c});
		  rightx.push_back({0.0,1.0,0.0,1.0});
		  rightx.push_back({1.0,0.0,1.0,0.0});

		  // Conserved variables covering the stencils
		  Vector<Real> u_iMinus2 = {arr(i-2,j,k,BY_LOCAL),arr(i-2,j,k,BZ_LOCAL),arr(i-2,j,k,EY_LOCAL),arr(i-2,j,k,EZ_LOCAL)};
		  Vector<Real> u_iMinus1 = {arr(i-1,j,k,BY_LOCAL),arr(i-1,j,k,BZ_LOCAL),arr(i-1,j,k,EY_LOCAL),arr(i-1,j,k,EZ_LOCAL)};
		  Vector<Real> u_i = {arr(i,j,k,BY_LOCAL),arr(i,j,k,BZ_LOCAL),arr(i,j,k,EY_LOCAL),arr(i,j,k,EZ_LOCAL)};
		  Vector<Real> u_iPlus1 = {arr(i+1,j,k,BY_LOCAL),arr(i+1,j,k,BZ_LOCAL),arr(i+1,j,k,EY_LOCAL),arr(i+1,j,k,EZ_LOCAL)};
		  Vector<Real> u_iPlus2 = {arr(i+2,j,k,BY_LOCAL),arr(i+2,j,k,BZ_LOCAL),arr(i+2,j,k,EY_LOCAL),arr(i+2,j,k,EZ_LOCAL)};
		  // Characteristic variables covering the stencils
		  Vector<Real> w_iMinus2,w_iMinus1,w_i,w_iPlus1,w_iPlus2;
		  for (int n=0; n<leftx.size(); n++)
		    {
		      w_iMinus2.push_back(dotProduct(leftx[n],u_iMinus2));
		      w_iMinus1.push_back(dotProduct(leftx[n],u_iMinus1));
		      w_i.push_back(dotProduct(leftx[n],u_i));
		      w_iPlus1.push_back(dotProduct(leftx[n],u_iPlus1));
		      w_iPlus2.push_back(dotProduct(leftx[n],u_iPlus2));
		    }

		  // Stencils of characteristic variables
		  Vector<Real> w1 = {w_iMinus2[0],w_iMinus1[0],w_i[0],
				     w_iPlus1[0],w_iPlus2[0]};
		  Vector<Real> w2 = {w_iMinus2[1],w_iMinus1[1],w_i[1],
				     w_iPlus1[1],w_iPlus2[1]};
		  Vector<Real> w3 = {w_iMinus2[2],w_iMinus1[2],w_i[2],
				     w_iPlus1[2],w_iPlus2[2]};
		  Vector<Real> w4 = {w_iMinus2[3],w_iMinus1[3],w_i[3],
				     w_iPlus1[3],w_iPlus2[3]};		  
		  
		  std::array<Real, 2> slopesw1 = WENO3_slope(w1);		      
		  std::array<Real, 2> slopesw2 = WENO3_slope(w2);
		  std::array<Real, 2> slopesw3 = WENO3_slope(w3);
		  std::array<Real, 2> slopesw4 = WENO3_slope(w4);

		  slopesBY(i,j,k,X) = dotProduct(rightx[0],{slopesw1[0],slopesw2[0],slopesw3[0],slopesw4[0]});
		  slopesBY(i,j,k,XX) = dotProduct(rightx[0],{slopesw1[1],slopesw2[1],slopesw3[1],slopesw4[1]});

		  slopesEY(i,j,k,X) = dotProduct(rightx[2],{slopesw1[0],slopesw2[0],slopesw3[0],slopesw4[0]});
		  slopesEY(i,j,k,XX) = dotProduct(rightx[2],{slopesw1[1],slopesw2[1],slopesw3[1],slopesw4[1]});

		  if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
		    {			      
		      {
			Real u_iMinus1 = arr(i-1,j,k,BY_LOCAL);
			Real u_i = arr(i,j,k,BY_LOCAL);		    
			Real u_iPlus1 = arr(i+1,j,k,BY_LOCAL);
			slopesBY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesBY(i,j,k,XX) = 0.0;
		      }
		      {
			Real u_iMinus1 = arr(i-1,j,k,EY_LOCAL);
			Real u_i = arr(i,j,k,EY_LOCAL);		    
			Real u_iPlus1 = arr(i+1,j,k,EY_LOCAL);
			slopesEY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesEY(i,j,k,XX) = 0.0;
		      }			      
		    }		  
		}
	    }
	}      
    }
  // Compute x-components slopes in y-direction
  for (MFIter mfi(slopes[BX_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());

      const Array4<Real> arr = S_EM_source[0].array(mfi);
      Array4<Real> slopesBX = slopes[BX_LOCAL].array(mfi);
      Array4<Real> slopesEX = slopes[EX_LOCAL].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  Vector<Vector<Real>> lefty, righty;
		  lefty.push_back({-0.5*c,0.0,0.0,0.5});
		  lefty.push_back({0.0,0.5*c,0.5,0.0});
		  lefty.push_back({0.5*c,0.0,0.0,0.5});
		  lefty.push_back({0.0,-0.5*c,0.5,0.0});

		  righty.push_back({-1.0/c,0.0,1.0/c,0.0});
		  righty.push_back({0.0,1.0/c,0.0,-1.0/c});
		  righty.push_back({0.0,1.0,0.0,1.0});
		  righty.push_back({1.0,0.0,1.0,0.0});

		  // Conserved variables covering the stencils
		  Vector<Real> u_iMinus2 = {arr(i,j-2,k,BX_LOCAL),arr(i,j-2,k,BZ_LOCAL),arr(i,j-2,k,EX_LOCAL),arr(i,j-2,k,EZ_LOCAL)};
		  Vector<Real> u_iMinus1 = {arr(i,j-1,k,BX_LOCAL),arr(i,j-1,k,BZ_LOCAL),arr(i,j-1,k,EX_LOCAL),arr(i,j-1,k,EZ_LOCAL)};
		  Vector<Real> u_i = {arr(i,j,k,BX_LOCAL),arr(i,j,k,BZ_LOCAL),arr(i,j,k,EX_LOCAL),arr(i,j,k,EZ_LOCAL)};
		  Vector<Real> u_iPlus1 = {arr(i,j+1,k,BX_LOCAL),arr(i,j+1,k,BZ_LOCAL),arr(i,j+1,k,EX_LOCAL),arr(i,j+1,k,EZ_LOCAL)};
		  Vector<Real> u_iPlus2 = {arr(i,j+2,k,BX_LOCAL),arr(i,j+2,k,BZ_LOCAL),arr(i,j+2,k,EX_LOCAL),arr(i,j+2,k,EZ_LOCAL)};
		  // Characteristic variables covering the stencils
		  Vector<Real> w_iMinus2,w_iMinus1,w_i,w_iPlus1,w_iPlus2;
		  for (int n=0; n<lefty.size(); n++)
		    {
		      w_iMinus2.push_back(dotProduct(lefty[n],u_iMinus2));
		      w_iMinus1.push_back(dotProduct(lefty[n],u_iMinus1));
		      w_i.push_back(dotProduct(lefty[n],u_i));
		      w_iPlus1.push_back(dotProduct(lefty[n],u_iPlus1));
		      w_iPlus2.push_back(dotProduct(lefty[n],u_iPlus2));
		    }

		  // Stencils of characteristic variables
		  Vector<Real> w1 = {w_iMinus2[0],w_iMinus1[0],w_i[0],
				     w_iPlus1[0],w_iPlus2[0]};
		  Vector<Real> w2 = {w_iMinus2[1],w_iMinus1[1],w_i[1],
				     w_iPlus1[1],w_iPlus2[1]};
		  Vector<Real> w3 = {w_iMinus2[2],w_iMinus1[2],w_i[2],
				     w_iPlus1[2],w_iPlus2[2]};
		  Vector<Real> w4 = {w_iMinus2[3],w_iMinus1[3],w_i[3],
				     w_iPlus1[3],w_iPlus2[3]};		  
		  
		  std::array<Real, 2> slopesw1 = WENO3_slope(w1);		      
		  std::array<Real, 2> slopesw2 = WENO3_slope(w2);
		  std::array<Real, 2> slopesw3 = WENO3_slope(w3);
		  std::array<Real, 2> slopesw4 = WENO3_slope(w4);

		  slopesBX(i,j,k,Y) = dotProduct(righty[0],{slopesw1[0],slopesw2[0],slopesw3[0],slopesw4[0]});
		  slopesBX(i,j,k,YY) = dotProduct(righty[0],{slopesw1[1],slopesw2[1],slopesw3[1],slopesw4[1]});
		     		      
		  slopesEX(i,j,k,Y) = dotProduct(righty[2],{slopesw1[0],slopesw2[0],slopesw3[0],slopesw4[0]});
		  slopesEX(i,j,k,YY) = dotProduct(righty[2],{slopesw1[1],slopesw2[1],slopesw3[1],slopesw4[1]});

		  if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
		    {
		      {
			Real u_iMinus1 = arr(i,j-1,k,BX_LOCAL);
			Real u_i = arr(i,j,k,BX_LOCAL);		    
			Real u_iPlus1 = arr(i,j+1,k,BX_LOCAL);
			slopesBX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesBX(i,j,k,YY) = 0.0;			      
		      }
		      {
			Real u_iMinus1 = arr(i,j-1,k,EX_LOCAL);
			Real u_i = arr(i,j,k,EX_LOCAL);		    
			Real u_iPlus1 = arr(i,j+1,k,EX_LOCAL);
			slopesEX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesEX(i,j,k,YY) = 0.0;			      
		      }			      
		    }
		}
	    }
	}      
    }


  // Coefficients for the magnetic and electric fields
  // For third order, there are 48 coeff.
  // i.e. a0,ax,ay,az,axx,...,b0,bx,...,c0,cx,...,czz
  // use 1 ghost cell
  int a0=0,ax=1,ay=2,az=3,axx=4,ayy=5,azz=6,axy=7,ayz=8,axz=9,axxx=10,axxy=11,axxz=12,axyy=13,axzz=14,axyz=15;
  int b0=16,bx=17,by=18,bz=19,bxx=20,byy=21,bzz=22,bxy=23,byz=24,bxz=25,byyy=26,bxyy=27,byyz=28,bxxy=29,byzz=30,bxyz=31;
  int c0=32,cx=33,cy=34,cz=35,cxx=36,cyy=37,czz=38,cxy=39,cyz=40,cxz=41,czzz=42,cxzz=43,cyzz=44,cxxz=45,cyyz=46,cxyz=47;
  MultiFab Bcoeff(grids, dmap, 48, 1);
  MultiFab Ecoeff(grids, dmap, 48, 1);

  // Compute coefficients for the field function
  for (MFIter mfi(Bcoeff, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      // coeff
      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      // data
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi);       

      // slopes
      const auto& slopesBX = slopes[BX_LOCAL].array(mfi);
      const auto& slopesBY = slopes[BY_LOCAL].array(mfi);
      const auto& slopesBZ = slopes[BZ_LOCAL].array(mfi);
      const auto& slopesEX = slopes[EX_LOCAL].array(mfi);
      const auto& slopesEY = slopes[EY_LOCAL].array(mfi);
      const auto& slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  // third order terms
		  Bc(i,j,k,ayy) = (slopesBX(i+1,j,k,YY)+slopesBX(i,j,k,YY))/2.0;
		  Bc(i,j,k,axyy) = slopesBX(i+1,j,k,YY)-slopesBX(i,j,k,YY);
		  Bc(i,j,k,azz) = 0.0;
		  Bc(i,j,k,axzz) = 0.0;
		  Bc(i,j,k,ayz) = 0.0;
		  Bc(i,j,k,axyz) = 0.0;
		  Bc(i,j,k,bxx) = (slopesBY(i,j+1,k,XX)+slopesBY(i,j,k,XX))/2.0;
		  Bc(i,j,k,bxxy) = slopesBY(i,j+1,k,XX)-slopesBY(i,j,k,XX);
		  Bc(i,j,k,bzz) = 0.0;
		  Bc(i,j,k,byzz) = 0.0;
		  Bc(i,j,k,bxz) = 0.0;
		  Bc(i,j,k,bxyz) = 0.0;
		  Bc(i,j,k,cxx) = slopesBZ(i,j,k,XX);
		  Bc(i,j,k,cxxz) = 0.0;
		  Bc(i,j,k,cyy) = slopesBZ(i,j,k,YY);
		  Bc(i,j,k,cyyz) = 0.0;
		  Bc(i,j,k,cxy) = slopesBZ(i,j,k,XY);
		  Bc(i,j,k,cxyz) = 0.0;
		  Bc(i,j,k,axxx) = -dx[0]/3.0*(Bc(i,j,k,bxxy)/dx[1]);
		  Bc(i,j,k,byyy) = -dx[1]/3.0*(Bc(i,j,k,axyy)/dx[0]);
		  Bc(i,j,k,czzz) = 0.0;
		  Bc(i,j,k,axxy) = 0.0;
		  Bc(i,j,k,bxyy) = 0.0;
		  Bc(i,j,k,byyz) = 0.0;
		  Bc(i,j,k,cyzz) = 0.0;
		  Bc(i,j,k,cxzz) = 0.0;
		  Bc(i,j,k,axxz) = 0.0;
		  // second order terms
		  Bc(i,j,k,ay) = (slopesBX(i+1,j,k,Y)+slopesBX(i,j,k,Y))/2.0 - Bc(i,j,k,axxy)/6.0;
		  Bc(i,j,k,axy) = (slopesBX(i+1,j,k,Y)-slopesBX(i,j,k,Y));
		  Bc(i,j,k,az) = 0.0;
		  Bc(i,j,k,axz) = 0.0;
		  Bc(i,j,k,bx) = (slopesBY(i,j+1,k,X)+slopesBY(i,j,k,X))/2.0 - Bc(i,j,k,bxyy)/6.0;
		  Bc(i,j,k,bxy) = (slopesBY(i,j+1,k,X)-slopesBY(i,j,k,X));
		  Bc(i,j,k,bz) = 0.0;
		  Bc(i,j,k,byz) = 0.0;
		  Bc(i,j,k,cx) = slopesBZ(i,j,k,X) - Bc(i,j,k,cxzz)/6.0;
		  Bc(i,j,k,cxz) = 0.0;
		  Bc(i,j,k,cy) = slopesBZ(i,j,k,Y) - Bc(i,j,k,cyzz)/6.0;
		  Bc(i,j,k,cyz) = 0.0;
		  Bc(i,j,k,axx) = -dx[0]/2.0*(Bc(i,j,k,bxy)/dx[1]); // +cxz/dx[2]
		  Bc(i,j,k,byy) = -dx[1]/2.0*(Bc(i,j,k,axy)/dx[0]); // ++cyz/dx[2]
		  Bc(i,j,k,czz) = 0.0;
		  Bc(i,j,k,a0) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - Bc(i,j,k,axx)/6.0;
		  Bc(i,j,k,ax) = (arrEM_X(i+1,j,k,BX_LOCAL)-arrEM_X(i,j,k,BX_LOCAL)) - Bc(i,j,k,axxx)/10.0;
		  Bc(i,j,k,b0) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - Bc(i,j,k,byy)/6.0;
		  Bc(i,j,k,by) = (arrEM_Y(i,j+1,k,BY_LOCAL)-arrEM_Y(i,j,k,BY_LOCAL)) - Bc(i,j,k,byyy)/10.0;
		  Bc(i,j,k,c0) = arr(i,j,k,BZ) - Bc(i,j,k,czz)/6.0;
		  Bc(i,j,k,cz) = 0.0;

		  // third order terms
		  Ec(i,j,k,ayy) = (slopesEX(i+1,j,k,YY)+slopesEX(i,j,k,YY))/2.0;
		  Ec(i,j,k,axyy) = slopesEX(i+1,j,k,YY)-slopesEX(i,j,k,YY);
		  Ec(i,j,k,azz) = 0.0;
		  Ec(i,j,k,axzz) = 0.0;
		  Ec(i,j,k,ayz) = 0.0;
		  Ec(i,j,k,axyz) = 0.0;
		  Ec(i,j,k,bxx) = (slopesEY(i,j+1,k,XX)+slopesEY(i,j,k,XX))/2.0;
		  Ec(i,j,k,bxxy) = slopesEY(i,j+1,k,XX)-slopesEY(i,j,k,XX);
		  Ec(i,j,k,bzz) = 0.0;
		  Ec(i,j,k,byzz) = 0.0;
		  Ec(i,j,k,bxz) = 0.0;
		  Ec(i,j,k,bxyz) = 0.0;
		  Ec(i,j,k,cxx) = slopesEZ(i,j,k,XX);
		  Ec(i,j,k,cxxz) = 0.0;
		  Ec(i,j,k,cyy) = slopesEZ(i,j,k,YY);
		  Ec(i,j,k,cyyz) = 0.0;
		  Ec(i,j,k,cxy) = slopesEZ(i,j,k,XY);
		  Ec(i,j,k,cxyz) = 0.0;
		  Ec(i,j,k,axxx) = -dx[0]/3.0*(Ec(i,j,k,bxxy)/dx[1]-slopesQ(i,j,k,XX));
		  Ec(i,j,k,byyy) = -dx[1]/3.0*(Ec(i,j,k,axyy)/dx[0]-slopesQ(i,j,k,YY));
		  Ec(i,j,k,czzz) = 0.0;
		  Ec(i,j,k,axxy) = dx[0]*dx[1]/(2.0*(dx[0]+dx[1]))*slopesQ(i,j,k,XY);
		  Ec(i,j,k,bxyy) = Ec(i,j,k,axxy);
		  Ec(i,j,k,byyz) = 0.0;
		  Ec(i,j,k,cyzz) = 0.0;
		  Ec(i,j,k,cxzz) = 0.0;
		  Ec(i,j,k,axxz) = 0.0;
		  // second order terms
		  Ec(i,j,k,ay) = (slopesEX(i+1,j,k,Y)+slopesEX(i,j,k,Y))/2.0 - Ec(i,j,k,axxy)/6.0;
		  Ec(i,j,k,axy) = (slopesEX(i+1,j,k,Y)-slopesEX(i,j,k,Y));
		  Ec(i,j,k,az) = 0.0;
		  Ec(i,j,k,axz) = 0.0;
		  Ec(i,j,k,bx) = (slopesEY(i,j+1,k,X)+slopesEY(i,j,k,X))/2.0 - Ec(i,j,k,bxyy)/6.0;
		  Ec(i,j,k,bxy) = (slopesEY(i,j+1,k,X)-slopesEY(i,j,k,X));
		  Ec(i,j,k,bz) = 0.0;
		  Ec(i,j,k,byz) = 0.0;
		  Ec(i,j,k,cx) = slopesEZ(i,j,k,X) - Ec(i,j,k,cxzz)/6.0;
		  Ec(i,j,k,cxz) = 0.0;
		  Ec(i,j,k,cy) = slopesEZ(i,j,k,Y) - Ec(i,j,k,cyzz)/6.0;
		  Ec(i,j,k,cyz) = 0.0;
		  Ec(i,j,k,axx) = -dx[0]/2.0*(Ec(i,j,k,bxy)/dx[1]-slopesQ(i,j,k,X)); // +cxz/dx[2]
		  Ec(i,j,k,byy) = -dx[1]/2.0*(Ec(i,j,k,axy)/dx[0]-slopesQ(i,j,k,Y)); // ++cyz/dx[2]
		  Ec(i,j,k,czz) = 0.0;
		  Ec(i,j,k,a0) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - Ec(i,j,k,axx)/6.0;
		  Ec(i,j,k,ax) = (arrEM_X(i+1,j,k,EX_LOCAL)-arrEM_X(i,j,k,EX_LOCAL)) - Ec(i,j,k,axxx)/10.0;
		  Ec(i,j,k,b0) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - Ec(i,j,k,byy)/6.0;
		  Ec(i,j,k,by) = (arrEM_Y(i,j+1,k,EY_LOCAL)-arrEM_Y(i,j,k,EY_LOCAL)) - Ec(i,j,k,byyy)/10.0;
		  Ec(i,j,k,c0) = arr(i,j,k,EZ) - Ec(i,j,k,czz)/6.0; //czz*dx[2]*dx[2]/6.0
		  Ec(i,j,k,cz) = 0.0;		  
		  
		}
	    }
	}
    }

  // Compute edge components of the EM fields    
  for (MFIter mfi(fluxesEM, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      // use old array, S_source should also work
      //const auto& arr = S_source.array(mfi);
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi); 
      const auto& fluxArrEM = fluxesEM.array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{		 
		  
		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
		  // LD state
		  x = dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = -1, jOffset = -1;
		  Vector<Real> EM_LD = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  		  
		  // RD state
		  x = -dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_RD = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // LU state
		  x = dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_LU = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // RU state
		  x = -dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_RU = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // If the reconstruction is consistent, then EyRU=EyRD
		  Real EyR = (EM_RU[EY_LOCAL]+EM_RD[EY_LOCAL])/2.0;
		  Real EyL = (EM_LU[EY_LOCAL]+EM_LD[EY_LOCAL])/2.0;
		  Real ExU = (EM_RU[EX_LOCAL]+EM_LU[EX_LOCAL])/2.0;
		  Real ExD = (EM_RD[EX_LOCAL]+EM_LD[EX_LOCAL])/2.0;
		  fluxArrEM(i,j,k,BZ_LOCAL) = 0.25*(EM_RU[BZ_LOCAL]+EM_RD[BZ_LOCAL]+EM_LU[BZ_LOCAL]+EM_LD[BZ_LOCAL])
		    - 0.5/c*(EyR-EyL) + 0.5/c*(ExU-ExD);
		  
		  Real ByR = (EM_RU[BY_LOCAL]+EM_RD[BY_LOCAL])/2.0;
		  Real ByL = (EM_LU[BY_LOCAL]+EM_LD[BY_LOCAL])/2.0;
		  Real BxU = (EM_RU[BX_LOCAL]+EM_LU[BX_LOCAL])/2.0;
		  Real BxD = (EM_RD[BX_LOCAL]+EM_LD[BX_LOCAL])/2.0;

		  fluxArrEM(i,j,k,EZ_LOCAL) = 0.25*(EM_RU[EZ_LOCAL]+EM_RD[EZ_LOCAL]+EM_LU[EZ_LOCAL]+EM_LD[EZ_LOCAL])
		    + 0.5*c*(ByR-ByL) - 0.5*c*(BxU-BxD);
		  		  
		}
	    }
	}       
    }

  // Update face-centred EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
    { 

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      int d_EM = (d==0) ? 1 : 0;
    
      // Loop over all the patches at this level
      for (MFIter mfi(S_EM_dest[d_EM], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();

	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  // Indexable arrays for the data, and the directional flux
	  // Based on the corner-centred definition of the flux array, the
	  // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	  //const auto& arr = S_source.array(mfi);
	  const auto& arrEM = S_EM_dest[d_EM].array(mfi);
	  const auto& arrEMOld = S_EM_source[d_EM].array(mfi);	  
	  const auto& arr = S_source.array(mfi);
	  const auto& fluxArrEM = fluxesEM.array(mfi);
	  const auto& fluxArr = fluxes[d_EM].array(mfi);

	  const Dim3 hiDomain = ubound(geom.Domain());
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      // 2D code; z-component updated using cell-centred scheme
		      arrEM(i,j,k,BX_LOCAL+d_EM) = arrEMOld(i,j,k,BX_LOCAL+d_EM) + std::pow(-1,d)*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EZ_LOCAL) - fluxArrEM(i,j,k,EZ_LOCAL));
		      arrEM(i,j,k,EX_LOCAL+d_EM) = arrEMOld(i,j,k,EX_LOCAL+d_EM) - std::pow(-1,d)*c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) - fluxArrEM(i,j,k,BZ_LOCAL));

		      // 3D code
		      /*arrEM(i,j,k,BX_LOCAL+(1+d)%3) += (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(2+d)%3));
			arrEM(i,j,k,BX_LOCAL+(2+d)%3) -= (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(1+d)%3));
			arrEM(i,j,k,EX_LOCAL+(1+d)%3) -= c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(2+d)%3));
			arrEM(i,j,k,EX_LOCAL+(2+d)%3) += c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(1+d)%3));
		      */
		      
		      // source terms
		      Real currentFace = r_i*fluxArr(i,j,k,RHO_I) + r_e*fluxArr(i,j,k,RHO_E);
		      arrEM(i,j,k,EX_LOCAL+d_EM) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;
		      
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

  // Compute cell-centred EM fields from face-centred
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrEM_X = S_EM_dest[0].array(mfi);
      const auto& arrEM_Y = S_EM_dest[1].array(mfi); 

      //const auto& Bc = Bcoeff.array(mfi);
      //const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

  		  arr(i,j,k,BX) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0;// - Bc(i,j,k,axx)/6.0;
  		  arr(i,j,k,BY) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0;// - Bc(i,j,k,byy)/6.0;
  		  arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0;// - Ec(i,j,k,axx)/6.0;
  		  arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0;// - Ec(i,j,k,byy)/6.0;
		  
  		}
  	    }
  	}       
    }

  // Only 2D
  // Compute edge components for y-components at x-faces
  for (MFIter mfi(fluxes[0], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      const auto& fluxArrX = fluxes[0].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;
		  /*
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // L state
  		  x = dx[0]/2.0, y = 0.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_L = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		   
  		  // R state
  		  x = -dx[0]/2.0, y = 0.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_R = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  */
		  // two gaussian points
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // L state
		  iOffset = -1, jOffset = 0;
  		  x = dx[0]/2.0, y = dx[1]/(2.0*std::sqrt(3.0));		  
		  Vector<Real> EM_L1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = dx[0]/2.0, y = -dx[1]/(2.0*std::sqrt(3.0));		  
		  Vector<Real> EM_L2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);		  
		   
  		  // R state
		  iOffset = 0, jOffset = 0;
  		  x = -dx[0]/2.0, y = dx[1]/(2.0*std::sqrt(3.0));		 
		  Vector<Real> EM_R1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = -dx[0]/2.0, y = -dx[1]/(2.0*std::sqrt(3.0));		 
		  Vector<Real> EM_R2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  Vector<Real> EM_L, EM_R;
		  for (int n=0; n<6; n++)
		    {
		      EM_L.push_back(0.5*(EM_L1[n]+EM_L2[n]));
		      EM_R.push_back(0.5*(EM_R1[n]+EM_R2[n]));
		    }
		  
		  // Note that this is not the flux, but the star states
  		  fluxArrX(i,j,k,BY) = 0.5*(EM_R[BY_LOCAL]+EM_L[BY_LOCAL])
		    + 0.5/c*(EM_R[EZ_LOCAL]-EM_L[EZ_LOCAL]);
  		  fluxArrX(i,j,k,EY) = 0.5*(EM_R[EY_LOCAL]+EM_L[EY_LOCAL])
		    - 0.5*c*(EM_R[BZ_LOCAL]-EM_L[BZ_LOCAL]);
  		}
  	    }
  	}
    }
  
  // Only 2D
  // Compute edge components for x-components at y-faces
  for (MFIter mfi(fluxes[1], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 		  

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;
		  	
		  /*// For 2D code
		  z = 0.0;
		  kOffset = 0;

  		  // D state
  		  x = 0.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_D = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
  		  // U state
  		  x = 0.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_U = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  */
		  // two gaussian points
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // D state
  		  iOffset = 0, jOffset = -1;
		  x = dx[0]/(2.0*std::sqrt(3.0)), y = dx[1]/2.0;		  
		  Vector<Real> EM_D1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = -dx[0]/(2.0*std::sqrt(3.0)), y = dx[1]/2.0;		  
		  Vector<Real> EM_D2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
  		  // U state
  		  iOffset = 0, jOffset = 0;
		  x = dx[0]/(2.0*std::sqrt(3.0)), y = -dx[1]/2.0;		  
		  Vector<Real> EM_U1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = -dx[0]/(2.0*std::sqrt(3.0)), y = -dx[1]/2.0;		  
		  Vector<Real> EM_U2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  Vector<Real> EM_D, EM_U;
		  for (int n=0; n<6; n++)
		    {
		      EM_D.push_back(0.5*(EM_D1[n]+EM_D2[n]));
		      EM_U.push_back(0.5*(EM_U1[n]+EM_U2[n]));
		    }

		  // Note that this is not the flux, but the star states 	  
  		  fluxArrY(i,j,k,BX) = 0.5*(EM_U[BX_LOCAL]+EM_D[BX_LOCAL])
		    - 0.5/c*(EM_U[EZ_LOCAL]-EM_D[EZ_LOCAL]);
  		  fluxArrY(i,j,k,EX) = 0.5*(EM_U[EX_LOCAL]+EM_D[EX_LOCAL])
		    + 0.5*c*(EM_U[BZ_LOCAL]-EM_D[BZ_LOCAL]); 
  		}
  	    }
  	}
    }

  // Update cell-centred z-components of EM fields
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrOld = S_source.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);

      Array4<Real> slopes = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  // Update cell-centred z-components becuause it is 2D code
  		  arr(i,j,k,BZ) = arrOld(i,j,k,BZ) - (dt / dx[0]) * (fluxArrX(i+1,j,k,EY) - fluxArrX(i,j,k,EY)) + (dt / dx[1]) * (fluxArrY(i,j+1,k,EX) - fluxArrY(i,j,k,EX));
  		  arr(i,j,k,EZ) = arrOld(i,j,k,EZ) + c*c*(dt / dx[0]) * (fluxArrX(i+1,j,k,BY) - fluxArrX(i,j,k,BY)) - c*c*(dt / dx[1]) * (fluxArrY(i,j+1,k,BX) - fluxArrY(i,j,k,BX));		  

  		  // source terms
  		  arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
		  // four gaussian points
		  /*Real data_charge = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
		  Real slopesX = slopes(i,j,k,X);
		  Real slopesXX = slopes(i,j,k,XX);
		  Real slopesY = slopes(i,j,k,Y);
		  Real slopesYY = slopes(i,j,k,YY);
		  Real slopesXY = slopes(i,j,k,XY);
		  // 1st quadrature
		  Real x = -0.5, y = 0.0;
		  Real Lx = x, Lxx = x*x - 1.0/12.0;
		  Real Ly = y, Lyy = y*y - 1.0/12.0;
		  Real data_charge1 = data_charge + slopesX*Lx + slopesXX*Lxx +
		    slopesY*Ly + slopesYY*Lyy + slopesXY*Lx*Ly;
		  // 2nd quadrature
		  x = 0.5, y = 0.0;
		  Lx = x, Lxx = x*x - 1.0/12.0;
		  Ly = y, Lyy = y*y - 1.0/12.0;
		  Real data_charge2 = data_charge + slopesX*Lx + slopesXX*Lxx +
		    slopesY*Ly + slopesYY*Lyy + slopesXY*Lx*Ly;
		  // 3nd quadrature
		  x = 0.0, y = 0.5;
		  Lx = x, Lxx = x*x - 1.0/12.0;
		  Ly = y, Lyy = y*y - 1.0/12.0;
		  Real data_charge3 = data_charge + slopesX*Lx + slopesXX*Lxx +
		    slopesY*Ly + slopesYY*Lyy + slopesXY*Lx*Ly;
		  // 2nd quadrature
		  x = 0.0, y = -0.5;
		  Lx = x, Lxx = x*x - 1.0/12.0;
		  Ly = y, Lyy = y*y - 1.0/12.0;
		  Real data_charge4 = data_charge + slopesX*Lx + slopesXX*Lxx +
		    slopesY*Ly + slopesYY*Lyy + slopesXY*Lx*Ly;
		  arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*-1.25*(data_charge1+data_charge2+
		
							   data_charge3+data_charge4);		  
		  */
  		}
  	    }
  	}
    }
  // We need to compute boundary conditions again after each update
  S_dest.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_dest, geom, bc);  

  // Compute face-centred EM fields from cell-centred
  for (int d = 0; d < AMREX_SPACEDIM ; d++)
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      for (MFIter mfi(S_EM_dest[d], true); mfi.isValid(); ++mfi)
	{
	  const Box& box = mfi.tilebox();
      
	  const Dim3 lo = lbound(box);
	  const Dim3 hi = ubound(box);
      
	  // Indexable arrays for the data, and the directional flux
	  // Based on the corner-centred definition of the flux array, the
	  // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	  const auto& arr = S_dest.array(mfi);
	  const auto& arrEM = S_EM_dest[d].array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {		 

		      arrEM(i,j,k,BZ_LOCAL) = (arr(i,j,k,BZ)+arr(i-iOffset,j-jOffset,k,BZ))/2.0;
		      arrEM(i,j,k,EZ_LOCAL) = (arr(i,j,k,EZ)+arr(i-iOffset,j-jOffset,k,EZ))/2.0;
		    
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

}
/*
void CAMReXmp::setBC_FCEM(Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest)
{
  if (!geom.isAllPeriodic())
    //amrex::Print() << "Not all periodic" << std::endl;
    {
      for (int d = 0; d < amrex::SpaceDim ; d++)  
	{

	  const int iOffset = ( d == 0 ? 1 : 0);
	  const int jOffset = ( d == 1 ? 1 : 0);
	  const int kOffset = ( d == 2 ? 1 : 0);

	  if (!geom.isPeriodic(0))
	    {
	      //amrex::Print() << "x dir not periodic" << std::endl;
	      for (MFIter mfi(S_EM_dest[d], true); mfi.isValid(); ++mfi)
		{
		  const Box& bx = mfi.tilebox();

		  const Dim3 lo = lbound(bx);
		  const Dim3 hi = ubound(bx);

		  const Dim3 hiDomain = ubound(geom.Domain());

		  const auto& arrEM = S_EM_dest[d].array(mfi);	      
		  
		  for(int k = lo.z; k <= hi.z; k++)
		    {
		      for(int j = lo.y; j <= hi.y; j++)
			{
			  for(int i = lo.x-NUM_GROW; i <= hi.x+NUM_GROW; i++)
			    {			      
			      if (i==0)
				{
				  for (int n=0; n<6; n++)
				    {
				      if ((bc_EM[n].lo(0) == BCType::reflect_odd) &&
					  (bx.ixType().nodeCentered(0)))
					arrEM(i,j,k,n) = 0.0;
				    }
				}

			      else if (i<0)
				{
				  for (int n=0; n<6; n++)
				    {
				      if (bc_EM[n].lo(0) == BCType::foextrap)
					arrEM(i,j,k,n) = arrEM(0,j,k,n);
				      else if (bc_EM[n].lo(0) == BCType::reflect_even)
					arrEM(i,j,k,n) =  (bx.ixType().nodeCentered(0)) ? arrEM(-i,j,k,n) : arrEM(-i-1,j,k,n);
				      else if (bc_EM[n].lo(0) == BCType::reflect_odd)
					arrEM(i,j,k,n) =  (bx.ixType().nodeCentered(0)) ? -arrEM(-i,j,k,n) : -arrEM(-i-1,j,k,n);
				    }
				}
			      else if (i==hiDomain.x+iOffset)
				{
				  for (int n=0; n<6; n++)
				    {
				      if ((bc_EM[n].hi(0) == BCType::reflect_odd) &&
					  (bx.ixType().nodeCentered(0)))
					arrEM(i,j,k,n) = 0.0;
				    }
				  
				}
			      else if (i>hiDomain.x+iOffset)
				{
				  int ihi = hiDomain.x+1;
				  for (int n=0; n<6; n++)
				    {
				      if (bc_EM[n].hi(0) == BCType::foextrap)
					arrEM(i,j,k,n) = arrEM(ihi,j,k,n);
				      else if (bc_EM[n].hi(0) == BCType::reflect_even)
					arrEM(i,j,k,n) =  (bx.ixType().nodeCentered(0)) ? arrEM(2*ihi-i,j,k,n) : arrEM(2*ihi-i-1,j,k,n);
				      else if (bc_EM[n].hi(0) == BCType::reflect_odd)
					arrEM(i,j,k,n) =  (bx.ixType().nodeCentered(0)) ? -arrEM(2*ihi-i,j,k,n) : -arrEM(2*ihi-i-1,j,k,n);
				    }
				} 	
			    }
			}
		    }
		}
	    }
	  if (!geom.isPeriodic(1))
	    {
	      //amrex::Print() << "y dir not periodic" << std::endl;
	      for (MFIter mfi(S_EM_dest[d], true); mfi.isValid(); ++mfi)
		{
		  const Box& bx = mfi.tilebox();

		  const Dim3 lo = lbound(bx);
		  const Dim3 hi = ubound(bx);

		  const Dim3 hiDomain = ubound(geom.Domain());

		  const auto& arrEM = S_EM_dest[d].array(mfi);	      
		  
		  for(int k = lo.z; k <= hi.z; k++)
		    {
		      for(int j = lo.y-NUM_GROW; j <= hi.y+NUM_GROW; j++)
			{
			  for(int i = lo.x; i <= hi.x; i++)
			    {			    				
			      if (j==0)
				{
				  for (int n=0; n<6; n++)
				    {
				      if ((bc_EM[n].lo(1) == BCType::reflect_odd) &&
					  (bx.ixType().nodeCentered(1)))
					arrEM(i,j,k,n) = 0.0;
				    }
				}

			      else if (j<0)
				{
				  for (int n=0; n<6; n++)
				    {
				      if (bc_EM[n].lo(1) == BCType::foextrap)
					arrEM(i,j,k,n) = arrEM(i,0,k,n);
				      else if (bc_EM[n].lo(1) == BCType::reflect_even)
					arrEM(i,j,k,n) =  (bx.ixType().nodeCentered(1)) ? arrEM(i,-j,k,n) : arrEM(i,-j-1,k,n);
				      else if (bc_EM[n].lo(1) == BCType::reflect_odd)
					arrEM(i,j,k,n) =  (bx.ixType().nodeCentered(1)) ? -arrEM(i,-j,k,n) : -arrEM(i,-j-1,k,n);
				    }
				}
			      else if (j==hiDomain.y+jOffset)
				{
				  for (int n=0; n<6; n++)
				    {
				      if ((bc_EM[n].hi(1) == BCType::reflect_odd) &&
					  (bx.ixType().nodeCentered(1)))
					arrEM(i,j,k,n) = 0.0;
				    }
				  
				}
			      else if (j>hiDomain.y+jOffset)
				{
				  int jhi = hiDomain.y+1;
				  for (int n=0; n<6; n++)
				    {

				      //if (i==62)
				      //std::cout << "Before " << bc_EM[n].hi(1) << " " << BCType::foextrap << " " << arrEM(i,j,k,n) << " " << arrEM(i,2*jhi-j,k,n) << " " << arrEM(i,2*jhi-j-1,k,n) << std::endl;

				      if (bc_EM[n].hi(1) == BCType::foextrap)
					arrEM(i,j,k,n) = arrEM(i,jhi,k,n);
				      else if (bc_EM[n].hi(1) == BCType::reflect_even)
					arrEM(i,j,k,n) =  (bx.ixType().nodeCentered(1)) ? arrEM(i,2*jhi-j,k,n) : arrEM(i,2*jhi-j-1,k,n);
				      else if (bc_EM[n].hi(1) == BCType::reflect_odd)
					arrEM(i,j,k,n) =  (bx.ixType().nodeCentered(1)) ? -arrEM(i,2*jhi-j,k,n) : -arrEM(i,2*jhi-j-1,k,n);
				      //if (i==62)
				      //std::cout << "After " << bc_EM[n].hi(1) << " " << BCType::foextrap << " " << arrEM(i,2*jhi-j,k,n) << " " << arrEM(i,2*jhi-j,k,n) << " " << arrEM(i,2*jhi-j-1,k,n) << std::endl;
				    }
				}
			    }
			}
		    }
		}	      
	    }
	}          
    }  

}
*/
void CAMReXmp::Projection(const Real* dx, Real time)
{
  MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  Array<MultiFab,AMREX_SPACEDIM> S_EM_input;
  S_EM_input[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_input[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2) 
  S_EM_input[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);
#endif

  std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_lobc;// = {LinOpBCType::Periodic, LinOpBCType::Periodic};
  std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_hibc;// = {LinOpBCType::Periodic, LinOpBCType::Periodic};

  Geometry const* gg = AMReX::top()->getDefaultGeometry();

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      if (gg->isPeriodic(idim))
	{
	  mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
	}
      else
	{
	  if ((bc_EM[BX_LOCAL].lo(idim) == FOEXTRAP
	       && bc_EM[BY_LOCAL].lo(idim) == FOEXTRAP
	       && bc_EM[BZ_LOCAL].lo(idim) == FOEXTRAP) ||
	      (bc_EM[BX_LOCAL].lo(idim) == HOEXTRAP
	       && bc_EM[BY_LOCAL].lo(idim) == HOEXTRAP
	       && bc_EM[BZ_LOCAL].lo(idim) == HOEXTRAP) ||
	      (bc_EM[BX_LOCAL].lo(idim) == REFLECT_EVEN
	       && bc_EM[BY_LOCAL].lo(idim) == REFLECT_EVEN
	       && bc_EM[BZ_LOCAL].lo(idim) == REFLECT_EVEN))
	    mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Neumann;
	  else
	    mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Dirichlet;//reflect_odd;
	}
    }
  
  bool agglomeration = true;
  bool consolidation = true;
  int max_coarsening_level = 30;
  
  LPInfo info;
  info.setAgglomeration(agglomeration);
  info.setConsolidation(consolidation);
  info.setMaxCoarseningLevel(max_coarsening_level);
  //info.setAgglomeration(1);
  //info.setConsolidation(1);
  //info.setMetricTerm(false);
  
  // Implicit solve using MLPoisson class
  MLPoisson linop({geom}, {grids}, {dmap}, info);
  linop.setMaxOrder(max_order);
  
  // Set boundary conditions for MLPoisson  
  linop.setDomainBC(mlmg_lobc, mlmg_hibc);
  
  MultiFab phi(S_input, amrex::make_alias, DIVE, 1);
  phi.setVal(0.0);

  // Set boundary conditions for the current patch 
  linop.setLevelBC(0,&phi);
  
  MultiFab Rhs;
  Rhs.define(grids, dmap, 1, 0);
  for (MFIter mfi(Rhs, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& rhs = Rhs.array(mfi);
      const auto& arr = S_input.array(mfi);
        
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{	    	    
  		  // Initialise to zero
  		  rhs(i,j,k,0) = -1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
  		}
  	    }
  	}    
    }
  for (int d = 0; d < amrex::SpaceDim ; d++)
    {
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      for (MFIter mfi(Rhs, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  Array4<Real> rhs = Rhs.array(mfi);
	  Array4<Real> arr = S_input.array(mfi);
	  Array4<Real> arrEM = S_EM_input[d].array(mfi);	
	
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      //rhs(i,j,k,0) += (arrEM(i+iOffset,j+jOffset,k+kOffset,EX_LOCAL+d)-arrEM(i,j,k,EX_LOCAL+d))/dx[d];
		      rhs(i,j,k,0) += (arr(i+iOffset,j+jOffset,k,EX+d)-arr(i-iOffset,j-jOffset,k,EX+d))/(2.0*dx[d]);
		      if (geom.Coord()==1 && d==0)
			{
			  const Real x = geom.ProbLo()[0]+(double(i)+0.5)*dx[0];
			  //const Real y = geom.ProbLo()[1]+(double(j)+0.5)*dx[1];
			  rhs(i,j,k,0) += arr(i,j,k,EX+d)/x;
			}
		    }
		}
	    }      
	}
    }
  
  MLMG mlmg(linop);
  
  mlmg.setMaxFmgIter(max_fmg_iter);
  mlmg.setVerbose(verbose);
  
  const Real phi_abs = soln_tol*Rhs.norm0();

  // Solve to get S^(n+1)
  mlmg.solve({&phi}, {&Rhs}, soln_tol, phi_abs);
  //mlmg.solve({&phi}, {&Rhs}, 1e-10, 0.0);

  // We need to compute boundary conditions again after each update
  S_input.FillBoundary(geom.periodicity());
  // Fill non-periodic physical boundaries
  FillDomainBoundary(S_input, geom, bc);

  /*for (int d = 0; d < amrex::SpaceDim ; d++)   
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
		      arrEM(i,j,k,EX_LOCAL+d) -= (arr(i,j,k,DIVE)-arr(i-iOffset,j-jOffset,k-kOffset,DIVE))/dx[d];
		      //arrEM(i,j,k,EX_LOCAL+d) = 0.5*(arr(i-iOffset,j-jOffset,k-kOffset,EX+d)+arr(i,j,k,EX+d)) - (arr(i,j,k,DIVE)-arr(i-iOffset,j-jOffset,k-kOffset,DIVE))/dx[d];
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
  FillPatch(*this, S_EM_dest[0], NUM_GROW, dt, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_int, S_EM_dest[1], 0, 0, 6, 0);
  FillPatch(*this, S_EM_dest[1], NUM_GROW, dt, EM_Y_Type, 0, 6);
#endif
  */
  // Compute cell-centred EM fields from face-centred
  for (MFIter mfi(S_input, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = S_input.array(mfi);
      const auto& arrEM_X = S_EM_input[0].array(mfi);
      const auto& arrEM_Y = S_EM_input[1].array(mfi); 

      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 
  		  //arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0;
  		  //arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0;
		  arr(i,j,k,EX) -= (arr(i+1,j,k,DIVE)-arr(i-1,j,k,DIVE))/(2.0*dx[0]);
		  arr(i,j,k,EY) -= (arr(i,j+1,k,DIVE)-arr(i,j-1,k,DIVE))/(2.0*dx[1]); 
  		}
  	    }
  	}       
    }

  MultiFab& S_new = get_new_data(Phi_Type);
  //MultiFab::Copy(S_new, S_input, 0, 0, NUM_STATE, 0);
  MultiFab::Copy(S_new, S_input, EX, EX, amrex::SpaceDim, 0);

  /*// We need to compute boundary conditions again after each update
  S_dest.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_dest, geom, bc);  
  */
}

std::function<bool(int,int,int,int,int)> markerFunction(const BCRec& bc, const Box& nddom){

  return [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int /*n*/)
    -> bool
	 {
	   //return true;
	   // Dirichlet or reflective
	   if ((bc.lo(0)==EXT_DIR || bc.hi(0)==EXT_DIR)
	       || (bc.lo(0)==REFLECT_ODD && bc.hi(0)==REFLECT_ODD))
	     {
	       if (i==nddom.smallEnd(0) || i==nddom.bigEnd(0))
		 return false;
	     }
	   if ((bc.lo(1)==EXT_DIR || bc.hi(1)==EXT_DIR)
	       || (bc.lo(1)==REFLECT_ODD && bc.hi(1)==REFLECT_ODD))
	     {
	       if (j==nddom.smallEnd(1) || j==nddom.bigEnd(1))
		 return false;
	     }
#if (AMREX_SPACEDIM > 2)
	   if ((bc.lo(2)==EXT_DIR || bc.hi(2)==EXT_DIR)
	       || (bc.lo(2)==REFLECT_ODD && bc.hi(2)==REFLECT_ODD))
	     {
	       if (k==nddom.smallEnd(2) || k==nddom.bigEnd(2))
		 return false;
	     }
#endif
	   return true;
	 };
}
std::function<void(int,int,int,int,int,Array4<HYPRE_Int const> const*,HYPRE_Int&,HYPRE_Int*,HYPRE_Real*)> fillerFunction(const BCRec& bc, const Box& nddom, const amrex::GpuArray<double,AMREX_SPACEDIM> dxi, Real dt){
  
  Real coeff = 0.5*0.5*c*c*dt*dt;
  
  GpuArray<HYPRE_Real,AMREX_SPACEDIM> fac
    {AMREX_D_DECL(static_cast<HYPRE_Real>(dxi[0]*dxi[0]),
		  static_cast<HYPRE_Real>(dxi[1]*dxi[1]),
		  static_cast<HYPRE_Real>(dxi[2]*dxi[2]))};


  HYPRE_Real fac0 = HYPRE_Real(-2.)*(AMREX_D_TERM(fac[0],+fac[1],+fac[2]));

  Geometry const* gg = AMReX::top()->getDefaultGeometry();

  // For variable n at (i,j,k) in Box boxno (local index), fill its row in
  // the matrix.
  // [in ] gid : gid[n] is the id for variable n at (i,j,k)
  // [out] ncols: # of columns in this row.
  // [out] cols: column indices in this row.
  // [out] mat : matrix elemens in this row.
  return [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int n,
			       Array4<HYPRE_Int const> const* gid,
			       HYPRE_Int& ncols, HYPRE_Int* cols,
			       HYPRE_Real* mat)
	 {
	   ncols = 0;
	   ////////////////////////////////////////////////////
	   // x-direction
	   ////////////////////////////////////////////////////
	   // Periodic
	   if (gg->isPeriodic(0))
	     {
	       if (i >= nddom.smallEnd(0)) {
		 cols[ncols] = gid[n](i-1,j,k);
		 mat [ncols] = -coeff*fac[0];
		 ++ncols;
	       }
	       if (i <= nddom.bigEnd(0)) {
		 cols[ncols] = gid[n](i+1,j,k);
		 mat [ncols] = -coeff*fac[0];
		 ++ncols;
	       }
	     }
	   // Dirichlet (and reflective)
	   else if ((bc.lo(0)==EXT_DIR && bc.hi(0)==EXT_DIR)//)
		    || (bc.lo(0)==REFLECT_ODD && bc.hi(0)==REFLECT_ODD))
	     {
	       if (i > nddom.smallEnd(0)+1) {
		 cols[ncols] = gid[n](i-1,j,k);
		 mat [ncols] = -coeff*fac[0];
		 ++ncols;
	       }
	       if (i < nddom.bigEnd(0)-1) {
		 cols[ncols] = gid[n](i+1,j,k);
		 mat [ncols] = -coeff*fac[0];
		 ++ncols;
	       }
	     }
	   // Neumann
	   else if ((bc.lo(0)==FOEXTRAP && bc.hi(0)==FOEXTRAP)
		    || (bc.lo(0)==HOEXTRAP && bc.hi(0)==HOEXTRAP)
		    || (bc.lo(0)==REFLECT_EVEN && bc.hi(0)==REFLECT_EVEN))
	     {
	       if (i > nddom.smallEnd(0) && i !=nddom.bigEnd(0)) {
		 cols[ncols] = gid[n](i-1,j,k);
		 mat [ncols] = -coeff*fac[0];
		 ++ncols;
	       }
	       if (i < nddom.bigEnd(0) && i !=nddom.smallEnd(0)) {
		 cols[ncols] = gid[n](i+1,j,k);
		 mat [ncols] = -coeff*fac[0];
		 ++ncols;
	       }
	       if (i == nddom.smallEnd(0)) {
		 cols[ncols] = gid[n](i+1,j,k);
		 mat [ncols] = -2*coeff*fac[0];
		 ++ncols;
	       }
	       if (i == nddom.bigEnd(0)) {
		 cols[ncols] = gid[n](i-1,j,k);
		 mat [ncols] = -2*coeff*fac[0];
		 ++ncols;
	       }
	     }
	   // Reflective
	   // Note that this is only when it is not nodal
	   // when it is nodal, it is 0 at the boundaries
	   // else if (bc.lo(0)==REFLECT_ODD && bc.hi(0)==REFLECT_ODD)
	   //   {
	   //     /*if (i > nddom.smallEnd(0)) {
	   // 	 cols[ncols] = gid[n](i-1,j,k);
	   // 	 mat [ncols] = -coeff*fac[0];
	   // 	 ++ncols;
	   //     }
	   //     if (i < nddom.bigEnd(0)) {
	   // 	 cols[ncols] = gid[n](i+1,j,k);
	   // 	 mat [ncols] = -coeff*fac[0];
	   // 	 ++ncols;
	   //     }
	   //     if (i == nddom.bigEnd(0) || i == nddom.smallEnd(0)) {
	   // 	 cols[ncols] = gid[n](i,j,k);
	   // 	 // opposite sign for reflective boundaries
	   // 	 mat [ncols] = coeff*fac[0];
	   // 	 ++ncols;
	   // 	 }*/
	   //     if (i > nddom.smallEnd(0) && i !=nddom.bigEnd(0)) {
	   // 	 cols[ncols] = gid[n](i-1,j,k);
	   // 	 mat [ncols] = -coeff*fac[0];
	   // 	 ++ncols;
	   //     }
	   //     if (i < nddom.bigEnd(0) && i !=nddom.smallEnd(0)) {
	   // 	 cols[ncols] = gid[n](i+1,j,k);
	   // 	 mat [ncols] = -coeff*fac[0];
	   // 	 ++ncols;
	   //     }
	   //     if (i == nddom.smallEnd(0)) {
	   // 	 cols[ncols] = gid[n](i+1,j,k);
	   // 	 mat [ncols] = -(4.0/3.0)*coeff*fac[0];
	   // 	 ++ncols;
	   // 	 cols[ncols] = gid[n](i,j,k);
	   // 	 mat [ncols] = 2.0*coeff*fac[0];
	   // 	 ++ncols;
	   //     }
	   //     if (i == nddom.bigEnd(0)) {
	   // 	 cols[ncols] = gid[n](i-1,j,k);
	   // 	 mat [ncols] = -(4.0/3.0)*coeff*fac[0];
	   // 	 ++ncols;
	   // 	 cols[ncols] = gid[n](i,j,k);
	   // 	 mat [ncols] = 2.0*coeff*fac[0];
	   // 	 ++ncols;
	   //     }	       
	   //   }
	   else
	     amrex::Abort("Unsupported BC in x-direction");
	   
	   ////////////////////////////////////////////////////
	   // y-direction
	   ////////////////////////////////////////////////////
	   // Periodic
	   if (gg->isPeriodic(1))
	     {
	       if (j >= nddom.smallEnd(1)) {
		 cols[ncols] = gid[n](i,j-1,k);
		 mat [ncols] = -coeff*fac[1];
		 ++ncols;
	       }
	       if (j <= nddom.bigEnd(1)) {
		 cols[ncols] = gid[n](i,j+1,k);
		 mat [ncols] = -coeff*fac[1];
		 ++ncols;
	       }
	     }
	   // Dirichlet
	   else if ((bc.lo(1)==EXT_DIR && bc.hi(1)==EXT_DIR)//)
		    || (bc.lo(1)==REFLECT_ODD && bc.hi(1)==REFLECT_ODD))
	     {
	       if (j > nddom.smallEnd(1)+1) {
		 cols[ncols] = gid[n](i,j-1,k);
		 mat [ncols] = -coeff*fac[1];
		 ++ncols;
	       }
	       if (j < nddom.bigEnd(1)-1) {
		 cols[ncols] = gid[n](i,j+1,k);
		 mat [ncols] = -coeff*fac[1];
		 ++ncols;
	       }
	     }
	   // Neumann
	   else if ((bc.lo(1)==FOEXTRAP && bc.hi(1)==FOEXTRAP)
		    || (bc.lo(1)==HOEXTRAP && bc.hi(1)==HOEXTRAP)
		    || (bc.lo(1)==REFLECT_EVEN && bc.hi(1)==REFLECT_EVEN))
	     {
	       if (j > nddom.smallEnd(1) && j !=nddom.bigEnd(1)) {
		 cols[ncols] = gid[n](i,j-1,k);
		 mat [ncols] = -coeff*fac[1];
		 ++ncols;
	       }
	       if (j < nddom.bigEnd(1) && j !=nddom.smallEnd(1)) {
		 cols[ncols] = gid[n](i,j+1,k);
		 mat [ncols] = -coeff*fac[1];
		 ++ncols;
	       }
	       if (j == nddom.smallEnd(1)) {
		 cols[ncols] = gid[n](i,j+1,k);
		 mat [ncols] = -2*coeff*fac[1];
		 ++ncols;
	       }
	       if (j == nddom.bigEnd(1)) {
		 cols[ncols] = gid[n](i,j-1,k);
		 mat [ncols] = -2*coeff*fac[1];
		 ++ncols;
	       }
	     }
	   // Reflective
	   // Note that this is only when it is not nodal
	   // when it is nodal, it is 0 at the boundaries
	   // else if (bc.lo(1)==REFLECT_ODD && bc.hi(1)==REFLECT_ODD)
	   //   {
	   //     /*if (j > nddom.smallEnd(1)) {
	   // 	 cols[ncols] = gid[n](i,j-1,k);
	   // 	 mat [ncols] = -coeff*fac[1];
	   // 	 ++ncols;
	   //     }
	   //     if (j < nddom.bigEnd(1)) {
	   // 	 cols[ncols] = gid[n](i,j+1,k);
	   // 	 mat [ncols] = -coeff*fac[1];
	   // 	 ++ncols;
	   //     }
	   //     if (j == nddom.bigEnd(1) || j == nddom.smallEnd(1)) {
	   // 	 cols[ncols] = gid[n](i,j,k);
	   // 	 // opposite sign for reflective boundaries
	   // 	 mat [ncols] = coeff*fac[1];
	   // 	 ++ncols;
	   //     }
	   //     */
	   //     if (j > nddom.smallEnd(1) && j !=nddom.bigEnd(1)) {
	   // 	 cols[ncols] = gid[n](i,j-1,k);
	   // 	 mat [ncols] = -coeff*fac[1];
	   // 	 ++ncols;
	   //     }
	   //     if (j < nddom.bigEnd(1) && j !=nddom.smallEnd(1)) {
	   // 	 cols[ncols] = gid[n](i,j+1,k);
	   // 	 mat [ncols] = -coeff*fac[1];
	   // 	 ++ncols;
	   //     }
	   //     if (j == nddom.smallEnd(1)) {
	   // 	 cols[ncols] = gid[n](i,j+1,k);
	   // 	 mat [ncols] = -(4.0/3.0)*coeff*fac[1];
	   // 	 ++ncols;
	   // 	 cols[ncols] = gid[n](i,j,k);
	   // 	 mat [ncols] = 2.0*coeff*fac[1];
	   // 	 ++ncols;
	   //     }
	   //     if (j == nddom.bigEnd(1)) {
	   // 	 cols[ncols] = gid[n](i,j-1,k);
	   // 	 mat [ncols] = -(4.0/3.0)*coeff*fac[1];
	   // 	 ++ncols;
	   // 	 cols[ncols] = gid[n](i,j,k);
	   // 	 mat [ncols] = 2.0*coeff*fac[1];
	   // 	 ++ncols;
	   //     }	       
	   //   }
	   else
	     amrex::Abort("Unsupported BC in y-direction");
#if (AMREX_SPACEDIM > 2)
	   ////////////////////////////////////////////////////
	   // z-direction
	   ////////////////////////////////////////////////////
	   // Periodic
	   if (gg->isPeriodic(2))
	     {
	       if (k >= nddom.smallEnd(2)) {
		 cols[ncols] = gid[n](i,j,k-1);
		 mat [ncols] = -coeff*fac[2];
		 ++ncols;
	       }
	       if (k <= nddom.bigEnd(2)) {
		 cols[ncols] = gid[n](i,j,k+1);
		 mat [ncols] = -coeff*fac[2];
		 ++ncols;
	       }
	     }
	   // Dirichlet
	   else if ((bc.lo(2)==EXT_DIR && bc.hi(2)==EXT_DIR)
		    || (bc.lo(2)==REFLECT_ODD && bc.hi(2)==REFLECT_ODD))
	     {
	       if (k > nddom.smallEnd(2)+1) {
		 cols[ncols] = gid[n](i,j,k-1);
		 mat [ncols] = -coeff*fac[2];
		 ++ncols;
	       }
	       if (k < nddom.bigEnd(2)-1) {
		 cols[ncols] = gid[n](i,j,k+1);
		 mat [ncols] = -coeff*fac[2];
		 ++ncols;
	       }
	     }
	   // Neumann
	   else if ((bc.lo(2)==FOEXTRAP && bc.hi(2)==FOEXTRAP)
		    || (bc.lo(2)==HOEXTRAP && bc.hi(2)==HOEXTRAP)
		    || (bc.lo(2)==REFLECT_EVEN && bc.hi(2)==REFLECT_EVEN))
	     {
	       if (k > nddom.smallEnd(2) && k !=nddom.bigEnd(2)) {
		 cols[ncols] = gid[n](i,j,k-1);
		 mat [ncols] = -coeff*fac[2];
		 ++ncols;
	       }
	       if (k < nddom.bigEnd(2) && k !=nddom.smallEnd(2)) {
		 cols[ncols] = gid[n](i,j,k+1);
		 mat [ncols] = -coeff*fac[2];
		 ++ncols;
	       }
	       if (k == nddom.smallEnd(2)) {
		 cols[ncols] = gid[n](i,j,k+1);
		 mat [ncols] = -2*coeff*fac[2];
		 ++ncols;
	       }
	       if (k == nddom.bigEnd(2)) {
		 cols[ncols] = gid[n](i,j,k-1);
		 mat [ncols] = -2*coeff*fac[2];
		 ++ncols;
	       }
	     }
	   // Reflective
	   // Note that this is only when it is not nodal
	   // when it is nodal, it is 0 at the boundaries
	   /*else if (bc.lo(2)==REFLECT_ODD && bc.hi(2)==REFLECT_ODD)
	     {
	       if (k > nddom.smallEnd(2)) {
		 cols[ncols] = gid[n](i,j,k-1);
		 mat [ncols] = -coeff*fac[2];
		 ++ncols;
	       }
	       if (k < nddom.bigEnd(2)) {
		 cols[ncols] = gid[n](i,j,k+1);
		 mat [ncols] = -coeff*fac[2];
		 ++ncols;
	       }
	       if (k == nddom.bigEnd(2) || k == nddom.smallEnd(2)) {
		 cols[ncols] = gid[n](i,j,k);
		 // opposite sign for reflective boundaries
		 mat [ncols] = coeff*fac[2];
		 ++ncols;
	       }
	       }*/
	   else
	     amrex::Abort("Unsupported BC in z-direction");
#endif
	   cols[ncols] = gid[n](i,j,k);
	   mat [ncols] = 1. - coeff*fac0;
	   ++ncols;
	 };
}
void CAMReXmp::implicitYeeMaxwellSolver(Array<MultiFab,AMREX_SPACEDIM>& S_EM_dest, Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& S_dest, MultiFab& S_source, const Real* dx, Real dt, Real time) 
{

  Real coeff = 0.5*0.5*c*c*dt*dt;

  IndexType xface(IntVect{AMREX_D_DECL(1,0,0)});
  IndexType yface(IntVect{AMREX_D_DECL(0,1,0)});
  IndexType edge(IntVect{AMREX_D_DECL(1,1,0)});

  auto const& nddom = amrex::surroundingNodes(geom.Domain());
  auto const& nddomX = amrex::surroundingNodes(geom.Domain(),0);
  auto const& nddomY = amrex::surroundingNodes(geom.Domain(),1);
   const auto dxi = geom.InvCellSizeArray();

  GpuArray<HYPRE_Real,AMREX_SPACEDIM> fac
    {AMREX_D_DECL(static_cast<HYPRE_Real>(dxi[0]*dxi[0]),
  		  static_cast<HYPRE_Real>(dxi[1]*dxi[1]),
  		  static_cast<HYPRE_Real>(dxi[2]*dxi[2]))};


  HYPRE_Real fac0 = HYPRE_Real(-2.)*(AMREX_D_TERM(fac[0],+fac[1],+fac[2]));

   if (bc_EM[BX_LOCAL].lo(1)==REFLECT_ODD || bc_EM[BY_LOCAL].lo(0)==REFLECT_ODD)
     amrex::Abort("Currently does not support fillerFunction reflec_odd for Bx in y-direction or By in x-direction");

   //functor that returns whether the variable n at (i,j,k) in Box boxno (local index) is valid
   auto marker = markerFunction(bc_EM[BZ_LOCAL],nddom);
   auto markerX = markerFunction(bc_EM[BY_LOCAL],nddomX);
   auto markerY = markerFunction(bc_EM[BX_LOCAL],nddomY);
  // Is variable n at (i,j,k) in Box boxno (local index) valid?
  // (i.e., not exactly on Dirichlet boundary)
  //  auto marker = [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int /*n*/)
  //   -> bool
  // 		{
  // 		  //return nddom.strictly_contains(i,j,k);
  // 		  return true;
  // 		};
  // auto markerX = [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int /*n*/)
  //   -> bool
  // 		{
  // 		  //return nddomX.strictly_contains(i,j,k);
  // 		  //return true;
  // 		  if (j==nddomX.smallEnd(1) || j==nddomX.bigEnd(1))
  // 		    return false;
  // 		  else
  // 		    return true;
  // 		};
  // auto markerY = [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int /*n*/)
  //   -> bool
  // 		{
  // 		  //return nddomY.strictly_contains(i,j,k);
  // 		  return true;
  // 		};  

  // functor that fills the row in the matrix A for variable n at (i,j,k) in Box boxno (local index)
  // using the CSR format
  auto filler = fillerFunction(bc_EM[BZ_LOCAL],nddom,dxi,dt);
  auto fillerX = fillerFunction(bc_EM[BY_LOCAL],nddomX,dxi,dt);
  auto fillerY = fillerFunction(bc_EM[BX_LOCAL],nddomY,dxi,dt);

//   auto filler = [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int n,
// 				      Array4<HYPRE_Int const> const* gid,
// 				      HYPRE_Int& ncols, HYPRE_Int* cols,
// 				      HYPRE_Real* mat)
// 		{
// 		  ncols = 0;
// 		  //if (i > nddom.smallEnd(0)+1) {
// 		  if (i >= nddom.smallEnd(0)) {
// 		    cols[ncols] = gid[n](i-1,j,k);
// 		    mat [ncols] = -coeff*fac[0];
// 		    ++ncols;
// 		  }
// 		  //if (i < nddom.bigEnd(0)-1) {
// 		  if (i <= nddom.bigEnd(0)) {
// 		    cols[ncols] = gid[n](i+1,j,k);
// 		    mat [ncols] = -coeff*fac[0];
// 		    ++ncols;
// 		  }
// 		  // Dirichlet
// 		  //if (j > nddom.smallEnd(1)+1) {
// 		  // Periodic
// 		  //if (j >= nddom.smallEnd(1)) {
// 		  // Neumann
// 		  if (j > nddom.smallEnd(1) && j!=nddom.bigEnd(1)) {		  
// 		    cols[ncols] = gid[n](i,j-1,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;
// 		  }
// 		  //if (j < nddom.bigEnd(1)-1) {
// 		  //if (j <= nddom.bigEnd(1)) {
// 		  if (j < nddom.bigEnd(1) && j!=nddom.smallEnd(1)) {
// 		    cols[ncols] = gid[n](i,j+1,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;
// 		  }
// 		  // Neumann
// 		  if (j==nddom.smallEnd(1)) {
// 		    cols[ncols] = gid[n](i,j+1,k);
// 		    mat [ncols] = -2*coeff*fac[1];
// 		    ++ncols;
// 		    /*cols[ncols] = gid[n](i,j+1,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;
// 		    cols[ncols] = gid[n](i,j,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;*/
// 		  }
// 		  if (j==nddom.bigEnd(1)) {
// 		    cols[ncols] = gid[n](i,j-1,k);
// 		    mat [ncols] = -2*coeff*fac[1];
// 		    ++ncols;
// 		    /*cols[ncols] = gid[n](i,j-1,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;
// 		    cols[ncols] = gid[n](i,j,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;*/
// 		  }
// #if (AMREX_SPACEDIM > 2)
// 		  if (k > nddom.smallEnd(2)+1) {
// 		    cols[ncols] = gid[n](i,j,k-1);
// 		    mat [ncols] = -coeff*fac[2];
// 		    ++ncols;
// 		  }
// 		  if (k < nddom.bigEnd(2)-1) {
// 		    cols[ncols] = gid[n](i,j,k+1);
// 		    mat [ncols] = -coeff*fac[2];
// 		    ++ncols;
// 		  }
// #endif
// 		  cols[ncols] = gid[n](i,j,k);
// 		  mat [ncols] = 1. - coeff*fac0;
// 		  ++ncols;
// 		};
//   auto fillerX = [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int n,
// 				      Array4<HYPRE_Int const> const* gid,
// 				      HYPRE_Int& ncols, HYPRE_Int* cols,
// 				      HYPRE_Real* mat)
// 		{
// 		  ncols = 0;
// 		  //if (i > nddomX.smallEnd(0)+1) {
// 		  if (i >= nddomX.smallEnd(0)) {
// 		    cols[ncols] = gid[n](i-1,j,k);
// 		    mat [ncols] = -coeff*fac[0];
// 		    ++ncols;
// 		  }
// 		  //if (i < nddomX.bigEnd(0)-1) {
// 		  if (i <= nddomX.bigEnd(0)) {
// 		    cols[ncols] = gid[n](i+1,j,k);
// 		    mat [ncols] = -coeff*fac[0];
// 		    ++ncols;
// 		  }
// 		  // Dirichlet
// 		  if (j > nddomX.smallEnd(1)+1) {
// 		  // Reflective
// 		  //if (j >= nddomX.smallEnd(1)+1) {
// 		  // Periodic
// 		  //if (j >= nddomX.smallEnd(1)) {
// 		    cols[ncols] = gid[n](i,j-1,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;
// 		  }
// 		  if (j < nddomX.bigEnd(1)-1) {
// 		    //if (j <= nddomX.bigEnd(1)-1) {
// 		  //if (j <= nddomX.bigEnd(1)) {
// 		    cols[ncols] = gid[n](i,j+1,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;
// 		  }
// 		  // Reflective
// 		  /*if (j == nddomX.bigEnd(1) || j == nddomX.smallEnd(1)) {
// 		    cols[ncols] = gid[n](i,j,k);
// 		    // opposite sign for reflective boundaries
// 		    mat [ncols] = coeff*fac[1];
// 		    ++ncols;
// 		    }*/
// #if (AMREX_SPACEDIM > 2)
// 		  if (k > nddomX.smallEnd(2)+1) {
// 		    cols[ncols] = gid[n](i,j,k-1);
// 		    mat [ncols] = -coeff*fac[2];
// 		    ++ncols;
// 		  }
// 		  if (k < nddomX.bigEnd(2)-1) {
// 		    cols[ncols] = gid[n](i,j,k+1);
// 		    mat [ncols] = -coeff*fac[2];
// 		    ++ncols;
// 		  }
// #endif
// 		  cols[ncols] = gid[n](i,j,k);
// 		  mat [ncols] = 1. - coeff*fac0;
// 		  ++ncols;
// 		};
//   auto fillerY = [=] AMREX_GPU_DEVICE (int /*boxno*/, int i, int j, int k, int n,
// 				      Array4<HYPRE_Int const> const* gid,
// 				      HYPRE_Int& ncols, HYPRE_Int* cols,
// 				      HYPRE_Real* mat)
// 		{
// 		  ncols = 0;
// 		  //if (i > nddomY.smallEnd(0)+1) {
// 		  if (i >= nddomY.smallEnd(0)) {
// 		    cols[ncols] = gid[n](i-1,j,k);
// 		    mat [ncols] = -coeff*fac[0];
// 		    ++ncols;
// 		  }
// 		  //if (i < nddomY.bigEnd(0)-1) {
// 		  if (i <= nddomY.bigEnd(0)) {
// 		    cols[ncols] = gid[n](i+1,j,k);
// 		    mat [ncols] = -coeff*fac[0];
// 		    ++ncols;
// 		  }
// 		  //if (j > nddomY.smallEnd(1)+1) {
// 		  //if (j >= nddomY.smallEnd(1)) {
// 		  if (j > nddomY.smallEnd(1) && j!=nddomY.bigEnd(1)) {
// 		    cols[ncols] = gid[n](i,j-1,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;
// 		  }
// 		  //if (j < nddomY.bigEnd(1)-1) {
// 		  //if (j <= nddomY.bigEnd(1)) {
// 		  if (j < nddomY.bigEnd(1) && j!=nddomY.smallEnd(1)) {
// 		    cols[ncols] = gid[n](i,j+1,k);
// 		    mat [ncols] = -coeff*fac[1];
// 		    ++ncols;
// 		  }
// 		  if (j==nddomY.smallEnd(1)) {
// 		    cols[ncols] = gid[n](i,j+1,k);
// 		    mat [ncols] = -2*coeff*fac[1];
// 		    ++ncols;		    
// 		  }
// 		  if (j==nddomY.bigEnd(1)) {
// 		    cols[ncols] = gid[n](i,j-1,k);
// 		    mat [ncols] = -2*coeff*fac[1];
// 		    ++ncols;
// 		  }
// #if (AMREX_SPACEDIM > 2)
// 		  if (k > nddomY.smallEnd(2)+1) {
// 		    cols[ncols] = gid[n](i,j,k-1);
// 		    mat [ncols] = -coeff*fac[2];
// 		    ++ncols;
// 		  }
// 		  if (k < nddomY.bigEnd(2)-1) {
// 		    cols[ncols] = gid[n](i,j,k+1);
// 		    mat [ncols] = -coeff*fac[2];
// 		    ++ncols;
// 		  }
// #endif
// 		  cols[ncols] = gid[n](i,j,k);
// 		  mat [ncols] = 1. - coeff*fac0;
// 		  ++ncols;
// 		};
  
  constexpr int max_stencil_size = 2*AMREX_SPACEDIM+1;
  HypreSolver<max_stencil_size> hypre_solverX
    ({xface}, IntVect(1), geom, grids, dmap,
     markerX, fillerX, 2);
  HypreSolver<max_stencil_size> hypre_solverY
    ({yface}, IntVect(1), geom, grids, dmap,
     markerY, fillerY, 2);
  HypreSolver<max_stencil_size> hypre_solverXY
    //({IndexType::TheNodeType()}, IntVect(1), geom, grids, dmap,
    ({edge}, IntVect(1), geom, grids, dmap,
     marker, filler, 2);
      
  MultiFab& S_EM_XY_new = get_new_data(EM_XY_Type);
  MultiFab S_EM_sourceEdge, S_EM_destEdge;
  S_EM_sourceEdge.define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 6, NUM_GROW);
  S_EM_destEdge.define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_sourceEdge, NUM_GROW, time, EM_XY_Type, 0, 6);
  FillPatch(*this, S_EM_destEdge, NUM_GROW, time, EM_XY_Type, 0, 6);

  MultiFab::Copy(S_dest, S_source, BX, BX, 6, NUM_GROW);
  MultiFab::Copy(S_EM_dest[0], S_EM_source[0], 0, 0, 6, NUM_GROW);
  MultiFab::Copy(S_EM_dest[1], S_EM_source[1], 0, 0, 6, NUM_GROW);

  /*LPInfo info;
  info.setAgglomeration(1);
  info.setConsolidation(1);
  info.setMetricTerm(false);
  
  // Implicit solve using MLABecLaplacian class
  //MLABecLaplacian mlabecX({geom}, {convert(grids,IntVect{AMREX_D_DECL(0,1,0)})}, {dmap}, info);
  MLABecLaplacian mlabecX({geom}, {grids}, {dmap}, info);
  mlabecX.setMaxOrder(max_order);
  //MLABecLaplacian mlabecY({geom}, {convert(grids,IntVect{AMREX_D_DECL(1,0,0)})}, {dmap}, info);
  MLABecLaplacian mlabecY({geom}, {grids}, {dmap}, info);
  mlabecY.setMaxOrder(max_order);
  //MLABecLaplacian mlabecZ({geom}, {convert(grids,IntVect{AMREX_D_DECL(1,1,0)})}, {dmap}, info);
  MLABecLaplacian mlabecZ({geom}, {grids}, {dmap}, info);
  mlabecZ.setMaxOrder(max_order);
    
  // Set boundary conditions for MLABecLaplacian  
  mlabecX.setDomainBC(mlmg_lobc_X, mlmg_hibc_X);
  mlabecY.setDomainBC(mlmg_lobc_Y, mlmg_hibc_Y);
  mlabecZ.setDomainBC(mlmg_lobc_Z, mlmg_hibc_Z);
  */
  MultiFab S_X(S_EM_dest[1], amrex::make_alias, BX_LOCAL, 1);
  MultiFab S_Y(S_EM_dest[0], amrex::make_alias, BY_LOCAL, 1);
  MultiFab S_Z(S_EM_destEdge, amrex::make_alias, BZ_LOCAL, 1);  
  
  //MultiFab S_X(S_dest, amrex::make_alias, BX, 1);
  //MultiFab S_Y(S_dest, amrex::make_alias, BY, 1);
  //MultiFab S_Z(S_dest, amrex::make_alias, BZ, 1);  

  /*// Set boundary conditions for the current patch 
  mlabecX.setLevelBC(0,&S_X);
  mlabecY.setLevelBC(0,&S_Y);
  mlabecZ.setLevelBC(0,&S_Z);  
  */
  // Set the RHS to be multiplied appropriately (by rho * c_v, i.e.\ acoef)
  // for(MFIter RHSmfi(Rhs,true); RHSmfi.isValid(); ++RHSmfi)
  // {
  //   const Box& box = RHSmfi.tilebox();
  //   Rhs[RHSmfi].mult(acoef[RHSmfi],box,0,0,1);
  // }
  std::array<MultiFab,3> Rhs;
  Rhs[0].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 1, NUM_GROW);
  Rhs[1].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 1, NUM_GROW);
  Rhs[2].define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 1, NUM_GROW);
  //Rhs[0].define(grids, dmap, 1, 0);
  //Rhs[1].define(grids, dmap, 1, 0);
  //Rhs[2].define(grids, dmap, 1, 0);

  for(MFIter mfi(Rhs[0], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs[0].array(mfi);
      const auto& arr = S_source.array(mfi);
      const auto& arrEM = S_EM_source[1].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  /////////////////////////////////////////////////////////////////////////
		  // 2D
		  Real dyEz = (arr(i,j,k,EZ)-arr(i,j-1,k,EZ))/dx[1];
		  Real dxdxBx = computeSecondDerivative(arrEM(i-1,j,k,BX_LOCAL), arrEM(i,j,k,BX_LOCAL), arrEM(i+1,j,k,BX_LOCAL), dx[0]);
		  Real dydyBx = computeSecondDerivative(arrEM(i,j-1,k,BX_LOCAL), arrEM(i,j,k,BX_LOCAL), arrEM(i,j+1,k,BX_LOCAL), dx[1]);
		  rhs(i,j,k) = arrEM(i,j,k,BX_LOCAL) + 0.25*c*c*dt*dt*(dxdxBx+dydyBx) - dt*dyEz;
		  
		  /*Real dyEz = computeDerivative(arr(i,j-1,k,EZ), arr(i,j+1,k,EZ), dx[1]);
		  Real dxdxBx = computeSecondDerivative(arr(i-1,j,k,BX), arr(i,j,k,BX), arr(i+1,j,k,BX), dx[0]);
		  Real dydyBx = computeSecondDerivative(arr(i,j-1,k,BX), arr(i,j,k,BX), arr(i,j+1,k,BX), dx[1]);
		  rhs(i,j,k) = arr(i,j,k,BX) + 0.25*c*c*dt*dt*(dxdxBx+dydyBx) - dt*dyEz;
		  */
		}
	    }
	}
    }
  for(MFIter mfi(Rhs[1], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
	  
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs[1].array(mfi);
      const auto& arr = S_source.array(mfi);
      const auto& arrEM = S_EM_source[0].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dxEz = (arr(i,j,k,EZ)-arr(i-1,j,k,EZ))/dx[0];
		  Real dxdxBy = computeSecondDerivative(arrEM(i-1,j,k,BY_LOCAL), arrEM(i,j,k,BY_LOCAL), arrEM(i+1,j,k,BY_LOCAL), dx[0]);
		  Real dydyBy = computeSecondDerivative(arrEM(i,j-1,k,BY_LOCAL), arrEM(i,j,k,BY_LOCAL), arrEM(i,j+1,k,BY_LOCAL), dx[1]);
		  rhs(i,j,k) = arrEM(i,j,k,BY_LOCAL) + 0.25*c*c*dt*dt*(dxdxBy+dydyBy) + dt*dxEz;

		  /*Real dxEz = computeDerivative(arr(i-1,j,k,EZ), arr(i+1,j,k,EZ), dx[0]);
		  Real dxdxBy = computeSecondDerivative(arr(i-1,j,k,BY), arr(i,j,k,BY), arr(i+1,j,k,BY), dx[0]);
		  Real dydyBy = computeSecondDerivative(arr(i,j-1,k,BY), arr(i,j,k,BY), arr(i,j+1,k,BY), dx[1]);		 
		  rhs(i,j,k) = arr(i,j,k,BY) + 0.25*c*c*dt*dt*(dxdxBy+dydyBy) + dt*dxEz;
		  */
		  /*if (rhs(i,j,k) != rhs(i,j,k))
		    //std::cout << "NaN value at " << i << " " << j << std::endl;
		    amrex::Abort("1");
		  */
		}
	    }
	}
    }
  for(MFIter mfi(Rhs[2], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
	  
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs[2].array(mfi);
      const auto& arr = S_source.array(mfi);
      const auto& arrEMX = S_EM_source[0].array(mfi);
      const auto& arrEMY = S_EM_source[1].array(mfi);
      const auto& arrEMXY = S_EM_sourceEdge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dyEx = (arrEMX(i,j,k,EX_LOCAL)-arrEMX(i,j-1,k,EX_LOCAL))/dx[1];
		  Real dxEy = (arrEMY(i,j,k,EY_LOCAL)-arrEMY(i-1,j,k,EY_LOCAL))/dx[0];
		  Real dxdxBz = computeSecondDerivative(arrEMXY(i-1,j,k,BZ_LOCAL), arrEMXY(i,j,k,BZ_LOCAL), arrEMXY(i+1,j,k,BZ_LOCAL), dx[0]);
		  Real dydyBz = computeSecondDerivative(arrEMXY(i,j-1,k,BZ_LOCAL), arrEMXY(i,j,k,BZ_LOCAL), arrEMXY(i,j+1,k,BZ_LOCAL), dx[1]);
		  rhs(i,j,k) = arrEMXY(i,j,k,BZ_LOCAL) + 0.25*c*c*dt*dt*(dxdxBz+dydyBz) - dt*(dxEy-dyEx);
		  //std::cout << i << " " << j << " " << rhs(i,j,k) << std::endl;		  

		  /*Real dyEx = computeDerivative(arr(i,j-1,k,EX), arr(i,j+1,k,EX), dx[1]);
		  Real dxEy = computeDerivative(arr(i-1,j,k,EY), arr(i+1,j,k,EY), dx[0]);
		  Real dxdxBz = computeSecondDerivative(arr(i-1,j,k,BZ), arr(i,j,k,BZ), arr(i+1,j,k,BZ), dx[0]);
		  Real dydyBz = computeSecondDerivative(arr(i,j-1,k,BZ), arr(i,j,k,BZ), arr(i,j+1,k,BZ), dx[1]);
		  rhs(i,j,k) = arr(i,j,k,BZ) + 0.25*c*c*dt*dt*(dxdxBz+dydyBz) - dt*(dxEy-dyEx);
		  */
		}
	    }
	}
    }

  MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);
  MultiFab& S_EM_XY_int = get_new_data(EM_XY_Type);

  // We need to compute boundary conditions again after each update
  Rhs[0].FillBoundary(geom.periodicity());
  Rhs[1].FillBoundary(geom.periodicity());
  Rhs[2].FillBoundary(geom.periodicity());
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(Rhs[0], geom, {bc_EM[0]});    
  FillDomainBoundary(Rhs[1], geom, {bc_EM[1]});
  FillDomainBoundary(Rhs[2], geom, {bc_EM[2]});

  MultiFab::Copy(S_EM_X_int, Rhs[1], 0, BY_LOCAL, 1, 0);
  FillPatch(*this, Rhs[1], NUM_GROW, time+dt, EM_X_Type, BY_LOCAL, 1);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_int, Rhs[0], 0, BX_LOCAL, 1, 0);
  FillPatch(*this, Rhs[0], NUM_GROW, time+dt, EM_Y_Type, BX_LOCAL, 1);
  MultiFab::Copy(S_EM_XY_int, Rhs[2], 0, BZ_LOCAL, 1, 0);
  FillPatch(*this, Rhs[2], NUM_GROW, time+dt, EM_XY_Type, BZ_LOCAL, 1);
#endif
  
  /*MLMG mlmgX(mlabecX);
  MLMG mlmgY(mlabecY);
  MLMG mlmgZ(mlabecZ);
  
  mlmgX.setMaxFmgIter(max_fmg_iter);
  mlmgX.setVerbose(verbose);
  mlmgY.setMaxFmgIter(max_fmg_iter);
  mlmgY.setVerbose(verbose);
  mlmgZ.setMaxFmgIter(max_fmg_iter);
  mlmgZ.setVerbose(verbose);
  
  const Real S_X_abs = soln_tol*Rhs[0].norm0();
  const Real S_Y_abs = soln_tol*Rhs[1].norm0();
  const Real S_Z_abs = soln_tol*Rhs[2].norm0();  
  */
  // Solve to get S^(n+1)
  //mlmgX.solve({&S_X}, {&Rhs[0]}, soln_tol, S_X_abs);
  //mlmgY.solve({&S_Y}, {&Rhs[1]}, soln_tol, S_Y_abs);
  //mlmgZ.solve({&S_Z}, {&Rhs[2]}, soln_tol, S_Z_abs);

  hypre_solverY.solve(Vector<MultiFab*>{&S_X}, Vector<MultiFab const*>{&Rhs[0]},
		       soln_tol, 0.0, 200);
  hypre_solverX.solve(Vector<MultiFab*>{&S_Y}, Vector<MultiFab const*>{&Rhs[1]},
		       soln_tol, 0.0, 200);
  hypre_solverXY.solve(Vector<MultiFab*>{&S_Z}, Vector<MultiFab const*>{&Rhs[2]},
  		     soln_tol, 0.0, 200);

  // MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
  // MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);
  // MultiFab& S_EM_XY_int = get_new_data(EM_XY_Type);
//   MultiFab::Copy(S_EM_X_int, S_EM_dest[0], BY_LOCAL, BY_LOCAL, 1, 0);
//   FillPatch(*this, S_EM_dest[0], NUM_GROW, time+dt, EM_X_Type, BY_LOCAL, 1);
// #if (AMREX_SPACEDIM >= 2)
//   MultiFab::Copy(S_EM_Y_int, S_EM_dest[1], BX_LOCAL, BX_LOCAL, 1, 0);
//   FillPatch(*this, S_EM_dest[1], NUM_GROW, time+dt, EM_Y_Type, BX_LOCAL, 1);
//   MultiFab::Copy(S_EM_XY_int, S_EM_destEdge, BZ_LOCAL, BZ_LOCAL, 1, 0);
//   FillPatch(*this, S_EM_destEdge, NUM_GROW, time+dt, EM_XY_Type, BZ_LOCAL, 1);
// #endif

  // We need to compute boundary conditions again after each update
  S_EM_dest[0].FillBoundary(geom.periodicity());
  S_EM_dest[1].FillBoundary(geom.periodicity());
  S_EM_destEdge.FillBoundary(geom.periodicity());
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM_dest[0], geom, bc_EM);    
  FillDomainBoundary(S_EM_dest[1], geom, bc_EM);
  FillDomainBoundary(S_EM_destEdge, geom, bc_EM);

  MultiFab::Copy(S_EM_X_int, S_EM_dest[0], 0, 0, 6, 0);
  FillPatch(*this, S_EM_dest[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_int, S_EM_dest[1], 0, 0, 6, 0);
  FillPatch(*this, S_EM_dest[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
  MultiFab::Copy(S_EM_XY_int, S_EM_destEdge, 0, 0, 6, 0);
  FillPatch(*this, S_EM_destEdge, NUM_GROW, time+dt, EM_XY_Type, 0, 6);
  //MultiFab::Copy(S_EM_XY_int, S_EM_destEdge, 0, 0, 6, 0);
#endif

  // MultiFab::Copy(S_dest, S_source, EX, EX, 3, 0);
  
  // for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
  //   {
  //     const Box& bx = mfi.tilebox();
      
  //     const Dim3 lo = lbound(bx);
  //     const Dim3 hi = ubound(bx);

  //     const auto& arr = S_dest.array(mfi);
  //     const auto& arr_old = S_source.array(mfi);

  //     for(int k = lo.z; k <= hi.z; k++)
  //     {
  //       for(int j = lo.y; j <= hi.y; j++)
  //       {
  //         for(int i = lo.x; i <= hi.x; i++)
  //         {	    
  // 	    // update electric field
  // 	    Real dxBy = computeDerivative(arr(i-1,j,k,BY),arr(i+1,j,k,BY),dx[0]);
  // 	    Real dxBz = computeDerivative(arr(i-1,j,k,BZ),arr(i+1,j,k,BZ),dx[0]);
  // 	    /////////////////////////////////////////////////////////////////////////
  // 	    // 2D	    
  // 	    Real dyBx = computeDerivative(arr(i,j-1,k,BX),arr(i,j+1,k,BX),dx[1]);
  // 	    Real dyBz = computeDerivative(arr(i,j-1,k,BZ),arr(i,j+1,k,BZ),dx[1]);
  // 	    // old data
  // 	    Real dxByOld = computeDerivative(arr_old(i-1,j,k,BY),arr_old(i+1,j,k,BY),dx[0]);
  // 	    Real dxBzOld = computeDerivative(arr_old(i-1,j,k,BZ),arr_old(i+1,j,k,BZ),dx[0]);
  // 	    Real dyBxOld = computeDerivative(arr_old(i,j-1,k,BX),arr_old(i,j+1,k,BX),dx[1]);
  // 	    Real dyBzOld = computeDerivative(arr_old(i,j-1,k,BZ),arr_old(i,j+1,k,BZ),dx[1]);
	    
  // 	    /*arr(i,j,k,EX) = arr_old(i,j,k,EX)
  // 	      + dt*(0.5*c*c*dyBz + 0.5*c*c*dyBzOld); 
  // 	    arr(i,j,k,EY) = arr_old(i,j,k,EY)
  // 	    + dt*(- 0.5*c*c*dxBz - 0.5*c*c*dxBzOld);*/
  // 	    arr(i,j,k,EZ) = arr_old(i,j,k,EZ)
  // 	      + dt*(0.5*c*c*(dxBy-dyBx) + 0.5*c*c*(dxByOld-dyBxOld));
  // 	  }
  // 	}
  //     }      
  //   }

  // // We need to compute boundary conditions again after each update 
  // S_dest.FillBoundary(geom.periodicity());
	    
  // // Fill non-periodic physical boundaries                         
  // FillDomainBoundary(S_dest, geom, bc);
  //return;
  /*
  for (int d = 0; d < amrex::SpaceDim ; d++)
    {
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
      
      for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const auto& arr = S_dest.array(mfi);
	  const auto& arr_old = S_source.array(mfi);

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      Real dxBy = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,BX+(1+d)%3),
						    arr(i+iOffset,j+jOffset,k+kOffset,BX+(1+d)%3),dx[d]);
		      Real dxBz = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,BX+(2+d)%3),
						    arr(i+iOffset,j+jOffset,k+kOffset,BX+(2+d)%3),dx[d]);
		      // when using implicit source treatment do not include the current

		      arr(i,j,k,EX+d) += 0.0;
		      arr(i,j,k,EX+(1+d)%3) -= 0.5*dt*c*c*dxBz;
		      arr(i,j,k,EX+(2+d)%3) += 0.5*dt*c*c*dxBy;

		      Real dxByOld = computeDerivative(arr_old(i-iOffset,j-jOffset,k-kOffset,BX+(1+d)%3),
						       arr_old(i+iOffset,j+jOffset,k+kOffset,BX+(1+d)%3),dx[d]);
		      Real dxBzOld = computeDerivative(arr_old(i-iOffset,j-jOffset,k-kOffset,BX+(2+d)%3),
						       arr_old(i+iOffset,j+jOffset,k+kOffset,BX+(2+d)%3),dx[d]);
		      
		      arr(i,j,k,EX+d) += 0.0;
		      arr(i,j,k,EX+(1+d)%3) -= 0.5*dt*c*c*dxBzOld;
		      arr(i,j,k,EX+(2+d)%3) += 0.5*dt*c*c*dxByOld;
			
		    }
		}
	    }
	}
    }

  // We need to compute boundary conditions again after each update
  S_dest.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_dest, geom, bc);  

  return;
  */
  // Update electric field
  for (MFIter mfi(S_EM_dest[0], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arrX = S_EM_dest[0].array(mfi);
      const auto& arrXY = S_EM_destEdge.array(mfi);
      const auto& arrX_old = S_EM_source[0].array(mfi);
      const auto& arrXY_old = S_EM_sourceEdge.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      const Real y = -0.5 + (double(j)+0.5) * dx[1];
	      for(int i = lo.x; i <= hi.x; i++)		
		{
		  const Real x = -0.5 + (double(i)) * dx[0];
		  //Real dyBz = computeDerivative(arrXY(i,j-1,k,BZ_LOCAL),
		  //				arrXY(i,j+1,k,BZ_LOCAL),dx[1]);
		  Real dyBz = (arrXY(i,j+1,k,BZ_LOCAL)-arrXY(i,j,k,BZ_LOCAL))/dx[1];
		  
		  //Real dyBzOld = computeDerivative(arrXY_old(i,j-1,k,BZ_LOCAL),
		  //				   arrXY_old(i,j+1,k,BZ_LOCAL),dx[1]);
		  Real dyBzOld = (arrXY_old(i,j+1,k,BZ_LOCAL)-arrXY_old(i,j,k,BZ_LOCAL))/dx[1];
		  
		  arrX(i,j,k,EX_LOCAL) = arrX_old(i,j,k,EX_LOCAL) + 0.5*dt*c*c*dyBz + 0.5*dt*c*c*dyBzOld;

		}
	    }
	}
    }
  for (MFIter mfi(S_EM_dest[1], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arrY = S_EM_dest[1].array(mfi);
      const auto& arrXY = S_EM_destEdge.array(mfi);
      const auto& arrY_old = S_EM_source[1].array(mfi);
      const auto& arrXY_old = S_EM_sourceEdge.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      const Real y = -0.5 + (double(j)) * dx[1];
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  const Real x = -0.5 + (double(i)+0.5) * dx[0];
		  //Real dxBz = computeDerivative(arrXY(i-1,j,k,BZ_LOCAL),
		  //				arrXY(i+1,j,k,BZ_LOCAL),dx[0]);
		  Real dxBz = (arrXY(i+1,j,k,BZ_LOCAL)-arrXY(i,j,k,BZ_LOCAL))/dx[0];

		  //Real dxBzOld = computeDerivative(arrXY_old(i-1,j,k,BZ_LOCAL),
		  //				   arrXY_old(i+1,j,k,BZ_LOCAL),dx[0]);
		  Real dxBzOld = (arrXY_old(i+1,j,k,BZ_LOCAL)-arrXY_old(i,j,k,BZ_LOCAL))/dx[0];
		  
		  arrY(i,j,k,EY_LOCAL) = arrY_old(i,j,k,EY_LOCAL) - 0.5*dt*c*c*dxBz - 0.5*dt*c*c*dxBzOld;
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

  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arr = S_dest.array(mfi);
      const auto& arr_old = S_source.array(mfi);
      const auto& arrX = S_EM_dest[0].array(mfi);
      const auto& arrY = S_EM_dest[1].array(mfi);
      const auto& arrX_old = S_EM_source[0].array(mfi);
      const auto& arrY_old = S_EM_source[1].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  //Real dxBy = computeDerivative(arrX(i-1,j,k,BY_LOCAL),
		  //				arrX(i+1,j,k,BY_LOCAL),dx[0]);
		  Real dxBy = (arrX(i+1,j,k,BY_LOCAL)-arrX(i,j,k,BY_LOCAL))/dx[0];
		  //Real dyBx = computeDerivative(arrY(i,j-1,k,BX_LOCAL),
		  //				arrY(i,j+1,k,BX_LOCAL),dx[1]);
		  Real dyBx = (arrY(i,j+1,k,BX_LOCAL)-arrY(i,j,k,BX_LOCAL))/dx[1];

		  //Real dxByOld = computeDerivative(arrX_old(i-1,j,k,BY_LOCAL),
		  //				   arrX_old(i+1,j,k,BY_LOCAL),dx[0]);
		  Real dxByOld = (arrX_old(i+1,j,k,BY_LOCAL)-arrX_old(i,j,k,BY_LOCAL))/dx[0];
		  //Real dyBxOld = computeDerivative(arrY_old(i,j-1,k,BX_LOCAL),
		  //				   arrY_old(i,j+1,k,BX_LOCAL),dx[1]);
		  Real dyBxOld = (arrY_old(i,j+1,k,BX_LOCAL)-arrY_old(i,j,k,BX_LOCAL))/dx[1];

		  arr(i,j,k,EZ) = arr_old(i,j,k,EZ) + 0.5*dt*c*c*dxBy + 0.5*dt*c*c*dxByOld - 0.5*dt*c*c*dyBx - 0.5*dt*c*c*dyBxOld;
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

  /*MultiFab::Copy(S_EM_X_int, S_EM_dest[0], 0, 0, 6, 0);
  FillPatch(*this, S_EM_dest[0], NUM_GROW, time+dt, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_int, S_EM_dest[1], 0, 0, 6, 0);
  FillPatch(*this, S_EM_dest[1], NUM_GROW, time+dt, EM_Y_Type, 0, 6);
#endif
  */
  // Compute cell-centred EM fields from Yee-grid EM fields
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrEM_X = S_EM_dest[0].array(mfi);
      const auto& arrEM_Y = S_EM_dest[1].array(mfi);
      const auto& arrEM_XY = S_EM_destEdge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 
		  
  		  arr(i,j,k,BX) = (arrEM_Y(i,j+1,k,BX_LOCAL)+arrEM_Y(i,j,k,BX_LOCAL))/2.0;
  		  arr(i,j,k,BY) = (arrEM_X(i+1,j,k,BY_LOCAL)+arrEM_X(i,j,k,BY_LOCAL))/2.0;
		  arr(i,j,k,BZ) = (arrEM_XY(i+1,j,k,BZ_LOCAL)+arrEM_XY(i,j+1,k,BZ_LOCAL)+
				   arrEM_XY(i+1,j+1,k,BZ_LOCAL)+arrEM_XY(i,j,k,BZ_LOCAL))/4.0;
  		  arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0;
  		  arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0;

  		}
  	    }
  	}       
    }
    
  // We need to compute boundary conditions again after each update
  S_dest.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_dest, geom, bc);  

  //MultiFab::Copy(S_EM_XY_new, S_EM_destEdge, BZ_LOCAL, BZ_LOCAL, 1, 0);
  MultiFab::Copy(S_EM_XY_new, S_EM_destEdge, 0, 0, 6, 0);
}

void CAMReXmp::elecFieldCellAve(Real time)
{

  MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  Array<MultiFab,AMREX_SPACEDIM> S_EM_input;
  S_EM_input[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_input[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2) 
  S_EM_input[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);
#endif

  // Update face-centred EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
    
      // Loop over all the patches at this level
      for (MFIter mfi(S_EM_input[d], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();

	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  // Indexable arrays for the data, and the directional flux
	  // Based on the corner-centred definition of the flux array, the
	  // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	  //const auto& arr = S_dest.array(mfi);
	  const auto& arrEM = S_EM_input[d].array(mfi);
	  const auto& arr = S_input.array(mfi);

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

  MultiFab& S_EM_X_int = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_int = get_new_data(EM_Y_Type);  
  MultiFab::Copy(S_EM_X_int, S_EM_input[0], 0, 0, 6, 0);
  //FillPatch(*this, S_EM_source[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_int, S_EM_input[1], 0, 0, 6, 0);
  //FillPatch(*this, S_EM_source[1], NUM_GROW, time, EM_Y_Type, 0, 6);
#endif

}
void CAMReXmp::MaxwellSolverFVTDVoid(Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& S_source, const Real* dx, Real dt)
{
  amrex::Abort("No Maxwell solver!!");
}
void CAMReXmp::MaxwellSolverFVTDTVD(Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& S_source, const Real* dx, Real dt)
{

  MultiFab& S_dest = get_new_data(Phi_Type);
  MultiFab& S_EM_X_dest = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_dest = get_new_data(EM_Y_Type);
  
  Array<MultiFab,AMREX_SPACEDIM> S_EM_dest;
  S_EM_dest[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  //FillPatch(*this, S_EM_dest[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2) 
  S_EM_dest[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  //FillPatch(*this, S_EM_dest[1], NUM_GROW, time, EM_Y_Type, 0, 6);
#endif

  // Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes[amrex::SpaceDim];
  
  // Define the appropriate size for the flux MultiFab.
  // Fluxes are defined at cell faces - this is taken care of by the
  // surroundingNodes(j) command, ensuring the size of the flux
  for (int j = 0; j < amrex::SpaceDim; j++)
  {
    BoxArray ba = S_source.boxArray();
    ba.surroundingNodes(j);
    fluxes[j].define(ba, dmap, NUM_STATE, 0);
  }

  // Set up a multifab that will contain the fluxes for the electromagnetic fields
  // It will be a nodal multifab (i.e. stored at the corner)
  MultiFab fluxesEM;
  BoxArray ba = S_EM_source[0].boxArray();
  ba.surroundingNodes(1);
  const DistributionMapping& dmX = S_EM_source[0].DistributionMap();
  fluxesEM.define(ba, dmX, 6, 0);

  // Note that fluxesEM and fluxes are different definitions
  // fluxesEM are EM star states, fluxes are fluxes
  
  // Useful when dealing with high-order reconstruction
  // To get first-order slope in y-direction for Bx: slopes[BX_LOCALL].array(i,j,k,Y)
  // Useful to define indices int X,Y,Z in second-order linear reconstruction
  // For third-order quadratic reconstruction use int X,Y,Z,XX,YY,ZZ,XY,YZ,XZ
  // Although for x-components fields, only exist Y,Z,YY,ZZ,YZ (see Balsara et al. 2017)
  
  Array<MultiFab,6> slopes;
  // number of slopes
  // 2 slopes for second-order linear reconstruction
  // Ex = Ex^0 + Ex^y*(y/Delta_y) + Ex^z*(z/Delta_z), where Ex^y and Ex^z are the slopes
  // For 2D code, only 1 slope for x- and y-compoenents
  // For clarity, will use 3 slopes for second-order, first one represents slope in x-direction
  // 9 slopes for third-order, and so on
  int nSlopes = 3;
  // indices for the slopes
  // for 2D, do not need Z index
  // for third order, use also XX,YY,ZZ,XY,YZ,XZ
  int X=0,Y=1,Z=2;
  slopes[BX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);  
  slopes[EX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);
#if (AMREX_SPACEDIM >= 2)
  slopes[BY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
#else  
  // For 1D code it is cell-centred
  slopes[BY_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
  
#if (AMREX_SPACEDIM != 3)  
  // For 2D code it is cell-centred
  slopes[BZ_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EZ_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
    
  // Three components, one for each direction
  MultiFab slopesCharge(grids, dmap, AMREX_SPACEDIM, 1);

  // Compute slopes of the charge densities
  for (int d = 0; d < AMREX_SPACEDIM ; d++)   
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      // Compute cell-centred density slopes
      for (MFIter mfi(slopesCharge, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
 	  const Dim3 hi = ubound(bx);

	  // uses old charge
	  const Array4<Real> arr = S_source.array(mfi);
	  Array4<Real> slopes = slopesCharge.array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y-1; j <= hi.y+1; j++)
		{
		  for(int i = lo.x-1; i <= hi.x+1; i++)
		    {
		      Real u_iMinus1 = (r_i*arr(i-iOffset,j-jOffset,k-kOffset,RHO_I)
					+ r_e*arr(i-iOffset,j-jOffset,k-kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_iPlus1 = (r_i*arr(i+iOffset,j+jOffset,k+kOffset,RHO_I)
				       + r_e*arr(i+iOffset,j+jOffset,k+kOffset,RHO_E))/(lambda_d*lambda_d*l_r);

		      // slopes for the current
		      slopes(i,j,k,d) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		    }
		}
	    }      
	}
    }
  
  // Compute cell-centred z-components slopes in x- and y-direction
  for (int d = 0; d < AMREX_SPACEDIM ; d++)
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
      
      for (MFIter mfi(slopes[BZ_LOCAL], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Dim3 hiDomain = ubound(geom.Domain());

	  const Array4<Real> arr = S_source.array(mfi);
	  Array4<Real> slopesBZ = slopes[BZ_LOCAL].array(mfi);
	  Array4<Real> slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y-1; j <= hi.y+1; j++)
		{
		  for(int i = lo.x-1; i <= hi.x+1; i++)
		    {

		      Real u_iMinus1 = arr(i-iOffset,j-jOffset,k-kOffset,BZ);
		      Real u_i = arr(i,j,k,BZ);
		      Real u_iPlus1 = arr(i+iOffset,j+jOffset,k+kOffset,BZ);
		      slopesBZ(i,j,k,d) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		      u_iMinus1 = arr(i-iOffset,j-jOffset,k-kOffset,EZ);
		      u_i = arr(i,j,k,EZ);
		      u_iPlus1 = arr(i+iOffset,j+jOffset,k+kOffset,EZ);
		      slopesEZ(i,j,k,d) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		    }
		}
	    }      
	}
    }
  
  // Compute y-components slopes in x-direction
  for (MFIter mfi(slopes[BY_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());
      
      const Array4<Real> arr = S_EM_source[1].array(mfi);
      Array4<Real> slopesBY = slopes[BY_LOCAL].array(mfi);
      Array4<Real> slopesEY = slopes[EY_LOCAL].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{

		  Real u_iMinus1 = arr(i-1,j,k,BY_LOCAL);
		  Real u_i = arr(i,j,k,BY_LOCAL);		    
		  Real u_iPlus1 = arr(i+1,j,k,BY_LOCAL);
		  slopesBY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		  u_iMinus1 = arr(i-1,j,k,EY_LOCAL);
		  u_i = arr(i,j,k,EY_LOCAL);		    
		  u_iPlus1 = arr(i+1,j,k,EY_LOCAL);
		  slopesEY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		}
	    }
	}      
    }
  // Compute x-components slopes in y-direction
  for (MFIter mfi(slopes[BX_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());

      const Array4<Real> arr = S_EM_source[0].array(mfi);
      Array4<Real> slopesBX = slopes[BX_LOCAL].array(mfi);
      Array4<Real> slopesEX = slopes[EX_LOCAL].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{

		  Real u_iMinus1 = arr(i,j-1,k,BX_LOCAL);
		  Real u_i = arr(i,j,k,BX_LOCAL);		    
		  Real u_iPlus1 = arr(i,j+1,k,BX_LOCAL);
		  slopesBX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		  u_iMinus1 = arr(i,j-1,k,EX_LOCAL);
		  u_i = arr(i,j,k,EX_LOCAL);		    
		  u_iPlus1 = arr(i,j+1,k,EX_LOCAL);
		  slopesEX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);

		}
	    }
	}      
    }

  // Coefficients for the magnetic and electric fields
  // For second order, there are 21 coeff.
  // i.e. a0,ax,ay,az,axx,...,b0,bx,...,c0,cx,...,czz
  // use 1 ghost cell
  int a0=0,ax=1,ay=2,az=3,axx=4,axy=5,axz=6;
  int b0=7,bx=8,by=9,bz=10,byy=11,bxy=12,byz=13;
  int c0=14,cx=15,cy=16,cz=17,czz=18,cxz=19,cyz=20;
  
  MultiFab Bcoeff(grids, dmap, 21, 1);
  MultiFab Ecoeff(grids, dmap, 21, 1);

  // Compute coefficients for the field function
  for (MFIter mfi(Bcoeff, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      // coeff
      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      // data
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi);       

      // slopes
      const auto& slopesBX = slopes[BX_LOCAL].array(mfi);
      const auto& slopesBY = slopes[BY_LOCAL].array(mfi);
      const auto& slopesBZ = slopes[BZ_LOCAL].array(mfi);
      const auto& slopesEX = slopes[EX_LOCAL].array(mfi);
      const auto& slopesEY = slopes[EY_LOCAL].array(mfi);
      const auto& slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{		 
		  
		  Bc(i,j,k,ay) = (slopesBX(i+1,j,k,Y)+slopesBX(i,j,k,Y))/2.0;
		  Bc(i,j,k,axy) = (slopesBX(i+1,j,k,Y)-slopesBX(i,j,k,Y));
		  Bc(i,j,k,az) = 0.0;
		  Bc(i,j,k,axz) = 0.0;
		  Bc(i,j,k,bx) = (slopesBY(i,j+1,k,X)+slopesBY(i,j,k,X))/2.0;
		  Bc(i,j,k,bxy) = (slopesBY(i,j+1,k,X)-slopesBY(i,j,k,X));
		  Bc(i,j,k,bz) = 0.0;
		  Bc(i,j,k,byz) = 0.0;
		  Bc(i,j,k,cx) = slopesBZ(i,j,k,X);
		  Bc(i,j,k,cxz) = 0.0;
		  Bc(i,j,k,cy) = slopesBZ(i,j,k,Y);
		  Bc(i,j,k,cyz) = 0.0;
		  Bc(i,j,k,axx) = -dx[0]/2.0*(Bc(i,j,k,bxy)/dx[1]); // +cxz/dx[2]
		  Bc(i,j,k,byy) = -dx[1]/2.0*(Bc(i,j,k,axy)/dx[0]); // ++cyz/dx[2]
		  Bc(i,j,k,czz) = 0.0;
		  Bc(i,j,k,a0) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - Bc(i,j,k,axx)/6.0;
		  Bc(i,j,k,ax) = (arrEM_X(i+1,j,k,BX_LOCAL)-arrEM_X(i,j,k,BX_LOCAL));
		  Bc(i,j,k,b0) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - Bc(i,j,k,byy)/6.0;
		  Bc(i,j,k,by) = (arrEM_Y(i,j+1,k,BY_LOCAL)-arrEM_Y(i,j,k,BY_LOCAL));
		  Bc(i,j,k,c0) = arr(i,j,k,BZ) - Bc(i,j,k,czz)/6.0;
		  Bc(i,j,k,cz) = 0.0;

		  Ec(i,j,k,ay) = (slopesEX(i+1,j,k,Y)+slopesEX(i,j,k,Y))/2.0;
		  Ec(i,j,k,axy) = (slopesEX(i+1,j,k,Y)-slopesEX(i,j,k,Y));
		  Ec(i,j,k,az) = 0.0;
		  Ec(i,j,k,axz) = 0.0;
		  Ec(i,j,k,bx) = (slopesEY(i,j+1,k,X)+slopesEY(i,j,k,X))/2.0;
		  Ec(i,j,k,bxy) = (slopesEY(i,j+1,k,X)-slopesEY(i,j,k,X));
		  Ec(i,j,k,bz) = 0.0;
		  Ec(i,j,k,byz) = 0.0;
		  Ec(i,j,k,cx) = slopesEZ(i,j,k,X);
		  Ec(i,j,k,cxz) = 0.0;
		  Ec(i,j,k,cy) = slopesEZ(i,j,k,Y);
		  Ec(i,j,k,cyz) = 0.0;
		  Ec(i,j,k,axx) = -dx[0]/2.0*(Ec(i,j,k,bxy)/dx[1]-slopesQ(i,j,k,0)); // +cxz/dx[2]
		  Ec(i,j,k,byy) = -dx[1]/2.0*(Ec(i,j,k,axy)/dx[0]-slopesQ(i,j,k,1)); // ++cyz/dx[2]
		  Ec(i,j,k,czz) = 0.0;
		  Ec(i,j,k,a0) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - Ec(i,j,k,axx)/6.0;
		  Ec(i,j,k,ax) = (arrEM_X(i+1,j,k,EX_LOCAL)-arrEM_X(i,j,k,EX_LOCAL));
		  Ec(i,j,k,b0) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - Ec(i,j,k,byy)/6.0;
		  Ec(i,j,k,by) = (arrEM_Y(i,j+1,k,EY_LOCAL)-arrEM_Y(i,j,k,EY_LOCAL));
		  Ec(i,j,k,c0) = arr(i,j,k,EZ) - Ec(i,j,k,czz)/6.0; //czz*dx[2]*dx[2]/6.0
		  Ec(i,j,k,cz) = 0.0;		  

		}
	    }
	}
    }

  // Compute edge components of the EM fields    
  for (MFIter mfi(fluxesEM, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      // use old array, S_dest should also work
      //const auto& arr = S_dest.array(mfi);
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi); 
      const auto& fluxArrEM = fluxesEM.array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);


      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{		 
		  
		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
		  // LD state
		  x = dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = -1, jOffset = -1;
		  Vector<Real> EM_LD = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  		  
		  // RD state
		  x = -dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_RD = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // LU state
		  x = dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_LU = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // RU state
		  x = -dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_RU = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // If the reconstruction is consistent, then EyRU=EyRD
		  Real EyR = (EM_RU[EY_LOCAL]+EM_RD[EY_LOCAL])/2.0;
		  Real EyL = (EM_LU[EY_LOCAL]+EM_LD[EY_LOCAL])/2.0;
		  Real ExU = (EM_RU[EX_LOCAL]+EM_LU[EX_LOCAL])/2.0;
		  Real ExD = (EM_RD[EX_LOCAL]+EM_LD[EX_LOCAL])/2.0;
		  fluxArrEM(i,j,k,BZ_LOCAL) = 0.25*(EM_RU[BZ_LOCAL]+EM_RD[BZ_LOCAL]+EM_LU[BZ_LOCAL]+EM_LD[BZ_LOCAL])
		    - 0.5/c*(EyR-EyL) + 0.5/c*(ExU-ExD);

		  Real ByR = (EM_RU[BY_LOCAL]+EM_RD[BY_LOCAL])/2.0;
		  Real ByL = (EM_LU[BY_LOCAL]+EM_LD[BY_LOCAL])/2.0;
		  Real BxU = (EM_RU[BX_LOCAL]+EM_LU[BX_LOCAL])/2.0;
		  Real BxD = (EM_RD[BX_LOCAL]+EM_LD[BX_LOCAL])/2.0;
		  fluxArrEM(i,j,k,EZ_LOCAL) = 0.25*(EM_RU[EZ_LOCAL]+EM_RD[EZ_LOCAL]+EM_LU[EZ_LOCAL]+EM_LD[EZ_LOCAL])
		    + 0.5*c*(ByR-ByL) - 0.5*c*(BxU-BxD);
		}
	    }
	}       
    }

  // Update face-centred EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      int d_EM = (d==0) ? 1 : 0;
    
      // Loop over all the patches at this level
      for (MFIter mfi(S_EM_source[d_EM], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();

	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  // Indexable arrays for the data, and the directional flux
	  // Based on the corner-centred definition of the flux array, the
	  // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	  //const auto& arr = S_dest.array(mfi);
	  const auto& arrEM = S_EM_dest[d_EM].array(mfi);
	  const auto& arrEMOld = S_EM_source[d_EM].array(mfi);
	  const auto& fluxArrEM = fluxesEM.array(mfi);
	  const auto& fluxArr = fluxes[d_EM].array(mfi);

	  const Dim3 hiDomain = ubound(geom.Domain());
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      // 2D code; z-component updated using cell-centred scheme
		      arrEM(i,j,k,BX_LOCAL+d_EM) = arrEMOld(i,j,k,BX_LOCAL+d_EM) + std::pow(-1,d)*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EZ_LOCAL) - fluxArrEM(i,j,k,EZ_LOCAL));
		      arrEM(i,j,k,EX_LOCAL+d_EM) = arrEMOld(i,j,k,EX_LOCAL+d_EM) - std::pow(-1,d)*c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) - fluxArrEM(i,j,k,BZ_LOCAL));

		      // 3D code
		      /*arrEM(i,j,k,BX_LOCAL+(1+d)%3) += (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(2+d)%3));
			arrEM(i,j,k,BX_LOCAL+(2+d)%3) -= (dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,EX_LOCAL+(1+d)%3));
			arrEM(i,j,k,EX_LOCAL+(1+d)%3) -= c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(2+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(2+d)%3));
			arrEM(i,j,k,EX_LOCAL+(2+d)%3) += c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BX_LOCAL+(1+d)%3) - fluxArrEM(i,j,k,BX_LOCAL+(1+d)%3));
		      */
		      
		    }
		}
	    }      
	}
    
    }

  // We need to compute boundary conditions again after each update
  /*S_EM_dest[0].FillBoundary(geom.periodicity());
  S_EM_dest[1].FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM_dest[0], geom, bc_EM);    
  FillDomainBoundary(S_EM_dest[1], geom, bc_EM);*/

  MultiFab::Copy(S_EM_X_dest, S_EM_dest[0], 0, 0, 6, 0);
  MultiFab::Copy(S_EM_Y_dest, S_EM_dest[1], 0, 0, 6, 0);
  
  // Compute cell-centred EM fields from face-centred
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrEM_X = S_EM_dest[0].array(mfi);
      const auto& arrEM_Y = S_EM_dest[1].array(mfi); 

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

  		  arr(i,j,k,BX) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0;// - Bc(i,j,k,axx)/6.0;
  		  arr(i,j,k,BY) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0;// - Bc(i,j,k,byy)/6.0;
  		  arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0;// - Ec(i,j,k,axx)/6.0;
  		  arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0;// - Ec(i,j,k,byy)/6.0;
		  
  		}
  	    }
  	}       
    }

  // Only 2D
  // Compute edge components for y-components at x-faces
  for (MFIter mfi(fluxes[0], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      const auto& fluxArrX = fluxes[0].array(mfi);
      //const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // L state
  		  x = dx[0]/2.0, y = 0.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_L = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		   
  		  // R state
  		  x = -dx[0]/2.0, y = 0.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_R = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // Note that this is not the flux, but the star states
  		  fluxArrX(i,j,k,BY) = 0.5*(EM_R[BY_LOCAL]+EM_L[BY_LOCAL])
		    + 0.5/c*(EM_R[EZ_LOCAL]-EM_L[EZ_LOCAL]);
  		  fluxArrX(i,j,k,EY) = 0.5*(EM_R[EY_LOCAL]+EM_L[EY_LOCAL])
		    - 0.5*c*(EM_R[BZ_LOCAL]-EM_L[BZ_LOCAL]);
  		}
  	    }
  	}
    }
  
  // Only 2D
  // Compute edge components for x-components at y-faces
  for (MFIter mfi(fluxes[1], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 		  

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;
		  	
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;

  		  // D state
  		  x = 0.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_D = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
  		  // U state
  		  x = 0.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_U = EM_linearFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // Note that this is not the flux, but the star states		  
  		  fluxArrY(i,j,k,BX) = 0.5*(EM_U[BX_LOCAL]+EM_D[BX_LOCAL])
		    - 0.5/c*(EM_U[EZ_LOCAL]-EM_D[EZ_LOCAL]);
  		  fluxArrY(i,j,k,EX) = 0.5*(EM_U[EX_LOCAL]+EM_D[EX_LOCAL])
		    + 0.5*c*(EM_U[BZ_LOCAL]-EM_D[BZ_LOCAL]); 
  		}
  	    }
  	}
    }

  // Update cell-centred z-components of EM fields
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrOld = S_source.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);
	        
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  // Update cell-centred z-components becuause it is 2D code
  		  arr(i,j,k,BZ) = arrOld(i,j,k,BZ) - (dt / dx[0]) * (fluxArrX(i+1,j,k,EY) - fluxArrX(i,j,k,EY)) + (dt / dx[1]) * (fluxArrY(i,j+1,k,EX) - fluxArrY(i,j,k,EX));
  		  arr(i,j,k,EZ) = arrOld(i,j,k,EZ) + c*c*(dt / dx[0]) * (fluxArrX(i+1,j,k,BY) - fluxArrX(i,j,k,BY)) - c*c*(dt / dx[1]) * (fluxArrY(i,j+1,k,BX) - fluxArrY(i,j,k,BX));
  		}
  	    }
  	}
    }
  // We need to compute boundary conditions again after each update
  /*S_dest.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_dest, geom, bc);  */
}
void CAMReXmp::MaxwellSolverFDTDCN(const Real* dx, Real dt, Real time) 
{

  // get multifabs references
  MultiFab& S_new = get_new_data(Phi_Type);
  MultiFab& S_EM_X_new = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_new = get_new_data(EM_Y_Type);
  MultiFab& S_EM_XY_new = get_new_data(EM_XY_Type);

  // input states and fill the data
  MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
  MultiFab S_output(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  FillPatch(*this, S_output, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  Array<MultiFab,AMREX_SPACEDIM> S_EM_input, S_EM_output;
  S_EM_input[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_input[0], NUM_GROW, time, EM_X_Type, 0, 6);
  S_EM_output[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_output[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2) 
  S_EM_input[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);
  S_EM_output[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_output[1], NUM_GROW, time, EM_Y_Type, 0, 6);
#endif

  MultiFab S_EM_edge_input, S_EM_edge_output;
  S_EM_edge_input.define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 6, NUM_GROW);
  S_EM_edge_output.define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_edge_input, NUM_GROW, time, EM_XY_Type, 0, 6);
  FillPatch(*this, S_EM_edge_output, NUM_GROW, time, EM_XY_Type, 0, 6);

  IndexType xface(IntVect{AMREX_D_DECL(1,0,0)});
  IndexType yface(IntVect{AMREX_D_DECL(0,1,0)});
  IndexType edge(IntVect{AMREX_D_DECL(1,1,0)});

  auto const& nddom = amrex::surroundingNodes(geom.Domain());
  auto const& nddomX = amrex::surroundingNodes(geom.Domain(),0);
  auto const& nddomY = amrex::surroundingNodes(geom.Domain(),1);
   const auto dxi = geom.InvCellSizeArray();

  GpuArray<HYPRE_Real,AMREX_SPACEDIM> fac
    {AMREX_D_DECL(static_cast<HYPRE_Real>(dxi[0]*dxi[0]),
  		  static_cast<HYPRE_Real>(dxi[1]*dxi[1]),
  		  static_cast<HYPRE_Real>(dxi[2]*dxi[2]))};


  HYPRE_Real fac0 = HYPRE_Real(-2.)*(AMREX_D_TERM(fac[0],+fac[1],+fac[2]));

   if (bc_EM[BX_LOCAL].lo(1)==REFLECT_ODD || bc_EM[BY_LOCAL].lo(0)==REFLECT_ODD)
     amrex::Abort("Currently does not support fillerFunction reflec_odd for Bx in y-direction or By in x-direction");

   //functor that returns whether the variable n at (i,j,k) in Box boxno (local index) is valid
   auto marker = markerFunction(bc_EM[BZ_LOCAL],nddom);
   auto markerX = markerFunction(bc_EM[BY_LOCAL],nddomX);
   auto markerY = markerFunction(bc_EM[BX_LOCAL],nddomY);
 
  // functor that fills the row in the matrix A for variable n at (i,j,k) in Box boxno (local index)
  // using the CSR format
  auto filler = fillerFunction(bc_EM[BZ_LOCAL],nddom,dxi,dt);
  auto fillerX = fillerFunction(bc_EM[BY_LOCAL],nddomX,dxi,dt);
  auto fillerY = fillerFunction(bc_EM[BX_LOCAL],nddomY,dxi,dt);
  
  constexpr int max_stencil_size = 2*AMREX_SPACEDIM+1;
  HypreSolver<max_stencil_size> hypre_solverX
    ({xface}, IntVect(1), geom, grids, dmap,
     markerX, fillerX, 2);
  HypreSolver<max_stencil_size> hypre_solverY
    ({yface}, IntVect(1), geom, grids, dmap,
     markerY, fillerY, 2);
  HypreSolver<max_stencil_size> hypre_solverXY
    //({IndexType::TheNodeType()}, IntVect(1), geom, grids, dmap,
    ({edge}, IntVect(1), geom, grids, dmap,
     marker, filler, 2);
      
  MultiFab S_X(S_EM_output[1], amrex::make_alias, BX_LOCAL, 1);
  MultiFab S_Y(S_EM_output[0], amrex::make_alias, BY_LOCAL, 1);
  MultiFab S_Z(S_EM_edge_output, amrex::make_alias, BZ_LOCAL, 1);  
  
  std::array<MultiFab,3> Rhs;
  Rhs[0].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 1, NUM_GROW);
  Rhs[1].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 1, NUM_GROW);
  Rhs[2].define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 1, NUM_GROW);

  for(MFIter mfi(Rhs[0], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs[0].array(mfi);
      const auto& arr = S_input.array(mfi);
      const auto& arrEM = S_EM_input[1].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dyEz = (arr(i,j,k,EZ)-arr(i,j-1,k,EZ))/dx[1];
		  Real dxdxBx = computeSecondDerivative(arrEM(i-1,j,k,BX_LOCAL), arrEM(i,j,k,BX_LOCAL), arrEM(i+1,j,k,BX_LOCAL), dx[0]);
		  Real dydyBx = computeSecondDerivative(arrEM(i,j-1,k,BX_LOCAL), arrEM(i,j,k,BX_LOCAL), arrEM(i,j+1,k,BX_LOCAL), dx[1]);
		  rhs(i,j,k) = arrEM(i,j,k,BX_LOCAL) + 0.25*c*c*dt*dt*(dxdxBx+dydyBx) - dt*dyEz;
		}
	    }
	}
    }
  for(MFIter mfi(Rhs[1], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
	  
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs[1].array(mfi);
      const auto& arr = S_input.array(mfi);
      const auto& arrEM = S_EM_input[0].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dxEz = (arr(i,j,k,EZ)-arr(i-1,j,k,EZ))/dx[0];
		  Real dxdxBy = computeSecondDerivative(arrEM(i-1,j,k,BY_LOCAL), arrEM(i,j,k,BY_LOCAL), arrEM(i+1,j,k,BY_LOCAL), dx[0]);
		  Real dydyBy = computeSecondDerivative(arrEM(i,j-1,k,BY_LOCAL), arrEM(i,j,k,BY_LOCAL), arrEM(i,j+1,k,BY_LOCAL), dx[1]);
		  rhs(i,j,k) = arrEM(i,j,k,BY_LOCAL) + 0.25*c*c*dt*dt*(dxdxBy+dydyBy) + dt*dxEz;
		}
	    }
	}
    }
  for(MFIter mfi(Rhs[2], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
	  
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs[2].array(mfi);
      const auto& arr = S_input.array(mfi);
      const auto& arrEMX = S_EM_input[0].array(mfi);
      const auto& arrEMY = S_EM_input[1].array(mfi);
      const auto& arrEMXY = S_EM_edge_input.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dyEx = (arrEMX(i,j,k,EX_LOCAL)-arrEMX(i,j-1,k,EX_LOCAL))/dx[1];
		  Real dxEy = (arrEMY(i,j,k,EY_LOCAL)-arrEMY(i-1,j,k,EY_LOCAL))/dx[0];
		  Real dxdxBz = computeSecondDerivative(arrEMXY(i-1,j,k,BZ_LOCAL), arrEMXY(i,j,k,BZ_LOCAL), arrEMXY(i+1,j,k,BZ_LOCAL), dx[0]);
		  Real dydyBz = computeSecondDerivative(arrEMXY(i,j-1,k,BZ_LOCAL), arrEMXY(i,j,k,BZ_LOCAL), arrEMXY(i,j+1,k,BZ_LOCAL), dx[1]);
		  rhs(i,j,k) = arrEMXY(i,j,k,BZ_LOCAL) + 0.25*c*c*dt*dt*(dxdxBz+dydyBz) - dt*(dxEy-dyEx);
		}
	    }
	}
    }

  // We need to compute boundary conditions again after each update
  /*Rhs[0].FillBoundary(geom.periodicity());
  Rhs[1].FillBoundary(geom.periodicity());
  Rhs[2].FillBoundary(geom.periodicity());
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(Rhs[0], geom, {bc_EM[0]});    
  FillDomainBoundary(Rhs[1], geom, {bc_EM[1]});
  FillDomainBoundary(Rhs[2], geom, {bc_EM[2]});
  */
  MultiFab::Copy(S_EM_X_new, Rhs[1], 0, BY_LOCAL, 1, 0);
  FillPatch(*this, Rhs[1], NUM_GROW, time, EM_X_Type, BY_LOCAL, 1);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_new, Rhs[0], 0, BX_LOCAL, 1, 0);
  FillPatch(*this, Rhs[0], NUM_GROW, time, EM_Y_Type, BX_LOCAL, 1);
  MultiFab::Copy(S_EM_XY_new, Rhs[2], 0, BZ_LOCAL, 1, 0);
  FillPatch(*this, Rhs[2], NUM_GROW, time, EM_XY_Type, BZ_LOCAL, 1);
#endif

  hypre_solverY.solve(Vector<MultiFab*>{&S_X}, Vector<MultiFab const*>{&Rhs[0]},
		       soln_tol, 0.0, 200);
  hypre_solverX.solve(Vector<MultiFab*>{&S_Y}, Vector<MultiFab const*>{&Rhs[1]},
		       soln_tol, 0.0, 200);
  hypre_solverXY.solve(Vector<MultiFab*>{&S_Z}, Vector<MultiFab const*>{&Rhs[2]},
  		     soln_tol, 0.0, 200);
  
  // We need to compute boundary conditions again after each update
  /*S_EM_output[0].FillBoundary(geom.periodicity());
  S_EM_output[1].FillBoundary(geom.periodicity());
  S_EM_edge_output.FillBoundary(geom.periodicity());
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM_output[0], geom, bc_EM);    
  FillDomainBoundary(S_EM_output[1], geom, bc_EM);
  FillDomainBoundary(S_EM_edge_output, geom, bc_EM);
  */
  MultiFab::Copy(S_EM_X_new, S_EM_output[0], 0, 0, 6, 0);
  FillPatch(*this, S_EM_output[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_new, S_EM_output[1], 0, 0, 6, 0);
  FillPatch(*this, S_EM_output[1], NUM_GROW, time, EM_Y_Type, 0, 6);
  MultiFab::Copy(S_EM_XY_new, S_EM_edge_output, 0, 0, 6, 0);
  FillPatch(*this, S_EM_edge_output, NUM_GROW, time, EM_XY_Type, 0, 6);
#endif

  // Update electric field
  for (MFIter mfi(S_EM_output[0], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arrX = S_EM_output[0].array(mfi);
      const auto& arrXY = S_EM_edge_output.array(mfi);
      const auto& arrX_old = S_EM_input[0].array(mfi);
      const auto& arrXY_old = S_EM_edge_input.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      const Real y = -0.5 + (double(j)+0.5) * dx[1];
	      for(int i = lo.x; i <= hi.x; i++)		
		{
		  const Real x = -0.5 + (double(i)) * dx[0];
		  Real dyBz = (arrXY(i,j+1,k,BZ_LOCAL)-arrXY(i,j,k,BZ_LOCAL))/dx[1];		  
		  Real dyBzOld = (arrXY_old(i,j+1,k,BZ_LOCAL)-arrXY_old(i,j,k,BZ_LOCAL))/dx[1];
		  
		  arrX(i,j,k,EX_LOCAL) = arrX_old(i,j,k,EX_LOCAL) + 0.5*dt*c*c*dyBz + 0.5*dt*c*c*dyBzOld;

		}
	    }
	}
    }
  for (MFIter mfi(S_EM_output[1], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arrY = S_EM_output[1].array(mfi);
      const auto& arrXY = S_EM_edge_output.array(mfi);
      const auto& arrY_old = S_EM_input[1].array(mfi);
      const auto& arrXY_old = S_EM_edge_input.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      const Real y = -0.5 + (double(j)) * dx[1];
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  const Real x = -0.5 + (double(i)+0.5) * dx[0];
		  Real dxBz = (arrXY(i+1,j,k,BZ_LOCAL)-arrXY(i,j,k,BZ_LOCAL))/dx[0];
		  Real dxBzOld = (arrXY_old(i+1,j,k,BZ_LOCAL)-arrXY_old(i,j,k,BZ_LOCAL))/dx[0];
		  
		  arrY(i,j,k,EY_LOCAL) = arrY_old(i,j,k,EY_LOCAL) - 0.5*dt*c*c*dxBz - 0.5*dt*c*c*dxBzOld;
		}
	    }
	}
    }

  // We need to compute boundary conditions again after each update
  /*S_EM_output[0].FillBoundary(geom.periodicity());
  S_EM_output[1].FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM_output[0], geom, bc_EM);    
  FillDomainBoundary(S_EM_output[1], geom, bc_EM);
  */
  for (MFIter mfi(S_output, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arr = S_output.array(mfi);
      const auto& arr_old = S_input.array(mfi);
      const auto& arrX = S_EM_output[0].array(mfi);
      const auto& arrY = S_EM_output[1].array(mfi);
      const auto& arrX_old = S_EM_input[0].array(mfi);
      const auto& arrY_old = S_EM_input[1].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dxBy = (arrX(i+1,j,k,BY_LOCAL)-arrX(i,j,k,BY_LOCAL))/dx[0];
		  Real dyBx = (arrY(i,j+1,k,BX_LOCAL)-arrY(i,j,k,BX_LOCAL))/dx[1];
		  Real dxByOld = (arrX_old(i+1,j,k,BY_LOCAL)-arrX_old(i,j,k,BY_LOCAL))/dx[0];
		  Real dyBxOld = (arrY_old(i,j+1,k,BX_LOCAL)-arrY_old(i,j,k,BX_LOCAL))/dx[1];

		  arr(i,j,k,EZ) = arr_old(i,j,k,EZ) + 0.5*dt*c*c*dxBy + 0.5*dt*c*c*dxByOld - 0.5*dt*c*c*dyBx - 0.5*dt*c*c*dyBxOld;
		}
	    }
	}
    }

  // We need to compute boundary conditions again after each update
  /*S_EM_output[0].FillBoundary(geom.periodicity());
  S_EM_output[1].FillBoundary(geom.periodicity());
  
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM_output[0], geom, bc_EM);    
  FillDomainBoundary(S_EM_output[1], geom, bc_EM);
  */
  MultiFab::Copy(S_EM_X_new, S_EM_output[0], 0, 0, 6, 0);
  FillPatch(*this, S_EM_output[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_new, S_EM_output[1], 0, 0, 6, 0);
  FillPatch(*this, S_EM_output[1], NUM_GROW, time, EM_Y_Type, 0, 6);
#endif

  // Compute cell-centred EM fields from Yee-grid EM fields
  for (MFIter mfi(S_output, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = S_output.array(mfi);
      const auto& arrEM_X = S_EM_output[0].array(mfi);
      const auto& arrEM_Y = S_EM_output[1].array(mfi);
      const auto& arrEM_XY = S_EM_edge_output.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 
		  
  		  arr(i,j,k,BX) = (arrEM_Y(i,j+1,k,BX_LOCAL)+arrEM_Y(i,j,k,BX_LOCAL))/2.0;
  		  arr(i,j,k,BY) = (arrEM_X(i+1,j,k,BY_LOCAL)+arrEM_X(i,j,k,BY_LOCAL))/2.0;
		  arr(i,j,k,BZ) = (arrEM_XY(i+1,j,k,BZ_LOCAL)+arrEM_XY(i,j+1,k,BZ_LOCAL)+
				   arrEM_XY(i+1,j+1,k,BZ_LOCAL)+arrEM_XY(i,j,k,BZ_LOCAL))/4.0;
  		  arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0;
  		  arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0;

  		}
  	    }
  	}       
    }
    
  /*// We need to compute boundary conditions again after each update
  S_output.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_output, geom, bc);  
  */
  
  MultiFab::Copy(S_new, S_output, BX, BX, 6, 0);
  MultiFab::Copy(S_EM_XY_new, S_EM_edge_output, 0, 0, 6, 0);
}
void CAMReXmp::MaxwellSolverFDTDCNAMReX(const Real* dx, Real dt, Real time) 
{

  // get multifabs references
  MultiFab& S_new = get_new_data(Phi_Type);
  MultiFab& S_EM_X_new = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_new = get_new_data(EM_Y_Type);
  MultiFab& S_EM_XY_new = get_new_data(EM_XY_Type);

  // input states and fill the data
  MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
  MultiFab S_output(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  FillPatch(*this, S_output, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  Array<MultiFab,AMREX_SPACEDIM> S_EM_input, S_EM_output;
  S_EM_input[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_input[0], NUM_GROW, time, EM_X_Type, 0, 6);
  S_EM_output[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_output[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2) 
  S_EM_input[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);
  S_EM_output[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_output[1], NUM_GROW, time, EM_Y_Type, 0, 6);
#endif

  MultiFab S_EM_edge_input, S_EM_edge_output;
  S_EM_edge_input.define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 6, NUM_GROW);
  S_EM_edge_output.define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 6, NUM_GROW);
  FillPatch(*this, S_EM_edge_input, NUM_GROW, time, EM_XY_Type, 0, 6);
  FillPatch(*this, S_EM_edge_output, NUM_GROW, time, EM_XY_Type, 0, 6);

  IndexType xface(IntVect{AMREX_D_DECL(1,0,0)});
  IndexType yface(IntVect{AMREX_D_DECL(0,1,0)});
  IndexType edge(IntVect{AMREX_D_DECL(1,1,0)});

  auto const& nddom = amrex::surroundingNodes(geom.Domain());
  auto const& nddomX = amrex::surroundingNodes(geom.Domain(),0);
  auto const& nddomY = amrex::surroundingNodes(geom.Domain(),1);
   const auto dxi = geom.InvCellSizeArray();

  GpuArray<HYPRE_Real,AMREX_SPACEDIM> fac
    {AMREX_D_DECL(static_cast<HYPRE_Real>(dxi[0]*dxi[0]),
  		  static_cast<HYPRE_Real>(dxi[1]*dxi[1]),
  		  static_cast<HYPRE_Real>(dxi[2]*dxi[2]))};


  HYPRE_Real fac0 = HYPRE_Real(-2.)*(AMREX_D_TERM(fac[0],+fac[1],+fac[2]));

   if (bc_EM[BX_LOCAL].lo(1)==REFLECT_ODD || bc_EM[BY_LOCAL].lo(0)==REFLECT_ODD)
     amrex::Abort("Currently does not support fillerFunction reflec_odd for Bx in y-direction or By in x-direction");

   //functor that returns whether the variable n at (i,j,k) in Box boxno (local index) is valid
   auto marker = markerFunction(bc_EM[BZ_LOCAL],nddom);
   auto markerX = markerFunction(bc_EM[BY_LOCAL],nddomX);
   auto markerY = markerFunction(bc_EM[BX_LOCAL],nddomY);
 
  // functor that fills the row in the matrix A for variable n at (i,j,k) in Box boxno (local index)
  // using the CSR format
  auto filler = fillerFunction(bc_EM[BZ_LOCAL],nddom,dxi,dt);
  auto fillerX = fillerFunction(bc_EM[BY_LOCAL],nddomX,dxi,dt);
  auto fillerY = fillerFunction(bc_EM[BX_LOCAL],nddomY,dxi,dt);
  
  constexpr int max_stencil_size = 2*AMREX_SPACEDIM+1;
  HypreSolver<max_stencil_size> hypre_solverX
    ({xface}, IntVect(1), geom, grids, dmap,
     markerX, fillerX, 2);
  HypreSolver<max_stencil_size> hypre_solverY
    ({yface}, IntVect(1), geom, grids, dmap,
     markerY, fillerY, 2);
  HypreSolver<max_stencil_size> hypre_solverXY
    //({IndexType::TheNodeType()}, IntVect(1), geom, grids, dmap,
    ({edge}, IntVect(1), geom, grids, dmap,
     marker, filler, 2);
      
  MultiFab S_X(S_EM_output[1], amrex::make_alias, BX_LOCAL, 1);
  MultiFab S_Y(S_EM_output[0], amrex::make_alias, BY_LOCAL, 1);
  MultiFab S_Z(S_EM_edge_output, amrex::make_alias, BZ_LOCAL, 1);  
  
  std::array<MultiFab,3> Rhs;
  Rhs[0].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 1, NUM_GROW);
  Rhs[1].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 1, NUM_GROW);
  Rhs[2].define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 1, NUM_GROW);

  for(MFIter mfi(Rhs[0], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs[0].array(mfi);
      const auto& arr = S_input.array(mfi);
      const auto& arrEM = S_EM_input[1].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dyEz = (arr(i,j,k,EZ)-arr(i,j-1,k,EZ))/dx[1];
		  Real dxdxBx = computeSecondDerivative(arrEM(i-1,j,k,BX_LOCAL), arrEM(i,j,k,BX_LOCAL), arrEM(i+1,j,k,BX_LOCAL), dx[0]);
		  Real dydyBx = computeSecondDerivative(arrEM(i,j-1,k,BX_LOCAL), arrEM(i,j,k,BX_LOCAL), arrEM(i,j+1,k,BX_LOCAL), dx[1]);
		  rhs(i,j,k) = arrEM(i,j,k,BX_LOCAL) + 0.25*c*c*dt*dt*(dxdxBx+dydyBx) - dt*dyEz;
		}
	    }
	}
    }
  for(MFIter mfi(Rhs[1], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
	  
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs[1].array(mfi);
      const auto& arr = S_input.array(mfi);
      const auto& arrEM = S_EM_input[0].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dxEz = (arr(i,j,k,EZ)-arr(i-1,j,k,EZ))/dx[0];
		  Real dxdxBy = computeSecondDerivative(arrEM(i-1,j,k,BY_LOCAL), arrEM(i,j,k,BY_LOCAL), arrEM(i+1,j,k,BY_LOCAL), dx[0]);
		  Real dydyBy = computeSecondDerivative(arrEM(i,j-1,k,BY_LOCAL), arrEM(i,j,k,BY_LOCAL), arrEM(i,j+1,k,BY_LOCAL), dx[1]);
		  rhs(i,j,k) = arrEM(i,j,k,BY_LOCAL) + 0.25*c*c*dt*dt*(dxdxBy+dydyBy) + dt*dxEz;
		}
	    }
	}
    }
  for(MFIter mfi(Rhs[2], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
	  
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& rhs = Rhs[2].array(mfi);
      const auto& arr = S_input.array(mfi);
      const auto& arrEMX = S_EM_input[0].array(mfi);
      const auto& arrEMY = S_EM_input[1].array(mfi);
      const auto& arrEMXY = S_EM_edge_input.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dyEx = (arrEMX(i,j,k,EX_LOCAL)-arrEMX(i,j-1,k,EX_LOCAL))/dx[1];
		  Real dxEy = (arrEMY(i,j,k,EY_LOCAL)-arrEMY(i-1,j,k,EY_LOCAL))/dx[0];
		  Real dxdxBz = computeSecondDerivative(arrEMXY(i-1,j,k,BZ_LOCAL), arrEMXY(i,j,k,BZ_LOCAL), arrEMXY(i+1,j,k,BZ_LOCAL), dx[0]);
		  Real dydyBz = computeSecondDerivative(arrEMXY(i,j-1,k,BZ_LOCAL), arrEMXY(i,j,k,BZ_LOCAL), arrEMXY(i,j+1,k,BZ_LOCAL), dx[1]);
		  rhs(i,j,k) = arrEMXY(i,j,k,BZ_LOCAL) + 0.25*c*c*dt*dt*(dxdxBz+dydyBz) - dt*(dxEy-dyEx);
		}
	    }
	}
    }

  // We need to compute boundary conditions again after each update
  /*Rhs[0].FillBoundary(geom.periodicity());
  Rhs[1].FillBoundary(geom.periodicity());
  Rhs[2].FillBoundary(geom.periodicity());
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(Rhs[0], geom, {bc_EM[0]});    
  FillDomainBoundary(Rhs[1], geom, {bc_EM[1]});
  FillDomainBoundary(Rhs[2], geom, {bc_EM[2]});
  */
  MultiFab::Copy(S_EM_X_new, Rhs[1], 0, BY_LOCAL, 1, 0);
  FillPatch(*this, Rhs[1], NUM_GROW, time, EM_X_Type, BY_LOCAL, 1);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_new, Rhs[0], 0, BX_LOCAL, 1, 0);
  FillPatch(*this, Rhs[0], NUM_GROW, time, EM_Y_Type, BX_LOCAL, 1);
  MultiFab::Copy(S_EM_XY_new, Rhs[2], 0, BZ_LOCAL, 1, 0);
  FillPatch(*this, Rhs[2], NUM_GROW, time, EM_XY_Type, BZ_LOCAL, 1);
#endif

  hypre_solverY.solve(Vector<MultiFab*>{&S_X}, Vector<MultiFab const*>{&Rhs[0]},
		       soln_tol, 0.0, 200);
  hypre_solverX.solve(Vector<MultiFab*>{&S_Y}, Vector<MultiFab const*>{&Rhs[1]},
		       soln_tol, 0.0, 200);
  //hypre_solverXY.solve(Vector<MultiFab*>{&S_Z}, Vector<MultiFab const*>{&Rhs[2]},
  //		     soln_tol, 0.0, 200);

  bool agglomeration = true;
  bool consolidation = true;
  int max_coarsening_level = 30;
  
  LPInfo info;
  info.setAgglomeration(agglomeration);
  info.setConsolidation(consolidation);
  info.setMaxCoarseningLevel(max_coarsening_level);

  //const auto tol_rel = Real(1.e-10);
  const auto tol_rel = Real(soln_tol);
  const auto tol_abs = Real(0.0);

  MLNodeABecLaplacian mlndabec({geom}, {grids}, {dmap}, info);

  setDomainBC(mlmg_lobc_Z, mlmg_hibc_Z, BZ);
  mlndabec.setDomainBC(mlmg_lobc_Z, mlmg_hibc_Z);

  Real ascalar = 1.0;
  Real bscalar = 0.5*0.5*c*c*dt*dt;
  mlndabec.setScalars(ascalar, bscalar);

  mlndabec.setACoeffs(0, 1.0);
  mlndabec.setBCoeffs(0, 1.0);

  int max_iter = 100;
  int bottom_verbose = 0;
  
  MLMG mlmg(mlndabec);
  mlmg.setMaxIter(max_iter);
  mlmg.setMaxFmgIter(max_fmg_iter);
  mlmg.setVerbose(verbose);
  mlmg.setBottomVerbose(bottom_verbose);  

  mlmg.solve({&S_Z}, {&Rhs[2]}, tol_rel, tol_abs);
  
  // We need to compute boundary conditions again after each update
  /*S_EM_output[0].FillBoundary(geom.periodicity());
  S_EM_output[1].FillBoundary(geom.periodicity());
  S_EM_edge_output.FillBoundary(geom.periodicity());
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM_output[0], geom, bc_EM);    
  FillDomainBoundary(S_EM_output[1], geom, bc_EM);
  FillDomainBoundary(S_EM_edge_output, geom, bc_EM);
  */
  MultiFab::Copy(S_EM_X_new, S_EM_output[0], 0, 0, 6, 0);
  FillPatch(*this, S_EM_output[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_new, S_EM_output[1], 0, 0, 6, 0);
  FillPatch(*this, S_EM_output[1], NUM_GROW, time, EM_Y_Type, 0, 6);
  MultiFab::Copy(S_EM_XY_new, S_EM_edge_output, 0, 0, 6, 0);
  FillPatch(*this, S_EM_edge_output, NUM_GROW, time, EM_XY_Type, 0, 6);
#endif

  // Update electric field
  for (MFIter mfi(S_EM_output[0], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arrX = S_EM_output[0].array(mfi);
      const auto& arrXY = S_EM_edge_output.array(mfi);
      const auto& arrX_old = S_EM_input[0].array(mfi);
      const auto& arrXY_old = S_EM_edge_input.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      const Real y = -0.5 + (double(j)+0.5) * dx[1];
	      for(int i = lo.x; i <= hi.x; i++)		
		{
		  const Real x = -0.5 + (double(i)) * dx[0];
		  Real dyBz = (arrXY(i,j+1,k,BZ_LOCAL)-arrXY(i,j,k,BZ_LOCAL))/dx[1];		  
		  Real dyBzOld = (arrXY_old(i,j+1,k,BZ_LOCAL)-arrXY_old(i,j,k,BZ_LOCAL))/dx[1];
		  
		  arrX(i,j,k,EX_LOCAL) = arrX_old(i,j,k,EX_LOCAL) + 0.5*dt*c*c*dyBz + 0.5*dt*c*c*dyBzOld;

		}
	    }
	}
    }
  for (MFIter mfi(S_EM_output[1], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arrY = S_EM_output[1].array(mfi);
      const auto& arrXY = S_EM_edge_output.array(mfi);
      const auto& arrY_old = S_EM_input[1].array(mfi);
      const auto& arrXY_old = S_EM_edge_input.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      const Real y = -0.5 + (double(j)) * dx[1];
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  const Real x = -0.5 + (double(i)+0.5) * dx[0];
		  Real dxBz = (arrXY(i+1,j,k,BZ_LOCAL)-arrXY(i,j,k,BZ_LOCAL))/dx[0];
		  Real dxBzOld = (arrXY_old(i+1,j,k,BZ_LOCAL)-arrXY_old(i,j,k,BZ_LOCAL))/dx[0];
		  
		  arrY(i,j,k,EY_LOCAL) = arrY_old(i,j,k,EY_LOCAL) - 0.5*dt*c*c*dxBz - 0.5*dt*c*c*dxBzOld;
		}
	    }
	}
    }

  // We need to compute boundary conditions again after each update
  /*S_EM_output[0].FillBoundary(geom.periodicity());
  S_EM_output[1].FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM_output[0], geom, bc_EM);    
  FillDomainBoundary(S_EM_output[1], geom, bc_EM);
  */
  for (MFIter mfi(S_output, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arr = S_output.array(mfi);
      const auto& arr_old = S_input.array(mfi);
      const auto& arrX = S_EM_output[0].array(mfi);
      const auto& arrY = S_EM_output[1].array(mfi);
      const auto& arrX_old = S_EM_input[0].array(mfi);
      const auto& arrY_old = S_EM_input[1].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real dxBy = (arrX(i+1,j,k,BY_LOCAL)-arrX(i,j,k,BY_LOCAL))/dx[0];
		  Real dyBx = (arrY(i,j+1,k,BX_LOCAL)-arrY(i,j,k,BX_LOCAL))/dx[1];
		  Real dxByOld = (arrX_old(i+1,j,k,BY_LOCAL)-arrX_old(i,j,k,BY_LOCAL))/dx[0];
		  Real dyBxOld = (arrY_old(i,j+1,k,BX_LOCAL)-arrY_old(i,j,k,BX_LOCAL))/dx[1];

		  arr(i,j,k,EZ) = arr_old(i,j,k,EZ) + 0.5*dt*c*c*dxBy + 0.5*dt*c*c*dxByOld - 0.5*dt*c*c*dyBx - 0.5*dt*c*c*dyBxOld;
		}
	    }
	}
    }

  // We need to compute boundary conditions again after each update
  /*S_EM_output[0].FillBoundary(geom.periodicity());
  S_EM_output[1].FillBoundary(geom.periodicity());
  
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM_output[0], geom, bc_EM);    
  FillDomainBoundary(S_EM_output[1], geom, bc_EM);
  */
  MultiFab::Copy(S_EM_X_new, S_EM_output[0], 0, 0, 6, 0);
  FillPatch(*this, S_EM_output[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2)
  MultiFab::Copy(S_EM_Y_new, S_EM_output[1], 0, 0, 6, 0);
  FillPatch(*this, S_EM_output[1], NUM_GROW, time, EM_Y_Type, 0, 6);
#endif

  // Compute cell-centred EM fields from Yee-grid EM fields
  for (MFIter mfi(S_output, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = S_output.array(mfi);
      const auto& arrEM_X = S_EM_output[0].array(mfi);
      const auto& arrEM_Y = S_EM_output[1].array(mfi);
      const auto& arrEM_XY = S_EM_edge_output.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 
		  
  		  arr(i,j,k,BX) = (arrEM_Y(i,j+1,k,BX_LOCAL)+arrEM_Y(i,j,k,BX_LOCAL))/2.0;
  		  arr(i,j,k,BY) = (arrEM_X(i+1,j,k,BY_LOCAL)+arrEM_X(i,j,k,BY_LOCAL))/2.0;
		  arr(i,j,k,BZ) = (arrEM_XY(i+1,j,k,BZ_LOCAL)+arrEM_XY(i,j+1,k,BZ_LOCAL)+
				   arrEM_XY(i+1,j+1,k,BZ_LOCAL)+arrEM_XY(i,j,k,BZ_LOCAL))/4.0;
  		  arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0;
  		  arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0;

  		}
  	    }
  	}       
    }
    
  /*// We need to compute boundary conditions again after each update
  S_output.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_output, geom, bc);  
  */
  
  MultiFab::Copy(S_new, S_output, BX, BX, 6, 0);
  MultiFab::Copy(S_EM_XY_new, S_EM_edge_output, 0, 0, 6, 0);
}
void CAMReXmp::MaxwellSolverFVTDWENO(Array<MultiFab,AMREX_SPACEDIM>& S_EM_source, MultiFab& S_source, const Real* dx, Real dt)
{

  MultiFab& S_dest = get_new_data(Phi_Type);
  MultiFab& S_EM_X_dest = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y_dest = get_new_data(EM_Y_Type);
  
  Array<MultiFab,AMREX_SPACEDIM> S_EM_dest;
  S_EM_dest[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
  //FillPatch(*this, S_EM_dest[0], NUM_GROW, time, EM_X_Type, 0, 6);
#if (AMREX_SPACEDIM >= 2) 
  S_EM_dest[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
  //FillPatch(*this, S_EM_dest[1], NUM_GROW, time, EM_Y_Type, 0, 6);
#endif

  // Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes[amrex::SpaceDim];
  
  // Define the appropriate size for the flux MultiFab.
  // Fluxes are defined at cell faces - this is taken care of by the
  // surroundingNodes(j) command, ensuring the size of the flux
  for (int j = 0; j < amrex::SpaceDim; j++)
  {
    BoxArray ba = S_source.boxArray();
    ba.surroundingNodes(j);
    fluxes[j].define(ba, dmap, NUM_STATE, 0);
  }

  // Set up a multifab that will contain the fluxes for the electromagnetic fields
  // It will be a nodal multifab (i.e. stored at the corner)
  MultiFab fluxesEM;
  BoxArray ba = S_EM_source[0].boxArray();
  ba.surroundingNodes(1);
  const DistributionMapping& dmX = S_EM_source[0].DistributionMap();
  fluxesEM.define(ba, dmX, 6, 0);

  // Note that fluxesEM and fluxes are different definitions
  // fluxesEM are EM star states, fluxes are fluxes
  
  // Useful when dealing with high-order reconstruction
  // To get first-order slope in y-direction for Bx: slopes[BX_LOCALL].array(i,j,k,Y)
  // Useful to define indices int X,Y,Z in second-order linear reconstruction
  // For third-order quadratic reconstruction use int X,Y,Z,XX,YY,ZZ,XY,YZ,XZ
  // Although for x-components fields, only exist Y,Z,YY,ZZ,YZ (see Balsara et al. 2017)
  
  Array<MultiFab,6> slopes;
  // number of slopes
  // 9 slopes for third-order parabolic reconstruction
  // Ex = Ex^0 + Ex^y*(y/Delta_y) + Ex^z*(z/Delta_z) + Ex^yy*((y/Delta_y)^2-1/12) + ...
  // where Ex^y, Ex^z, Ex^yy are the slopes
  // For 2D code, only 2 slopes for x- and y-compoenents
  // For clarity, will use 9 slopes for third-order, first one represents slope in x-direction
  int nSlopes = 9;
  // indices for the slopes
  // for 2D, do not need Z indices
  int X=0,Y=1,Z=2,XX=3,YY=4,ZZ=5,XY=6,YZ=7,XZ=8;
  slopes[BX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);  
  slopes[EX_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, nSlopes, 1);
#if (AMREX_SPACEDIM >= 2)
  slopes[BY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, nSlopes, 1);
#else  
  // For 1D code it is cell-centred
  slopes[BY_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EY_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
  
#if (AMREX_SPACEDIM != 3)  
  // For 2D code it is cell-centred
  slopes[BZ_LOCAL].define(grids, dmap, nSlopes, 1);
  slopes[EZ_LOCAL].define(grids, dmap, nSlopes, 1);
#endif
  
  // Three components, one for each direction
  MultiFab slopesCharge(grids, dmap, nSlopes, 1);

  // Compute slopes of the charge densities  
  for (MFIter mfi(slopesCharge, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());
      
      // uses old charge
      const Array4<Real> arr = S_source.array(mfi);
      Array4<Real> slopes = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  Vector<Real> data_i = get_data_stencil(arr, i, j, k, 2, 0, 0, RHO_I);
		  Vector<Real> data_e = get_data_stencil(arr, i, j, k, 2, 0, 0, RHO_E);
		  Vector<Real> data_charge = get_charge_scaled(data_i,data_e);
		  std::array<Real, 2> slopesX = WENO3_slope(data_charge);
		  slopes(i,j,k,X) = slopesX[0];
		  slopes(i,j,k,XX) = slopesX[1];

		  data_i = get_data_stencil(arr, i, j, k, 0, 2, 0, RHO_I);
		  data_e = get_data_stencil(arr, i, j, k, 0, 2, 0, RHO_E);
		  data_charge = get_charge_scaled(data_i,data_e);
		  std::array<Real, 2> slopesY = WENO3_slope(data_charge);
		  slopes(i,j,k,Y) = slopesY[0];
		  slopes(i,j,k,YY) = slopesY[1];
		  
		  data_i = get_data_stencil(arr, i, j, k, 1, 1, 0, RHO_I);
		  data_e = get_data_stencil(arr, i, j, k, 1, 1, 0, RHO_E);
		  data_charge = get_charge_scaled(data_i,data_e);
		  Real slopesCross = WENO3_slopeCross(data_charge,
						      {slopes(i,j,k,X),slopes(i,j,k,XX),slopes(i,j,k,Y),slopes(i,j,k,YY)});
		  slopes(i,j,k,XY) = slopesCross;
		  
		  if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
		    {
		      int iOffset = 1, jOffset = 0, kOffset = 0;
		      Real u_iMinus1 = (r_i*arr(i-iOffset,j-jOffset,k-kOffset,RHO_I)
					+ r_e*arr(i-iOffset,j-jOffset,k-kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_iPlus1 = (r_i*arr(i+iOffset,j+jOffset,k+kOffset,RHO_I)
				       + r_e*arr(i+iOffset,j+jOffset,k+kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      slopes(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
		      slopes(i,j,k,XX) = 0.0, slopes(i,j,k,XY) = 0.0; 
		    }
		  if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
		    {
		      int iOffset = 0, jOffset = 1, kOffset = 0;
		      Real u_iMinus1 = (r_i*arr(i-iOffset,j-jOffset,k-kOffset,RHO_I)
					+ r_e*arr(i-iOffset,j-jOffset,k-kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_i = (r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E))/(lambda_d*lambda_d*l_r);
		      Real u_iPlus1 = (r_i*arr(i+iOffset,j+jOffset,k+kOffset,RHO_I)
				       + r_e*arr(i+iOffset,j+jOffset,k+kOffset,RHO_E))/(lambda_d*lambda_d*l_r);
		      slopes(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
		      slopes(i,j,k,YY) = 0.0, slopes(i,j,k,XY) = 0.0; 			  
		    }
		}
	    }
	}      
    }
  
  // Compute cell-centred z-components slopes in x- and y-direction
  //for (int d = 0; d < AMREX_SPACEDIM ; d++)
    {

      //const int iOffset = ( d == 0 ? 1 : 0);
      //const int jOffset = ( d == 1 ? 1 : 0);
      //const int kOffset = ( d == 2 ? 1 : 0);
      
      for (MFIter mfi(slopes[BZ_LOCAL], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const Dim3 hiDomain = ubound(geom.Domain());
	  
	  const Array4<Real> arr = S_source.array(mfi);
	  Array4<Real> slopesBZ = slopes[BZ_LOCAL].array(mfi);
	  Array4<Real> slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y-1; j <= hi.y+1; j++)
		{
		  for(int i = lo.x-1; i <= hi.x+1; i++)
		    {
		      // x slopes
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, BZ);
		      std::array<Real, 2> slopesX = WENO3_slope(dataX);
		      slopesBZ(i,j,k,X) = slopesX[0];
		      slopesBZ(i,j,k,XX) = slopesX[1];

		      dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, EZ);
		      slopesX = WENO3_slope(dataX);
		      slopesEZ(i,j,k,X) = slopesX[0];
		      slopesEZ(i,j,k,XX) = slopesX[1];
		      
		      // y slopes
		      dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, BZ);
		      slopesX = WENO3_slope(dataX);
		      slopesBZ(i,j,k,Y) = slopesX[0];
		      slopesBZ(i,j,k,YY) = slopesX[1];
		     		      
		      dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, EZ);
		      slopesX = WENO3_slope(dataX);
		      slopesEZ(i,j,k,Y) = slopesX[0];
		      slopesEZ(i,j,k,YY) = slopesX[1];
		      
		      // cross slopes
		      Vector<Real> dataXY = get_data_stencil(arr, i, j, k, 1, 1, 0, BZ);
		      Real slopesCross = WENO3_slopeCross(dataXY,
							  {slopesBZ(i,j,k,X),slopesBZ(i,j,k,XX),slopesBZ(i,j,k,Y),slopesBZ(i,j,k,YY)});
		      slopesBZ(i,j,k,XY) = slopesCross;

		      dataXY = get_data_stencil(arr, i, j, k, 1, 1, 0, EZ);
		      slopesCross = WENO3_slopeCross(dataXY,
						     {slopesEZ(i,j,k,X),slopesEZ(i,j,k,XX),slopesEZ(i,j,k,Y),slopesEZ(i,j,k,YY)});
		      slopesEZ(i,j,k,XY) = slopesCross;		      		      

		      if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
			{			      
			  {
			    Real u_iMinus1 = arr(i-1,j,k,BZ);
			    Real u_i = arr(i,j,k,BZ);
			    Real u_iPlus1 = arr(i+1,j,k,BZ);
			    slopesBZ(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesBZ(i,j,k,XX) = 0.0, slopesBZ(i,j,k,XY) = 0.0;
			  }
			  {				  
			    Real u_iMinus1 = arr(i-1,j,k,EZ);
			    Real u_i = arr(i,j,k,EZ);
			    Real u_iPlus1 = arr(i+1,j,k,EZ);
			    slopesEZ(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesEZ(i,j,k,XX) = 0.0, slopesEZ(i,j,k,XY) = 0.0;
			  }			      
			}
		      if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
			{
			  {			
			    Real u_iMinus1 = arr(i,j-1,k,BZ);
			    Real u_i = arr(i,j,k,BZ);
			    Real u_iPlus1 = arr(i,j+1,k,BZ);
			    slopesBZ(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesBZ(i,j,k,YY) = 0.0, slopesBZ(i,j,k,XY) = 0.0;
			    
			  }			  
			  {
			    Real u_iMinus1 = arr(i,j-1,k,EZ);
			    Real u_i = arr(i,j,k,EZ);
			    Real u_iPlus1 = arr(i,j+1,k,EZ);
			    slopesEZ(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			    slopesEZ(i,j,k,YY) = 0.0, slopesEZ(i,j,k,XY) = 0.0;
			  }			      
			}
		    }
		}
	    }      
	}
    }
  
  // Compute y-components slopes in x-direction
  for (MFIter mfi(slopes[BY_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());

      const Array4<Real> arr = S_EM_source[1].array(mfi);
      Array4<Real> slopesBY = slopes[BY_LOCAL].array(mfi);
      Array4<Real> slopesEY = slopes[EY_LOCAL].array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  Vector<Real> dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, BY_LOCAL);
		  std::array<Real, 2> slopesX = WENO3_slope(dataX);
		  slopesBY(i,j,k,X) = slopesX[0];
		  slopesBY(i,j,k,XX) = slopesX[1];

		  dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, EY_LOCAL);
		  slopesX = WENO3_slope(dataX);
		  slopesEY(i,j,k,X) = slopesX[0];
		  slopesEY(i,j,k,XX) = slopesX[1];
		  
		  if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
		    {			      
		      {
			Real u_iMinus1 = arr(i-1,j,k,BY_LOCAL);
			Real u_i = arr(i,j,k,BY_LOCAL);		    
			Real u_iPlus1 = arr(i+1,j,k,BY_LOCAL);
			slopesBY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesBY(i,j,k,XX) = 0.0;
		      }
		      {
			Real u_iMinus1 = arr(i-1,j,k,EY_LOCAL);
			Real u_i = arr(i,j,k,EY_LOCAL);		    
			Real u_iPlus1 = arr(i+1,j,k,EY_LOCAL);
			slopesEY(i,j,k,X) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesEY(i,j,k,XX) = 0.0;
		      }			      
		    }
		}
	    }
	}      
    }
  // Compute x-components slopes in y-direction
  for (MFIter mfi(slopes[BX_LOCAL], true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const Dim3 hiDomain = ubound(geom.Domain());

      const Array4<Real> arr = S_EM_source[0].array(mfi);
      Array4<Real> slopesBX = slopes[BX_LOCAL].array(mfi);
      Array4<Real> slopesEX = slopes[EX_LOCAL].array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  Vector<Real> dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, BX_LOCAL);
		  std::array<Real, 2> slopesX = WENO3_slope(dataX);
		  slopesBX(i,j,k,Y) = slopesX[0];
		  slopesBX(i,j,k,YY) = slopesX[1];

		  dataX = get_data_stencil(arr, i, j, k, 0, 2, 0, EX_LOCAL);
		  slopesX = WENO3_slope(dataX);
		  slopesEX(i,j,k,Y) = slopesX[0];
		  slopesEX(i,j,k,YY) = slopesX[1];
		  
		  if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
		    {
		      {
			Real u_iMinus1 = arr(i,j-1,k,BX_LOCAL);
			Real u_i = arr(i,j,k,BX_LOCAL);		    
			Real u_iPlus1 = arr(i,j+1,k,BX_LOCAL);
			slopesBX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesBX(i,j,k,YY) = 0.0;			      
		      }
		      {
			Real u_iMinus1 = arr(i,j-1,k,EX_LOCAL);
			Real u_i = arr(i,j,k,EX_LOCAL);		    
			Real u_iPlus1 = arr(i,j+1,k,EX_LOCAL);
			slopesEX(i,j,k,Y) = TVD_slope(u_iMinus1,u_i,u_iPlus1);
			slopesEX(i,j,k,YY) = 0.0;			      
		      }			      
		    }
		}
	    }
	}      
    }


  // Coefficients for the magnetic and electric fields
  // For third order, there are 48 coeff.
  // i.e. a0,ax,ay,az,axx,...,b0,bx,...,c0,cx,...,czz
  // use 1 ghost cell
  int a0=0,ax=1,ay=2,az=3,axx=4,ayy=5,azz=6,axy=7,ayz=8,axz=9,axxx=10,axxy=11,axxz=12,axyy=13,axzz=14,axyz=15;
  int b0=16,bx=17,by=18,bz=19,bxx=20,byy=21,bzz=22,bxy=23,byz=24,bxz=25,byyy=26,bxyy=27,byyz=28,bxxy=29,byzz=30,bxyz=31;
  int c0=32,cx=33,cy=34,cz=35,cxx=36,cyy=37,czz=38,cxy=39,cyz=40,cxz=41,czzz=42,cxzz=43,cyzz=44,cxxz=45,cyyz=46,cxyz=47;
  MultiFab Bcoeff(grids, dmap, 48, 1);
  MultiFab Ecoeff(grids, dmap, 48, 1);

  // Compute coefficients for the field function
  for (MFIter mfi(Bcoeff, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      // coeff
      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      // data
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi);       

      // slopes
      const auto& slopesBX = slopes[BX_LOCAL].array(mfi);
      const auto& slopesBY = slopes[BY_LOCAL].array(mfi);
      const auto& slopesBZ = slopes[BZ_LOCAL].array(mfi);
      const auto& slopesEX = slopes[EX_LOCAL].array(mfi);
      const auto& slopesEY = slopes[EY_LOCAL].array(mfi);
      const auto& slopesEZ = slopes[EZ_LOCAL].array(mfi);
      
      const auto& slopesQ = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y-1; j <= hi.y+1; j++)
	    {
	      for(int i = lo.x-1; i <= hi.x+1; i++)
		{
		  // third order terms
		  Bc(i,j,k,ayy) = (slopesBX(i+1,j,k,YY)+slopesBX(i,j,k,YY))/2.0;
		  Bc(i,j,k,axyy) = slopesBX(i+1,j,k,YY)-slopesBX(i,j,k,YY);
		  Bc(i,j,k,azz) = 0.0;
		  Bc(i,j,k,axzz) = 0.0;
		  Bc(i,j,k,ayz) = 0.0;
		  Bc(i,j,k,axyz) = 0.0;
		  Bc(i,j,k,bxx) = (slopesBY(i,j+1,k,XX)+slopesBY(i,j,k,XX))/2.0;
		  Bc(i,j,k,bxxy) = slopesBY(i,j+1,k,XX)-slopesBY(i,j,k,XX);
		  Bc(i,j,k,bzz) = 0.0;
		  Bc(i,j,k,byzz) = 0.0;
		  Bc(i,j,k,bxz) = 0.0;
		  Bc(i,j,k,bxyz) = 0.0;
		  Bc(i,j,k,cxx) = slopesBZ(i,j,k,XX);
		  Bc(i,j,k,cxxz) = 0.0;
		  Bc(i,j,k,cyy) = slopesBZ(i,j,k,YY);
		  Bc(i,j,k,cyyz) = 0.0;
		  Bc(i,j,k,cxy) = slopesBZ(i,j,k,XY);
		  Bc(i,j,k,cxyz) = 0.0;
		  Bc(i,j,k,axxx) = -dx[0]/3.0*(Bc(i,j,k,bxxy)/dx[1]);
		  Bc(i,j,k,byyy) = -dx[1]/3.0*(Bc(i,j,k,axyy)/dx[0]);
		  Bc(i,j,k,czzz) = 0.0;
		  Bc(i,j,k,axxy) = 0.0;
		  Bc(i,j,k,bxyy) = 0.0;
		  Bc(i,j,k,byyz) = 0.0;
		  Bc(i,j,k,cyzz) = 0.0;
		  Bc(i,j,k,cxzz) = 0.0;
		  Bc(i,j,k,axxz) = 0.0;
		  // second order terms
		  Bc(i,j,k,ay) = (slopesBX(i+1,j,k,Y)+slopesBX(i,j,k,Y))/2.0 - Bc(i,j,k,axxy)/6.0;
		  Bc(i,j,k,axy) = (slopesBX(i+1,j,k,Y)-slopesBX(i,j,k,Y));
		  Bc(i,j,k,az) = 0.0;
		  Bc(i,j,k,axz) = 0.0;
		  Bc(i,j,k,bx) = (slopesBY(i,j+1,k,X)+slopesBY(i,j,k,X))/2.0 - Bc(i,j,k,bxyy)/6.0;
		  Bc(i,j,k,bxy) = (slopesBY(i,j+1,k,X)-slopesBY(i,j,k,X));
		  Bc(i,j,k,bz) = 0.0;
		  Bc(i,j,k,byz) = 0.0;
		  Bc(i,j,k,cx) = slopesBZ(i,j,k,X) - Bc(i,j,k,cxzz)/6.0;
		  Bc(i,j,k,cxz) = 0.0;
		  Bc(i,j,k,cy) = slopesBZ(i,j,k,Y) - Bc(i,j,k,cyzz)/6.0;
		  Bc(i,j,k,cyz) = 0.0;
		  Bc(i,j,k,axx) = -dx[0]/2.0*(Bc(i,j,k,bxy)/dx[1]); // +cxz/dx[2]
		  Bc(i,j,k,byy) = -dx[1]/2.0*(Bc(i,j,k,axy)/dx[0]); // ++cyz/dx[2]
		  Bc(i,j,k,czz) = 0.0;
		  Bc(i,j,k,a0) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0 - Bc(i,j,k,axx)/6.0;
		  Bc(i,j,k,ax) = (arrEM_X(i+1,j,k,BX_LOCAL)-arrEM_X(i,j,k,BX_LOCAL)) - Bc(i,j,k,axxx)/10.0;
		  Bc(i,j,k,b0) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0 - Bc(i,j,k,byy)/6.0;
		  Bc(i,j,k,by) = (arrEM_Y(i,j+1,k,BY_LOCAL)-arrEM_Y(i,j,k,BY_LOCAL)) - Bc(i,j,k,byyy)/10.0;
		  Bc(i,j,k,c0) = arr(i,j,k,BZ) - Bc(i,j,k,czz)/6.0;
		  Bc(i,j,k,cz) = 0.0;

		  // third order terms
		  Ec(i,j,k,ayy) = (slopesEX(i+1,j,k,YY)+slopesEX(i,j,k,YY))/2.0;
		  Ec(i,j,k,axyy) = slopesEX(i+1,j,k,YY)-slopesEX(i,j,k,YY);
		  Ec(i,j,k,azz) = 0.0;
		  Ec(i,j,k,axzz) = 0.0;
		  Ec(i,j,k,ayz) = 0.0;
		  Ec(i,j,k,axyz) = 0.0;
		  Ec(i,j,k,bxx) = (slopesEY(i,j+1,k,XX)+slopesEY(i,j,k,XX))/2.0;
		  Ec(i,j,k,bxxy) = slopesEY(i,j+1,k,XX)-slopesEY(i,j,k,XX);
		  Ec(i,j,k,bzz) = 0.0;
		  Ec(i,j,k,byzz) = 0.0;
		  Ec(i,j,k,bxz) = 0.0;
		  Ec(i,j,k,bxyz) = 0.0;
		  Ec(i,j,k,cxx) = slopesEZ(i,j,k,XX);
		  Ec(i,j,k,cxxz) = 0.0;
		  Ec(i,j,k,cyy) = slopesEZ(i,j,k,YY);
		  Ec(i,j,k,cyyz) = 0.0;
		  Ec(i,j,k,cxy) = slopesEZ(i,j,k,XY);
		  Ec(i,j,k,cxyz) = 0.0;
		  Ec(i,j,k,axxx) = -dx[0]/3.0*(Ec(i,j,k,bxxy)/dx[1]-slopesQ(i,j,k,XX));
		  Ec(i,j,k,byyy) = -dx[1]/3.0*(Ec(i,j,k,axyy)/dx[0]-slopesQ(i,j,k,YY));
		  Ec(i,j,k,czzz) = 0.0;
		  Ec(i,j,k,axxy) = dx[0]*dx[1]/(2.0*(dx[0]+dx[1]))*slopesQ(i,j,k,XY);
		  Ec(i,j,k,bxyy) = Ec(i,j,k,axxy);
		  Ec(i,j,k,byyz) = 0.0;
		  Ec(i,j,k,cyzz) = 0.0;
		  Ec(i,j,k,cxzz) = 0.0;
		  Ec(i,j,k,axxz) = 0.0;
		  // second order terms
		  Ec(i,j,k,ay) = (slopesEX(i+1,j,k,Y)+slopesEX(i,j,k,Y))/2.0 - Ec(i,j,k,axxy)/6.0;
		  Ec(i,j,k,axy) = (slopesEX(i+1,j,k,Y)-slopesEX(i,j,k,Y));
		  Ec(i,j,k,az) = 0.0;
		  Ec(i,j,k,axz) = 0.0;
		  Ec(i,j,k,bx) = (slopesEY(i,j+1,k,X)+slopesEY(i,j,k,X))/2.0 - Ec(i,j,k,bxyy)/6.0;
		  Ec(i,j,k,bxy) = (slopesEY(i,j+1,k,X)-slopesEY(i,j,k,X));
		  Ec(i,j,k,bz) = 0.0;
		  Ec(i,j,k,byz) = 0.0;
		  Ec(i,j,k,cx) = slopesEZ(i,j,k,X) - Ec(i,j,k,cxzz)/6.0;
		  Ec(i,j,k,cxz) = 0.0;
		  Ec(i,j,k,cy) = slopesEZ(i,j,k,Y) - Ec(i,j,k,cyzz)/6.0;
		  Ec(i,j,k,cyz) = 0.0;
		  Ec(i,j,k,axx) = -dx[0]/2.0*(Ec(i,j,k,bxy)/dx[1]-slopesQ(i,j,k,X)); // +cxz/dx[2]
		  Ec(i,j,k,byy) = -dx[1]/2.0*(Ec(i,j,k,axy)/dx[0]-slopesQ(i,j,k,Y)); // ++cyz/dx[2]
		  Ec(i,j,k,czz) = 0.0;
		  Ec(i,j,k,a0) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0 - Ec(i,j,k,axx)/6.0;
		  Ec(i,j,k,ax) = (arrEM_X(i+1,j,k,EX_LOCAL)-arrEM_X(i,j,k,EX_LOCAL)) - Ec(i,j,k,axxx)/10.0;
		  Ec(i,j,k,b0) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0 - Ec(i,j,k,byy)/6.0;
		  Ec(i,j,k,by) = (arrEM_Y(i,j+1,k,EY_LOCAL)-arrEM_Y(i,j,k,EY_LOCAL)) - Ec(i,j,k,byyy)/10.0;
		  Ec(i,j,k,c0) = arr(i,j,k,EZ) - Ec(i,j,k,czz)/6.0; //czz*dx[2]*dx[2]/6.0
		  Ec(i,j,k,cz) = 0.0;		  

		}
	    }
	}
    }

  // Compute edge components of the EM fields    
  for (MFIter mfi(fluxesEM, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      const Dim3 hiDomain = ubound(geom.Domain());
      
      // use old array, S_source should also work
      //const auto& arr = S_source.array(mfi);
      const auto& arr = S_source.array(mfi);
      const auto& arrEM_X = S_EM_source[0].array(mfi);
      const auto& arrEM_Y = S_EM_source[1].array(mfi); 
      const auto& fluxArrEM = fluxesEM.array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{		 
		  
		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
		  // LD state
		  x = dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = -1, jOffset = -1;
		  Vector<Real> EM_LD = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  		  
		  // RD state
		  x = -dx[0]/2.0, y = dx[1]/2.0;
		  iOffset = 0, jOffset = -1;
		  Vector<Real> EM_RD = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  // LU state
		  x = dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = -1, jOffset = 0;
		  Vector<Real> EM_LU = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // RU state
		  x = -dx[0]/2.0, y = -dx[1]/2.0;
		  iOffset = 0, jOffset = 0;
		  Vector<Real> EM_RU = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
		  // If the reconstruction is consistent, then EyRU=EyRD
		  Real EyR = (EM_RU[EY_LOCAL]+EM_RD[EY_LOCAL])/2.0;
		  Real EyL = (EM_LU[EY_LOCAL]+EM_LD[EY_LOCAL])/2.0;
		  Real ExU = (EM_RU[EX_LOCAL]+EM_LU[EX_LOCAL])/2.0;
		  Real ExD = (EM_RD[EX_LOCAL]+EM_LD[EX_LOCAL])/2.0;
		  fluxArrEM(i,j,k,BZ_LOCAL) = 0.25*(EM_RU[BZ_LOCAL]+EM_RD[BZ_LOCAL]+EM_LU[BZ_LOCAL]+EM_LD[BZ_LOCAL])
		    - 0.5/c*(EyR-EyL) + 0.5/c*(ExU-ExD);

		  Real ByR = (EM_RU[BY_LOCAL]+EM_RD[BY_LOCAL])/2.0;
		  Real ByL = (EM_LU[BY_LOCAL]+EM_LD[BY_LOCAL])/2.0;
		  Real BxU = (EM_RU[BX_LOCAL]+EM_LU[BX_LOCAL])/2.0;
		  Real BxD = (EM_RD[BX_LOCAL]+EM_LD[BX_LOCAL])/2.0;

		  fluxArrEM(i,j,k,EZ_LOCAL) = 0.25*(EM_RU[EZ_LOCAL]+EM_RD[EZ_LOCAL]+EM_LU[EZ_LOCAL]+EM_LD[EZ_LOCAL])
		    + 0.5*c*(ByR-ByL) - 0.5*c*(BxU-BxD);
		  		  
		}
	    }
	}       
    }

  // Update face-centred EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
    {

      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      int d_EM = (d==0) ? 1 : 0;
    
      // Loop over all the patches at this level
      for (MFIter mfi(S_EM_dest[d_EM], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();

	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  // Indexable arrays for the data, and the directional flux
	  // Based on the corner-centred definition of the flux array, the
	  // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	  //const auto& arr = S_source.array(mfi);
	  const auto& arrEM = S_EM_dest[d_EM].array(mfi);
	  const auto& arrEMOld = S_EM_source[d_EM].array(mfi);	  
	  const auto& arr = S_source.array(mfi);
	  const auto& fluxArrEM = fluxesEM.array(mfi);
	  const auto& fluxArr = fluxes[d_EM].array(mfi);

	  const Dim3 hiDomain = ubound(geom.Domain());
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      // 2D code; z-component updated using cell-centred scheme
		      arrEM(i,j,k,BX_LOCAL+d_EM) = arrEMOld(i,j,k,BX_LOCAL+d_EM) + std::pow(-1,d)*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, EZ_LOCAL) - fluxArrEM(i,j,k,EZ_LOCAL));
		      arrEM(i,j,k,EX_LOCAL+d_EM) = arrEMOld(i,j,k,EX_LOCAL+d_EM) - std::pow(-1,d)*c*c*(dt / dx[d]) * (fluxArrEM(i+iOffset, j+jOffset, k+kOffset, BZ_LOCAL) - fluxArrEM(i,j,k,BZ_LOCAL));
		      
		    }
		}
	    }      
	}          
    }

  // We need to compute boundary conditions again after each update
  /*S_EM_dest[0].FillBoundary(geom.periodicity());
  S_EM_dest[1].FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_EM_dest[0], geom, bc_EM);    
  FillDomainBoundary(S_EM_dest[1], geom, bc_EM);*/

  MultiFab::Copy(S_EM_X_dest, S_EM_dest[0], 0, 0, 6, 0);
  MultiFab::Copy(S_EM_Y_dest, S_EM_dest[1], 0, 0, 6, 0);
  
  // Compute cell-centred EM fields from face-centred
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the corner-centred definition of the flux array, the
      // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrEM_X = S_EM_dest[0].array(mfi);
      const auto& arrEM_Y = S_EM_dest[1].array(mfi); 

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 
  		  arr(i,j,k,BX) = (arrEM_X(i+1,j,k,BX_LOCAL)+arrEM_X(i,j,k,BX_LOCAL))/2.0;// - Bc(i,j,k,axx)/6.0;
  		  arr(i,j,k,BY) = (arrEM_Y(i,j+1,k,BY_LOCAL)+arrEM_Y(i,j,k,BY_LOCAL))/2.0;// - Bc(i,j,k,byy)/6.0;
  		  arr(i,j,k,EX) = (arrEM_X(i+1,j,k,EX_LOCAL)+arrEM_X(i,j,k,EX_LOCAL))/2.0;// - Ec(i,j,k,axx)/6.0;
  		  arr(i,j,k,EY) = (arrEM_Y(i,j+1,k,EY_LOCAL)+arrEM_Y(i,j,k,EY_LOCAL))/2.0;// - Ec(i,j,k,byy)/6.0;
		  
  		}
  	    }
  	}       
    }

  // Only 2D
  // Compute edge components for y-components at x-faces
  for (MFIter mfi(fluxes[0], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      
      const auto& fluxArrX = fluxes[0].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // 2 gaussian points
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // L state
		  iOffset = -1, jOffset = 0;
  		  x = dx[0]/2.0, y = dx[1]/(2.0*std::sqrt(3.0));		  
		  Vector<Real> EM_L1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = dx[0]/2.0, y = -dx[1]/(2.0*std::sqrt(3.0));		  
		  Vector<Real> EM_L2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);		  
		   
  		  // R state
		  iOffset = 0, jOffset = 0;
  		  x = -dx[0]/2.0, y = dx[1]/(2.0*std::sqrt(3.0));		 
		  Vector<Real> EM_R1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = -dx[0]/2.0, y = -dx[1]/(2.0*std::sqrt(3.0));		 
		  Vector<Real> EM_R2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  Vector<Real> EM_L, EM_R;
		  for (int n=0; n<6; n++)
		    {
		      EM_L.push_back(0.5*(EM_L1[n]+EM_L2[n]));
		      EM_R.push_back(0.5*(EM_R1[n]+EM_R2[n]));
		    }
		  
		  // Note that this is not the flux, but the star states
  		  fluxArrX(i,j,k,BY) = 0.5*(EM_R[BY_LOCAL]+EM_L[BY_LOCAL])
		    + 0.5/c*(EM_R[EZ_LOCAL]-EM_L[EZ_LOCAL]);
  		  fluxArrX(i,j,k,EY) = 0.5*(EM_R[EY_LOCAL]+EM_L[EY_LOCAL])
		    - 0.5*c*(EM_R[BZ_LOCAL]-EM_L[BZ_LOCAL]);
  		}
  	    }
  	}
    }
  
  // Only 2D
  // Compute edge components for x-components at y-faces
  for (MFIter mfi(fluxes[1], true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();
      
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);

      const auto& fluxArrY = fluxes[1].array(mfi);

      const auto& Bc = Bcoeff.array(mfi);
      const auto& Ec = Ecoeff.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{		 		  

		  Real x,y,z;
		  int iOffset,jOffset,kOffset;

		  // 2 gaussian points
		  // For 2D code
		  z = 0.0;
		  kOffset = 0;
		  
  		  // D state
  		  iOffset = 0, jOffset = -1;
		  x = dx[0]/(2.0*std::sqrt(3.0)), y = dx[1]/2.0;		  
		  Vector<Real> EM_D1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = -dx[0]/(2.0*std::sqrt(3.0)), y = dx[1]/2.0;		  
		  Vector<Real> EM_D2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  
  		  // U state
  		  iOffset = 0, jOffset = 0;
		  x = dx[0]/(2.0*std::sqrt(3.0)), y = -dx[1]/2.0;		  
		  Vector<Real> EM_U1 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);
		  x = -dx[0]/(2.0*std::sqrt(3.0)), y = -dx[1]/2.0;		  
		  Vector<Real> EM_U2 = EM_quadraticFunc(Bc,Ec,i+iOffset,j+jOffset,k+kOffset,x,y,z,dx);

		  Vector<Real> EM_D, EM_U;
		  for (int n=0; n<6; n++)
		    {
		      EM_D.push_back(0.5*(EM_D1[n]+EM_D2[n]));
		      EM_U.push_back(0.5*(EM_U1[n]+EM_U2[n]));
		    }
		  
		  // Note that this is not the flux, but the star states 	  
  		  fluxArrY(i,j,k,BX) = 0.5*(EM_U[BX_LOCAL]+EM_D[BX_LOCAL])
		    - 0.5/c*(EM_U[EZ_LOCAL]-EM_D[EZ_LOCAL]);
  		  fluxArrY(i,j,k,EX) = 0.5*(EM_U[EX_LOCAL]+EM_D[EX_LOCAL])
		    + 0.5*c*(EM_U[BZ_LOCAL]-EM_D[BZ_LOCAL]); 
  		}
  	    }
  	}
    }

  // Update cell-centred z-components of EM fields
  for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrOld = S_source.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);

      Array4<Real> slopes = slopesCharge.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
  	{
  	  for(int j = lo.y; j <= hi.y; j++)
  	    {
  	      for(int i = lo.x; i <= hi.x; i++)
  		{
  		  // Update cell-centred z-components becuause it is 2D code
  		  arr(i,j,k,BZ) = arrOld(i,j,k,BZ) - (dt / dx[0]) * (fluxArrX(i+1,j,k,EY) - fluxArrX(i,j,k,EY)) + (dt / dx[1]) * (fluxArrY(i,j+1,k,EX) - fluxArrY(i,j,k,EX));
  		  arr(i,j,k,EZ) = arrOld(i,j,k,EZ) + c*c*(dt / dx[0]) * (fluxArrX(i+1,j,k,BY) - fluxArrX(i,j,k,BY)) - c*c*(dt / dx[1]) * (fluxArrY(i,j+1,k,BX) - fluxArrY(i,j,k,BX));		  

  		}
  	    }
  	}
    }
  // We need to compute boundary conditions again after each update
  /*S_dest.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_dest, geom, bc);  */

}
