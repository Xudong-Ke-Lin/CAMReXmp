#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

// Implicit Maxwell Solver
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>

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

		      if (d_EM==1 && i==256 && (j>252))
			std::cout << arrEM(i,j,k,EX_LOCAL+d_EM) << " " << currentFace << std::endl;
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
void CAMReXmp::hyperbolicMaxwellSolverDivFreeSubCycle(Array<MultiFab,AMREX_SPACEDIM>& S_EM, MultiFab& fluxesEM, MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{

  Real dt_final = dt;
  Real dt_finalSource = dt;
  // assumes a constant cfl=0.2 for the Magnetic reconnection test
  dt = 0.2*std::min(dx[0],dx[1])/c;
  //dt = cfl*std::min(dx[0],dx[1])/c; 
  Real dt_current = dt;
  Real dt_currentSource = dt;

  // supports subcycling
  do{
    amrex::Print() << "Updating hyperbolic part, dt=" << dt << ", current dt=" << dt_current << " and final dt=" << dt_final << std::endl;

    // Compute transverse face-centred EM values from cell-centred  
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
	    // Based on the vertex-centred definition of the flux array, the
	    // data array runs from e.g. [0,N] and the flux array from [0,N+1]
	    const auto& arr = Sborder.array(mfi);
	    const auto& arrEM = S_EM[d].array(mfi);

	    for(int k = lo.z; k <= hi.z; k++)
	      {
		for(int j = lo.y; j <= hi.y; j++)
		  {
		    for(int i = lo.x; i <= hi.x; i++)
		      {
	    
			// Calculate transverse facial components from cell-centred components
			// transverse values, e.g. By, Bz, Ey and Ez for the x-face
			Vector<Real> transverseEM = Maxwell_transverse_comp(arr, i, j, k, iOffset, jOffset, kOffset);	    
			arrEM(i,j,k,BX_LOCAL+(1+d)%3) = transverseEM[BX_LOCAL+(1+d)%3];
			arrEM(i,j,k,BX_LOCAL+(2+d)%3) = transverseEM[BX_LOCAL+(2+d)%3];
			arrEM(i,j,k,EX_LOCAL+(1+d)%3) = transverseEM[EX_LOCAL+(1+d)%3];
			arrEM(i,j,k,EX_LOCAL+(2+d)%3) = transverseEM[EX_LOCAL+(2+d)%3];	   

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
			//if (d_EM==1 && (i==253 && j==253))
			//std::cout << currentFace << std::endl;


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

    if (dt_current+dt > dt_final)
      {
	// if at final time step, stop
	if (dt_current==dt_final)
	  break;
	// at the final time step, use the remaining dt
	dt = std::abs(dt_final - dt_current);      
      }
    dt_current += dt;
  }  while (dt_current <= dt_final);


  return;

  // if commented out, only one source term step
  //dt = cfl*std::min(dx[0],dx[1])/c;
  dt = dt_final;
  dt_currentSource = dt_final;
  
  do{

    amrex::Print() << "Updating source part, dt=" << dt << ", current dt=" << dt_currentSource << " and final dt=" << dt_finalSource << std::endl;

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
	    
			// source terms
			Real currentFace = r_i*fluxArr(i,j,k,RHO_I) + r_e*fluxArr(i,j,k,RHO_E);
			arrEM(i,j,k,EX_LOCAL+d_EM) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;

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

    if (dt_currentSource+dt > dt_finalSource)
      {
	// if at final time step, stop
	if (dt_currentSource==dt_finalSource)
	  break;
	// at the final time step, use the remaining dt
	dt = std::abs(dt_finalSource - dt_currentSource);      
      }
    dt_currentSource += dt;
  }  while (dt_currentSource <= dt_finalSource);
}
void CAMReXmp::hyperbolicMaxwellSolverDivFreeSubCycleIneff(Array<MultiFab,AMREX_SPACEDIM>& S_EM, MultiFab& fluxesEM, MultiFab& S0, MultiFab (&fluxes0)[AMREX_SPACEDIM], MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{

  Real dt_final = dt;
  dt = cfl*std::min(dx[0],dx[1])/c;
  Real dt_current = dt;

  // supports subcycling
  do{

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
	    // Based on the vertex-centred definition of the flux array, the
	    // data array runs from e.g. [0,N] and the flux array from [0,N+1]
	    const auto& arr = Sborder.array(mfi);
	    const auto& arrEM = S_EM[d].array(mfi);

	    for(int k = lo.z; k <= hi.z; k++)
	      {
		for(int j = lo.y; j <= hi.y; j++)
		  {
		    for(int i = lo.x; i <= hi.x; i++)
		      {
	    
			// Calculate transverse facial components from cell-centred components
			// transverse values, e.g. By, Bz, Ey and Ez for the x-face
			Vector<Real> transverseEM = Maxwell_transverse_comp(arr, i, j, k, iOffset, jOffset, kOffset);	    
			arrEM(i,j,k,BX_LOCAL+(1+d)%3) = transverseEM[BX_LOCAL+(1+d)%3];
			arrEM(i,j,k,BX_LOCAL+(2+d)%3) = transverseEM[BX_LOCAL+(2+d)%3];
			arrEM(i,j,k,EX_LOCAL+(1+d)%3) = transverseEM[EX_LOCAL+(1+d)%3];
			arrEM(i,j,k,EX_LOCAL+(2+d)%3) = transverseEM[EX_LOCAL+(2+d)%3];	   

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
			Vector<Real> flux = MUSCL_Hancock_HLLC_flux(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);

			for(int n=0; n<NUM_STATE_FLUID; n++)
			  {		
			    fluxArr(i,j,k,n) = flux[n];		
			  }

			// Fluxes for the z-components of the EM fields because it is 2D code
			Vector<Real> fluxEM = Maxwell_flux_Godunov(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);
			fluxArr(i,j,k,BZ) = fluxEM[BZ_LOCAL];
			fluxArr(i,j,k,EZ) = fluxEM[EZ_LOCAL];

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
			// Update fluid variables
			for(int n=0; n<NUM_STATE_FLUID; n++)
			  {
			    // Conservative update formula
			    arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
			  }
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
    
    // Compute cell-centred Ez source terms
    for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
      {
	const Box& bx = mfi.tilebox();
      
	const Dim3 lo = lbound(bx);
	const Dim3 hi = ubound(bx);

	const Dim3 hiDomain = ubound(geom.Domain());

	Array4<Real> arr = Sborder.array(mfi);
	Array4<Real> arr0 = S0.array(mfi);
	Array4<Real> arrEMX = S_EM[0].array(mfi);
	Array4<Real> arrEMY = S_EM[1].array(mfi);

	for(int k = lo.z; k <= hi.z; k++)
	  {
	    for(int j = lo.y; j <= hi.y; j++)
	      {
		for(int i = lo.x; i <= hi.x; i++)
		  {		    
		    if (dt_current!=dt_final)
		      arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
		    else
		      {
			arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr0(i,j,k,MOMZ_I) + r_e*arr0(i,j,k,MOMZ_E));
		      }
		    
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
	    const auto& fluxArr0 = fluxes0[d_EM].array(mfi);

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

	    
			// source terms
			Real currentFace = r_i*fluxArr(i,j,k,RHO_I) + r_e*fluxArr(i,j,k,RHO_E);
			if (dt_current!=dt_final)
			  arrEM(i,j,k,EX_LOCAL+d_EM) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;
			else
			  {
			    amrex::Print() << currentFace << " " << r_i*fluxArr0(i,j,k,RHO_I) + r_e*fluxArr0(i,j,k,RHO_E) << std::endl;
			     currentFace = r_i*fluxArr0(i,j,k,RHO_I) + r_e*fluxArr0(i,j,k,RHO_E);
			     arrEM(i,j,k,EX_LOCAL+d_EM) -= dt*1.0/(lambda_d*lambda_d*l_r)*currentFace;
			  }
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
    
    // We need to compute boundary conditions again after each update
    S_EM[0].FillBoundary(geom.periodicity());
    S_EM[1].FillBoundary(geom.periodicity());
     
    // added by 2020D 
    // Fill non-periodic physical boundaries                          
    FillDomainBoundary(S_EM[0], geom, bc_EM);    
    FillDomainBoundary(S_EM[1], geom, bc_EM);  

    // Compute cell-centred EM fields from face-centred
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
		    arr(i,j,k,BX) = 0.5*(arrEMX(i,j,k,BX_LOCAL)+arrEMX(i+1,j,k,BX_LOCAL));
		    arr(i,j,k,BY) = 0.5*(arrEMY(i,j,k,BY_LOCAL)+arrEMY(i,j+1,k,BY_LOCAL));
		    arr(i,j,k,EX) = 0.5*(arrEMX(i,j,k,EX_LOCAL)+arrEMX(i+1,j,k,EX_LOCAL));
		    arr(i,j,k,EY) = 0.5*(arrEMY(i,j,k,EY_LOCAL)+arrEMY(i,j+1,k,EY_LOCAL));
	    
		    /*
		    if (i==1 && bc[MOMX_I].lo(0) == BCType::reflect_odd)
		      arr(0,j,k,EX_LOCAL) = arr(1,j,k,EX_LOCAL);
		    if (i==hiDomain.x && bc[MOMX_I].hi(0) == BCType::reflect_odd)
		      arr(hiDomain.x,j,k,EX_LOCAL) = arr(hiDomain.x-1,j,k,EX_LOCAL);
		    if (j==1 && bc[MOMY_I].lo(1) == BCType::reflect_odd)
		      arr(i,0,k,EY_LOCAL) = arr(i,1,k,EY_LOCAL);
		    if (j==hiDomain.y && bc[MOMY_I].hi(1) == BCType::reflect_odd)
		      arr(i,hiDomain.y,k,EY_LOCAL) = arr(i,hiDomain.y-1,k,EY_LOCAL);
		    */
		    // Ez source terms
		    //arr(i,j,k,EZ) = arr(i,j,k,EZ) - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,MOMZ_I) + r_e*arr(i,j,k,MOMZ_E));
		  }
	      }
	  }      
      }
    // We need to compute boundary conditions again after each update       
    Sborder.FillBoundary(geom.periodicity());
	    
    // Fill non-periodic physical boundaries          
    FillDomainBoundary(Sborder, geom, bc);
  
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

    
    if (dt_current+dt > dt_final)
      {
	// if at final time step, stop
	if (dt_current==dt_final)
	  break;
	// at the final time step, use the remaining dt
	dt = std::abs(dt_final - dt_current);      
      }
    dt_current += dt;
  }  while (dt_current <= dt_final);

  MultiFab::Copy(S0, Sborder, BX, BX, 6, NUM_GROW);
  //MultiFab::Copy(S0, Sborder, EZ, EZ, 1, NUM_GROW);
  
}
