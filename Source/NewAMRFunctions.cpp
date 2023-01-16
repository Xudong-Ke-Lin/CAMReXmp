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

// Eigen library
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

// 9x9 identity matrix
MatrixXd identity = MatrixXd::Identity(9, 9);

void CAMReXmp::sourceUpdate(MultiFab& Sborder, const Real* dx, Real dt)
{
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
		  (this->*sourceUpdateWithChosenMethod)(arr,i,j,k,dt);		  
		}
	    }
	}      
    }
  // We need to compute boundary conditions again after each update                           
  Sborder.FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries                           
  FillDomainBoundary(Sborder, geom, bc);
  
}

void CAMReXmp::sourceUpdateEX(Array4<Real>& arr, int i, int j, int k, Real dt)
{
  
  // define conserved variables                     
  Real rho_i = arr(i,j,k,0);
  Real momX_i = arr(i,j,k,1);
  Real momY_i = arr(i,j,k,2);
  Real momZ_i = arr(i,j,k,MOMZ_I);
  Real E_i = arr(i,j,k,ENER_I);
  Real rho_e = arr(i,j,k,RHO_E);
  Real momX_e = arr(i,j,k,MOMX_E);
  Real momY_e = arr(i,j,k,MOMY_E);
  Real momZ_e = arr(i,j,k,MOMZ_E);
  Real E_e = arr(i,j,k,ENER_E);
  Real B_x = arr(i,j,k,BX);
  Real B_y = arr(i,j,k,BY);
  Real B_z = arr(i,j,k,BZ);
  Real E_x = arr(i,j,k,EX);
  Real E_y = arr(i,j,k,EY);
  Real E_z = arr(i,j,k,EZ);
		  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});
		  
  arr(i,j,k,0) = rho_i;
  arr(i,j,k,1) = momX_i + dt*r_i*rho_i/l_r*(E_x + B_z*v_y_i - B_y*v_z_i);
  arr(i,j,k,2) = momY_i + dt*r_i*rho_i/l_r*(E_y + B_x*v_z_i - B_z*v_x_i);
  arr(i,j,k,3) = momZ_i + dt*r_i*rho_i/l_r*(E_z + B_y*v_x_i - B_x*v_y_i);
  arr(i,j,k,ENER_I) = E_i + dt*r_i*rho_i/l_r*(E_x*v_x_i + E_y*v_y_i + E_z*v_z_i);
  arr(i,j,k,RHO_E) = rho_e;
  arr(i,j,k,MOMX_E) = momX_e + dt*r_e*rho_e/l_r*(E_x + B_z*v_y_e - B_y*v_z_e);
  arr(i,j,k,MOMY_E) = momY_e + dt*r_e*rho_e/l_r*(E_y + B_x*v_z_e - B_z*v_x_e);
  arr(i,j,k,MOMZ_E) = momZ_e + dt*r_e*rho_e/l_r*(E_z + B_y*v_x_e - B_x*v_y_e);
  arr(i,j,k,ENER_E) = E_e + dt*r_e*rho_e/l_r*(E_x*v_x_e + E_y*v_y_e + E_z*v_z_e);
  
  if (MaxwellMethod=="HYP")
    {
      arr(i,j,k,EX) = E_x - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_x_i + r_e*rho_e*v_x_e);
      arr(i,j,k,EY) = E_y - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_y_i + r_e*rho_e*v_y_e);
      arr(i,j,k,EZ) = E_z - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_z_i + r_e*rho_e*v_z_e);
    }
}
void CAMReXmp::sourceUpdateIM(Array4<Real>& arr, int i, int j, int k, Real dt)
{
  
  // define conserved variables                     
  Real rho_i = arr(i,j,k,0);
  Real momX_i = arr(i,j,k,1);
  Real momY_i = arr(i,j,k,2);
  Real momZ_i = arr(i,j,k,MOMZ_I);
  Real E_i = arr(i,j,k,ENER_I);
  Real rho_e = arr(i,j,k,RHO_E);
  Real momX_e = arr(i,j,k,MOMX_E);
  Real momY_e = arr(i,j,k,MOMY_E);
  Real momZ_e = arr(i,j,k,MOMZ_E);
  Real E_e = arr(i,j,k,ENER_E);
  Real B_x = arr(i,j,k,BX);
  Real B_y = arr(i,j,k,BY);
  Real B_z = arr(i,j,k,BZ);
  Real E_x = arr(i,j,k,EX);
  Real E_y = arr(i,j,k,EY);
  Real E_z = arr(i,j,k,EZ);
		  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});

  // matrix for momentum and electric field source terms update
  MatrixXd matrix = MatrixXd::Constant(9,9,0.0);
  matrix(0,1) = r_i*B_z/l_r, matrix(0,2) = -r_i*B_y/l_r, matrix(0,6) = r_i*rho_i/l_r;
  matrix(1,0) = -r_i*B_z/l_r, matrix(1,2) = r_i*B_x/l_r, matrix(1,7) = r_i*rho_i/l_r;
  matrix(2,0) = r_i*B_y/l_r, matrix(2,1) = -r_i*B_x/l_r, matrix(2,8) = r_i*rho_i/l_r;
  matrix(3,4) = r_e*B_z/l_r, matrix(3,5) = -r_e*B_y/l_r, matrix(3,6) = r_e*rho_e/l_r;
  matrix(4,3) = -r_e*B_z/l_r, matrix(4,5) = r_e*B_x/l_r, matrix(4,7) = r_e*rho_e/l_r;
  matrix(5,3) = r_e*B_y/l_r, matrix(5,4) = -r_e*B_x/l_r, matrix(5,8) = r_e*rho_e/l_r;
  matrix(6,0) = -r_i/(lambda_d*lambda_d*l_r), matrix(6,3) = -r_e/(lambda_d*lambda_d*l_r);
  matrix(7,1) = -r_i/(lambda_d*lambda_d*l_r), matrix(7,4) = -r_e/(lambda_d*lambda_d*l_r);
  matrix(8,2) = -r_i/(lambda_d*lambda_d*l_r), matrix(8,5) = -r_e/(lambda_d*lambda_d*l_r);

  // calculate the matrix used in to invert
  MatrixXd finalMatrix = identity - dt*matrix;
  
  // LU decomposition, useful to calculate the inverse
  Eigen::FullPivLU<MatrixXd> lu_matrix(finalMatrix);
  // calculate inverse
  MatrixXd inverse_matrix = lu_matrix.inverse();
  VectorXd source_var(15);
  source_var << momX_i,momY_i,momZ_i,momX_e,momY_e,momZ_e,E_x,E_y,E_z,rho_i,rho_e,B_x,B_y,B_z,dt;
  VectorXd source_var_new = inverse_matrix*source_var;
  
  arr(i,j,k,0) = rho_i;
  arr(i,j,k,1) = source_var_new(0);
  arr(i,j,k,2) = source_var_new(1);
  arr(i,j,k,3) = source_var_new(2);

  arr(i,j,k,RHO_E) = rho_e;
  arr(i,j,k,MOMX_E) = source_var_new(3);
  arr(i,j,k,MOMY_E) = source_var_new(4);
  arr(i,j,k,MOMZ_E) = source_var_new(5);

  arr(i,j,k,BX) = B_x;
  arr(i,j,k,BY) = B_y;
  arr(i,j,k,BZ) = B_z;
  arr(i,j,k,EX) = source_var_new(6);
  arr(i,j,k,EY) = source_var_new(7);
  arr(i,j,k,EZ) = source_var_new(8);
  
  arr(i,j,k,ENER_I) = E_i + dt*r_i/l_r*(arr(i,j,k,EX)*arr(i,j,k,1) + arr(i,j,k,EY)*arr(i,j,k,2) + arr(i,j,k,EZ)*arr(i,j,k,3));
  arr(i,j,k,ENER_E) = E_e + dt*r_e/l_r*(arr(i,j,k,EX)*arr(i,j,k,MOMX_E) + arr(i,j,k,EY)*arr(i,j,k,MOMY_E) + arr(i,j,k,EZ)*arr(i,j,k,MOMZ_E));

}
void CAMReXmp::sourceUpdateANEX(Array4<Real>& arr, int i, int j, int k, Real dt)
{
  
  // define conserved variables                     
  Real rho_i = arr(i,j,k,0);
  Real momX_i = arr(i,j,k,1);
  Real momY_i = arr(i,j,k,2);
  Real momZ_i = arr(i,j,k,MOMZ_I);
  Real E_i = arr(i,j,k,ENER_I);
  Real rho_e = arr(i,j,k,RHO_E);
  Real momX_e = arr(i,j,k,MOMX_E);
  Real momY_e = arr(i,j,k,MOMY_E);
  Real momZ_e = arr(i,j,k,MOMZ_E);
  Real E_e = arr(i,j,k,ENER_E);
  Real B_x = arr(i,j,k,BX);
  Real B_y = arr(i,j,k,BY);
  Real B_z = arr(i,j,k,BZ);
  Real E_x = arr(i,j,k,EX);
  Real E_y = arr(i,j,k,EY);
  Real E_z = arr(i,j,k,EZ);
		  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});
  
  // useful constants
  Real B_squared = get_magnitude_squared(B_x, B_y, B_z);
  Real B = std::sqrt(B_squared);
  Real B_three_halves = B*B_squared;
  Real EdotB = E_x*B_x + E_y*B_y + E_z*B_z;
  Real v_idotB = v_x_i*B_x + v_y_i*B_y + v_z_i*B_z;
  Real v_edotB = v_x_e*B_x + v_y_e*B_y + v_z_e*B_z;

  Vector<Real> EcrossB = cross_product({E_x,E_y,E_z},{B_x,B_y,B_z});
  Vector<Real> v_icrossB = cross_product({v_x_i,v_y_i,v_z_i},{B_x,B_y,B_z});
  Vector<Real> v_ecrossB = cross_product({v_x_e,v_y_e,v_z_e},{B_x,B_y,B_z});
  Vector<Real> Bcross_EcrossB = cross_product({B_x,B_y,B_z},EcrossB);
  Vector<Real> Bcross_v_icrossB = cross_product({B_x,B_y,B_z},v_icrossB);
  Vector<Real> Bcross_v_ecrossB = cross_product({B_x,B_y,B_z},v_ecrossB);
  
  // exact solution of new velocities
  Vector<Real> v_i_new(3,0.0);
  Vector<Real> v_e_new(3,0.0);
  for (int dir=0; dir<3; dir++){
    v_i_new[dir] = v_exact_dir(B_squared, B, B_three_halves, EdotB, v_idotB,
			     EcrossB[dir], v_icrossB[dir], Bcross_EcrossB[dir], Bcross_v_icrossB[dir],
			     arr(i,j,k,BX+dir), dt, r_i);
    v_e_new[dir] = v_exact_dir(B_squared, B, B_three_halves, EdotB, v_edotB,
                             EcrossB[dir], v_ecrossB[dir], Bcross_EcrossB[dir], Bcross_v_ecrossB[dir],
                             arr(i,j,k,BX+dir), dt, r_e);
  }
  
  arr(i,j,k,0) = rho_i;
  arr(i,j,k,1) = rho_i*v_i_new[0];
  arr(i,j,k,2) = rho_i*v_i_new[1];
  arr(i,j,k,3) = rho_i*v_i_new[2];
  
  arr(i,j,k,RHO_E) = rho_e;
  arr(i,j,k,MOMX_E) = rho_e*v_e_new[0];
  arr(i,j,k,MOMY_E) = rho_e*v_e_new[1];
  arr(i,j,k,MOMZ_E) = rho_e*v_e_new[2];
  
  arr(i,j,k,ENER_I) = E_i + dt*r_i*rho_i/l_r*(E_x*v_i_new[0] + E_y*v_i_new[1] + E_z*v_i_new[2]);
  arr(i,j,k,ENER_E) = E_e + dt*r_e*rho_e/l_r*(E_x*v_e_new[0] + E_y*v_e_new[1] + E_z*v_e_new[2]);

  if (MaxwellMethod=="HYP")
    {
      arr(i,j,k,EX) = E_x - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_i_new[0] + r_e*rho_e*v_e_new[0]);
      arr(i,j,k,EY) = E_y - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_i_new[1] + r_e*rho_e*v_e_new[1]);
      arr(i,j,k,EZ) = E_z - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_i_new[2] + r_e*rho_e*v_e_new[2]);
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
void CAMReXmp::implicitMaxwellSolver(MultiFab& Sborder, const Real* dx, Real dt)
{
  
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
  Real b = c*c*dt*dt;

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
  computeRhs(Rhs, J, Sborder, dx, dt);

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
  implicitMaxwellSolverElecFieldUpdate(Sborder, J, dx, dt);
}
void CAMReXmp::implicitMaxwellSolverElecFieldUpdate(MultiFab& Sborder, const amrex::MultiFab& J, const Real* dx, Real dt)
{
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
	  const auto& J_nPlusHalf = J.array(mfi); 

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      Real dxBy = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,BX+(1+d)%3),arr(i+iOffset,j+jOffset,k+kOffset,BX+(1+d)%3),dx[d]);
		      Real dxBz = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,BX+(2+d)%3),arr(i+iOffset,j+jOffset,k+kOffset,BX+(2+d)%3),dx[d]);

		      // when using implicit source treatment do not include the current
		      if (sourceMethod=="IM")
			{
			  arr(i,j,k,EX+d) = arr(i,j,k,EX)/(Real)amrex::SpaceDim;
			  arr(i,j,k,EX+(1+d)%3) = arr(i,j,k,EX+(1+d)%3)/(Real)amrex::SpaceDim - dt*c*c*dxBz;
			  arr(i,j,k,EX+(2+d)%3) = arr(i,j,k,EX+(2+d)%3)/(Real)amrex::SpaceDim + dt*c*c*dxBy;
			}
		      else
			{
			  arr(i,j,k,EX+d) = arr(i,j,k,EX)/(Real)amrex::SpaceDim + dt*(-J_nPlusHalf(i,j,k,0)/((Real)amrex::SpaceDim*lambda_d*lambda_d*l_r));
			  arr(i,j,k,EX+(1+d)%3) = arr(i,j,k,EX+(1+d)%3)/(Real)amrex::SpaceDim + dt*(-J_nPlusHalf(i,j,k,1)/((Real)amrex::SpaceDim*lambda_d*lambda_d*l_r) - c*c*dxBz);
			  arr(i,j,k,EX+(2+d)%3) = arr(i,j,k,EX+(2+d)%3)/(Real)amrex::SpaceDim + dt*(-J_nPlusHalf(i,j,k,2)/((Real)amrex::SpaceDim*lambda_d*lambda_d*l_r) + c*c*dxBy);
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
		  //J_nPlusHalf(i,j,k,0) = computeJx_nPlusHalf(arr, i, j, k, dx[0], dx[1], dt);
		  //J_nPlusHalf(i,j,k,1) = computeJy_nPlusHalf(arr, i, j, k, dx[0], dx[1], dt);
		  //J_nPlusHalf(i,j,k,2) = computeJz_nPlusHalf(arr, i, j, k, dx[0], dx[1], dt);
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
  Rhs[0].define(grids, dmap, 1, 0);
  Rhs[1].define(grids, dmap, 1, 0);
  Rhs[2].define(grids, dmap, 1, 0);
  
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
		      // contribution for Rhs-Y (or corresponding cylic relations)
		      Real dxJz_nPlusHalf = computeDerivative(J_nPlusHalf(i-iOffset,j-jOffset,k-kOffset,(2+d)%3), J_nPlusHalf(i+iOffset,j+jOffset,k+kOffset,(2+d)%3), dx[d]);
		      Real dxEz = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,EX+(2+d)%3), arr(i+iOffset,j+jOffset,k+kOffset,EX+(2+d)%3), dx[d]);
		      // contribution for Rhs-Z (or corresponding cylic relations) 
		      Real dxJy_nPlusHalf = computeDerivative(J_nPlusHalf(i-iOffset,j-jOffset,k-kOffset,(1+d)%3), J_nPlusHalf(i+iOffset,j+jOffset,k+kOffset,(1+d)%3), dx[d]);
		      Real dxEy = computeDerivative(arr(i-iOffset,j-jOffset,k-kOffset,EX+(1+d)%3), arr(i+iOffset,j+jOffset,k+kOffset,EX+(1+d)%3), dx[d]);

		      // use (Real)amrex::SpaceDim to avoid overcounting in multiple dimensions
		      rhsArray[d] = arr(i,j,k,BX+d)/(Real)amrex::SpaceDim;
		      rhsArray[(1+d)%3] = arr(i,j,k,BX+(1+d)%3)/(Real)amrex::SpaceDim + dt*dxEz - dt*dt*dxJz_nPlusHalf/(lambda_d*lambda_d*l_r);
		      rhsArray[(2+d)%3] = arr(i,j,k,BX+(2+d)%3)/(Real)amrex::SpaceDim - dt*dxEy + dt*dt*dxJy_nPlusHalf/(lambda_d*lambda_d*l_r);
		      // update rhs variables
		      rhsX(i,j,k) = rhsArray[0];
		      rhsY(i,j,k) = rhsArray[1];
		      rhsZ(i,j,k) = rhsArray[2];
		      
		    }
		}
	    }
	}
    }
}

// For debugging purposes - this will print all the data from a multifab
void CAMReXmp::printComponents(MultiFab& mfIn)
{

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
		  //Real dyJz_nPlusHalf = computeDerivative(J_nPlusHalf(i,j-1,k,2), J_nPlusHalf(i,j+1,k,2), dy);
		  //Real dyEz = computeDerivative(arr(i,j-1,k,EZ), arr(i,j+1,k,EZ), dy);
		  //Real dxdxBx = computeSecondDerivative(arr(i-1,j,k,BX), arr(i,j,k,BX), arr(i+1,j,k,BX), dx);
		  //Real dydyBx = computeSecondDerivative(arr(i,j-1,k,BX), arr(i,j,k,BX), arr(i,j+1,k,BX), dx);		  
		  //rhs(i,j,k) = arr(i,j,k,BX) - dt*dyEz + dt*dt*dyJz_nPlusHalf/(lambda_d*lambda_d*l_r);
		  /////////////////////////////////////////////////////////////////////////
		  // 1D
		  rhs(i,j,k) = arr(i,j,k,BX);
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
		  //Real dyJx_nPlusHalf = computeDerivative(J_nPlusHalf(i,j-1,k,0), J_nPlusHalf(i,j+1,k,0), dy);
		  //Real dyEx = computeDerivative(arr(i,j-1,k,EX), arr(i,j+1,k,EX), dy);
		  ///////////////////////////////////////////////////////////////////////// 
		  Real dxJy_nPlusHalf = computeDerivative(J_nPlusHalf(i-1,j,k,1), J_nPlusHalf(i+1,j,k,1), dx);
		  Real dxEy = computeDerivative(arr(i-1,j,k,EY), arr(i+1,j,k,EY), dx);
		  //Real dxdxBz = computeSecondDerivative(arr(i-1,j,k,BZ), arr(i,j,k,BZ), arr(i+1,j,k,BZ), dx);
		  //Real dydyBz = computeSecondDerivative(arr(i,j-1,k,BZ), arr(i,j,k,BZ), arr(i,j+1,k,BZ), dx);
		  /////////////////////////////////////////////////////////////////////////
		  // 2D
		  //rhs(i,j,k) = arr(i,j,k,BZ) - dt*(dxEy-dyEx) + dt*dt*(dxJy_nPlusHalf-dyJx_nPlusHalf)/(lambda_d*lambda_d*l_r);
		  /////////////////////////////////////////////////////////////////////////
		  // 1D
		  rhs(i,j,k) = arr(i,j,k,BZ) - dt*dxEy + dt*dt*dxJy_nPlusHalf/(lambda_d*lambda_d*l_r); 
		}
	    }
	}
    }
}
