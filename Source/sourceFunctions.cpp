#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

//using namespace amrex;

// Eigen library
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

// 3x3 identity matrix
MatrixXd identity = MatrixXd::Identity(9, 9);

#include "sourceFunctionsStiff.H"

void CAMReXmp::sourceUpdate(Real dt, Real time)
{

  MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  MFIter::allowMultipleMFIters(true);
  
  for (MFIter mfi(S_input, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      
      Array4<Real> arr = S_input.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  (this->*sourceUpdateWithChosenMethod)(arr, i, j, k, dt);		  
		}
	    }
	}      
    }

  MultiFab& S_new = get_new_data(Phi_Type);
  MultiFab::Copy(S_new, S_input, 0, 0, NUM_STATE, 0);
  //FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

}
void CAMReXmp::sourceUpdateCyl(const Real* dx, Real dt, Real time)
{

  MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
  FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  if (MaxwellOrder!=0)
    {
      for (MFIter mfi(S_input, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);
      
	  Array4<Real> arr = S_input.array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      const Real x = geom.ProbLo()[0]+(double(i)+0.5)*dx[0];

		      //Real B_x = arr(i,j,k,BX);
		      //Real B_y = arr(i,j,k,BY);
		      Real B_z = arr(i,j,k,BZ);
		      //Real E_x = arr(i,j,k,EX);
		      //Real E_y = arr(i,j,k,EY);
		      Real E_z = arr(i,j,k,EZ);
		      
		      // z-component (y-) needs a theta-component (z-) update
		      // use negative sign because theta-component is reversed??
		      arr(i,j,k,BY) += -dt*(-E_z/x);
		      arr(i,j,k,EY) += -dt*(c*c*B_z/x);
		    }
		}
	    }      
	}

      MultiFab& S_new = get_new_data(Phi_Type);
      MultiFab::Copy(S_new, S_input, 0, 0, NUM_STATE, 0); 
      //MultiFab::Copy(S_new, S_input, BY, BY, 1, 0);
      //MultiFab::Copy(S_new, S_input, EY, EY, 1, 0);
      FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

#if (AMREX_SPACEDIM >= 2)       
      Array<MultiFab,AMREX_SPACEDIM> S_EM_input;
      S_EM_input[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
      FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);
  
      for (MFIter mfi(S_EM_input[1], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);
      
	  // Indexable arrays for the data, and the directional flux
	  // Based on the corner-centred definition of the flux array, the
	  // data array runs from e.g. [0,N+1] and the flux array from [-1,N+1]
	  //const auto& arr = S_dest.array(mfi);
	  const auto& arrEM = S_EM_input[1].array(mfi);
	  const auto& arr = S_input.array(mfi);

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      arrEM(i,j,k,BY_LOCAL) = 0.5*(arr(i,j-1,k,BY)+arr(i,j,k,BY));
		      // electric field is updated using elecFieldCellAve() in the next step
		      arrEM(i,j,k,EY_LOCAL) = 0.5*(arr(i,j-1,k,EY)+arr(i,j,k,EY));
		    }
		}
	    }      
	}

      MultiFab& S_EM_Y_new = get_new_data(EM_Y_Type);
      MultiFab::Copy(S_EM_Y_new, S_EM_input[1], BY_LOCAL, BY_LOCAL, 1, 0);
      MultiFab::Copy(S_EM_Y_new, S_EM_input[1], EY_LOCAL, EY_LOCAL, 1, 0);
#endif
    }

  if (fluidOrder!=0)
    {

      MFIter::allowMultipleMFIters(true);
  
      for (MFIter mfi(S_input, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
      
	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);
      
	  Array4<Real> arr = S_input.array(mfi);
      
	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      const Real x = geom.ProbLo()[0]+(double(i)+0.5)*dx[0];

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
		  
		      // define primitive variables
		      Real v_x_i = momX_i/rho_i;
		      Real v_y_i = momY_i/rho_i;
		      //Real v_z_i = momZ_i/rho_i;
		      Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
		      Real v_x_e = momX_e/rho_e;
		      Real v_y_e = momY_e/rho_e;
		      //Real v_z_e = momZ_e/rho_e;
		      Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});
		      
		      arr(i,j,k,0) -= dt*(rho_i*v_x_i/x);
		      arr(i,j,k,1) -= dt*(rho_i*v_x_i*v_x_i/x);
		      arr(i,j,k,2) -= dt*(rho_i*v_x_i*v_y_i/x);
		      arr(i,j,k,ENER_I) -= dt*((E_i+p_i)*v_x_i/x);
    
		      arr(i,j,k,RHO_E) -= dt*(rho_e*v_x_e/x);
		      arr(i,j,k,MOMX_E) -= dt*(rho_e*v_x_e*v_x_e/x);
		      arr(i,j,k,MOMY_E) -= dt*(rho_e*v_x_e*v_y_e/x);
		      arr(i,j,k,ENER_E) -= dt*((E_e+p_e)*v_x_e/x);
		    }
		}
	    }      
	}

      MultiFab& S_new = get_new_data(Phi_Type);
      MultiFab::Copy(S_new, S_input, 0, 0, NUM_STATE_FLUID, 0);
    }  
}
void CAMReXmp::sourceUpdateVoid(Array4<Real>& arr, int i, int j, int k, Real dt)
{
  amrex::Abort("No source terms!!");
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
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;

  //arr(i,j,k,EZ) = E_z - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_z_i + r_e*rho_e*v_z_e);
  //E_z = arr(i,j,k,EZ);
		  
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
  
  //if (MaxwellMethod=="HYP")
    {
      //arr(i,j,k,EX) = E_x - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_x_i + r_e*rho_e*v_x_e);
      //arr(i,j,k,EY) = E_y - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_y_i + r_e*rho_e*v_y_e);
      //arr(i,j,k,EZ) = E_z - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_z_i + r_e*rho_e*v_z_e);
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
  VectorXd source_var(9);
  source_var << momX_i,momY_i,momZ_i,momX_e,momY_e,momZ_e,E_x,E_y,E_z;//,rho_i,rho_e,B_x,B_y,B_z,dt;
  VectorXd source_var_new = inverse_matrix*source_var;  
  /*
  // matrix for momentum source terms update
  MatrixXd matrix_i = MatrixXd::Constant(3,3,0.0);
  MatrixXd matrix_e = MatrixXd::Constant(3,3,0.0);
  matrix_i(0,1) = B_z/l_r, matrix_i(0,2) = -B_y/l_r; //, matrix(0,6) = rho_i/l_r;
  matrix_i(1,0) = -B_z/l_r, matrix_i(1,2) = B_x/l_r; //, matrix(1,7) = rho_i/l_r;
  matrix_i(2,0) = B_y/l_r, matrix_i(2,1) = -B_x/l_r; //, matrix(2,8) = rho_i/l_r;
  matrix_e(0,1) = -m*B_z/l_r, matrix_e(0,2) = m*B_y/l_r; //, matrix(3,6) = -m*rho_e/l_r;
  matrix_e(1,0) = m*B_z/l_r, matrix_e(1,2) = -m*B_x/l_r; //, matrix(4,7) = -m*rho_e/l_r;
  matrix_e(2,0) = -m*B_y/l_r, matrix_e(2,1) = m*B_x/l_r; //, matrix(5,8) = -m*rho_e/l_r;

  // calculate the matrix used in to invert
  MatrixXd finalMatrix_i = identity - dt*matrix_i;
  MatrixXd finalMatrix_e = identity - dt*matrix_e;
  
  // LU decomposition, useful to calculate the inverse
  //Eigen::FullPivLU<MatrixXd> lu_matrix_i(finalMatrix_i);
  //Eigen::FullPivLU<MatrixXd> lu_matrix_e(finalMatrix_e);
  // calculate inverse
  //MatrixXd inverse_matrix_i = lu_matrix_i.inverse();
  //MatrixXd inverse_matrix_e = lu_matrix_e.inverse();
  MatrixXd inverse_matrix_i = finalMatrix_i.inverse();
  MatrixXd inverse_matrix_e = finalMatrix_e.inverse();
  
  VectorXd source_var_i(3), source_var_e(3);
  source_var_i << momX_i+dt*E_x*rho_i/l_r,momY_i+dt*E_y*rho_i/l_r,momZ_i+dt*E_z*rho_i/l_r;
  source_var_e << momX_e-m*dt*E_x*rho_e/l_r,momY_e-m*dt*E_y*rho_e/l_r,momZ_e-m*dt*E_z*rho_e/l_r;
  VectorXd source_var_new_i = inverse_matrix_i*source_var_i;
  VectorXd source_var_new_e = inverse_matrix_e*source_var_e;
  */
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
void CAMReXmp::sourceUpdateIMMidpoint(Array4<Real>& arr, int i, int j, int k, Real dt)
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
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;

  // get kinetic energies
  Real v_i = get_magnitude(v_x_i, v_y_i, v_z_i);
  Real v_e = get_magnitude(v_x_e, v_y_e, v_z_e);
  Real kin_i = 0.5*rho_i*v_i*v_i;
  Real kin_e = 0.5*rho_e*v_e*v_e;
  
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
  MatrixXd finalMatrix = identity - 0.5*dt*matrix;

  // LU decomposition, useful to calculate the inverse
  Eigen::FullPivLU<MatrixXd> lu_matrix(finalMatrix);
  // calculate inverse
  MatrixXd inverse_matrix = lu_matrix.inverse();
  VectorXd source_var(9);
  source_var << momX_i,momY_i,momZ_i,momX_e,momY_e,momZ_e,E_x,E_y,E_z;//,rho_i,rho_e,B_x,B_y,B_z,dt;
  VectorXd source_var_new = inverse_matrix*source_var;  

  // flux function
  Vector<Real> functionInt(NUM_STATE,0.0);
  functionInt[0] = rho_i;
  functionInt[1] = source_var_new(0);
  functionInt[2] = source_var_new(1);
  functionInt[3] = source_var_new(2);

  functionInt[RHO_E] = rho_e;
  functionInt[MOMX_E] = source_var_new(3);
  functionInt[MOMY_E] = source_var_new(4);
  functionInt[MOMZ_E] = source_var_new(5);

  functionInt[BX] = B_x;
  functionInt[BY] = B_y;
  functionInt[BZ] = B_z;
  functionInt[EX] = source_var_new(6);
  functionInt[EY] = source_var_new(7);
  functionInt[EZ] = source_var_new(8);
  
  functionInt[ENER_I] = E_i;
  functionInt[ENER_E] = E_e;
  
  for (int n=0; n<NUM_STATE; n++)
    arr(i,j,k,n) = 2.0*functionInt[n]-arr(i,j,k,n);

  // get new kinetic energies
  Real mom_i_new = get_magnitude(arr(i,j,k,MOMX_I), arr(i,j,k,MOMY_I), arr(i,j,k,MOMZ_I));
  Real mom_e_new = get_magnitude(arr(i,j,k,MOMX_E), arr(i,j,k,MOMY_E), arr(i,j,k,MOMZ_E));
  Real kin_i_new = 0.5*mom_i_new*mom_i_new/rho_i;
  Real kin_e_new = 0.5*mom_e_new*mom_e_new/rho_e;

  arr(i,j,k,ENER_I) += kin_i_new-kin_i;
  arr(i,j,k,ENER_E) += kin_e_new-kin_e;

#if (AMREX_SPACEDIM >= 2)
  if (MaxwellDivMethod=="HDC")
    {
      Real psi_e = arr(i,j,k,DIVE);
      arr(i,j,k,DIVE) = psi_e + ce*dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i + r_e*rho_e);
    }
#endif  
  
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
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  
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

  // Resistivity
  /*Real currentX = r_i*rho_i*v_i_new[0] + r_e*rho_e*v_e_new[0];
  Real currentY = r_i*rho_i*v_i_new[1] + r_e*rho_e*v_e_new[1];
  Real currentZ = r_i*rho_i*v_i_new[2] + r_e*rho_e*v_e_new[2];
  
  arr(i,j,k,1) -= dt/l_r*eta*r_i*rho_i*currentX;
  arr(i,j,k,2) -= dt/l_r*eta*r_i*rho_i*currentY;
  arr(i,j,k,3) -= dt/l_r*eta*r_i*rho_i*currentZ;
  
  arr(i,j,k,MOMX_E) -= dt/l_r*eta*r_e*rho_e*currentX;
  arr(i,j,k,MOMY_E) -= dt/l_r*eta*r_e*rho_e*currentY;
  arr(i,j,k,MOMZ_E) -= dt/l_r*eta*r_e*rho_e*currentZ;
  */
  //if (MaxwellMethod=="HYP")
    {
      //arr(i,j,k,EX) = E_x - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_i_new[0] + r_e*rho_e*v_e_new[0]);
      //arr(i,j,k,EY) = E_y - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_i_new[1] + r_e*rho_e*v_e_new[1]);
      //arr(i,j,k,EZ) = E_z - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_i_new[2] + r_e*rho_e*v_e_new[2]);
      //arr(i,j,k,EX) = E_x - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_x_i + r_e*rho_e*v_x_e);
      //arr(i,j,k,EY) = E_y - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_y_i + r_e*rho_e*v_y_e);
      //arr(i,j,k,EZ) = E_z - dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_z_i + r_e*rho_e*v_z_e);


      /*
      arr(i,j,k,EX) -= dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_i_new[0] + r_e*rho_e*v_e_new[0]);
      arr(i,j,k,EY) -= dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_i_new[1] + r_e*rho_e*v_e_new[1]);
      arr(i,j,k,EZ) -= dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i*v_i_new[2] + r_e*rho_e*v_e_new[2]);

      arr(i,j,k,EX) = 0.5*E_x + 0.5*arr(i,j,k,EX);
      arr(i,j,k,EY) = 0.5*E_y + 0.5*arr(i,j,k,EY);
      arr(i,j,k,EZ) = 0.5*E_z + 0.5*arr(i,j,k,EZ);
      */
      
    }

}
void CAMReXmp::sourceUpdateStiff(Array4<Real>& arr, int i, int j, int k, Real dt)
{

  // define conserved variables                     
  Real rho_i = arr(i,j,k,0);
  Real momX_i = arr(i,j,k,1);
  Real momY_i = arr(i,j,k,2);
  Real momZ_i = arr(i,j,k,MOMZ_I);
  //Real E_i = arr(i,j,k,ENER_I);
  Real rho_e = arr(i,j,k,RHO_E);
  Real momX_e = arr(i,j,k,MOMX_E);
  Real momY_e = arr(i,j,k,MOMY_E);
  Real momZ_e = arr(i,j,k,MOMZ_E);
  //Real E_e = arr(i,j,k,ENER_E);
  /*Real B_x = arr(i,j,k,BX);
  Real B_y = arr(i,j,k,BY);
  Real B_z = arr(i,j,k,BZ);
  Real E_x = arr(i,j,k,EX);
  Real E_y = arr(i,j,k,EY);
  Real E_z = arr(i,j,k,EZ);
  */
  // define velocities
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;

  // define input
  /*std::vector<double> u_i = {rho_i,v_x_i,v_y_i,v_z_i,E_i,
			     B_x,B_y,B_z,E_x,E_y,E_z};

  size_t num_of_steps = stiffSolver(u_i, r_i, l_r, dt);

  // define input
  std::vector<double> u_e = {rho_e,v_x_e,v_y_e,v_z_e,E_e,
			     B_x,B_y,B_z,E_x,E_y,E_z};
  size_t num_of_steps_e = stiffSolver(u_e, r_e, l_r, dt);

  arr(i,j,k,1) = rho_i*u_i[1];
  arr(i,j,k,2) = rho_i*u_i[2];
  arr(i,j,k,3) = rho_i*u_i[3];
  arr(i,j,k,ENER_I) = u_i[4];
  
  arr(i,j,k,MOMX_E) = rho_e*u_e[1];
  arr(i,j,k,MOMY_E) = rho_e*u_e[2];
  arr(i,j,k,MOMZ_E) = rho_e*u_e[3];
  arr(i,j,k,ENER_E) = u_e[4];
  */

  // get kinetic energies
  Real v_i = get_magnitude(v_x_i, v_y_i, v_z_i);
  Real v_e = get_magnitude(v_x_e, v_y_e, v_z_e);
  Real kin_i = 0.5*rho_i*v_i*v_i;
  Real kin_e = 0.5*rho_e*v_e*v_e;

  std::vector<double> u = get_data_zone(arr,i,j,k,0,NUM_STATE);

  size_t num_of_steps = stiffSolverFull(u, r_i, r_e, l_r, lambda_d, dt);
  for (int n=0; n<NUM_STATE; n++)
    arr(i,j,k,n) = u[n];

  // get new kinetic energies
  Real mom_i_new = get_magnitude(arr(i,j,k,MOMX_I), arr(i,j,k,MOMY_I), arr(i,j,k,MOMZ_I));
  Real mom_e_new = get_magnitude(arr(i,j,k,MOMX_E), arr(i,j,k,MOMY_E), arr(i,j,k,MOMZ_E));
  Real kin_i_new = 0.5*mom_i_new*mom_i_new/rho_i;
  Real kin_e_new = 0.5*mom_e_new*mom_e_new/rho_e;

  arr(i,j,k,ENER_I) += kin_i_new-kin_i;
  arr(i,j,k,ENER_E) += kin_e_new-kin_e;

#if (AMREX_SPACEDIM >= 2)
  if (MaxwellDivMethod=="HDC")
    {
      Real psi_e = arr(i,j,k,DIVE);
      arr(i,j,k,DIVE) = psi_e + ce*dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*rho_i + r_e*rho_e);
    }
#endif    
}
