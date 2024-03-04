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
MatrixXd identityMat = MatrixXd::Identity(9, 9);

void CAMReXmp::sourceUpdateIMMidpointGOL(Array4<Real>& arr, int i, int j, int k, Real dt)
{

  // define conserved variables                     
  Real rho = arr(i,j,k,0);
  Real momX = arr(i,j,k,1);
  Real momY = arr(i,j,k,2);
  Real momZ = arr(i,j,k,MOMZ_I);
  Real E = arr(i,j,k,ENER_I);
  Real J_x = arr(i,j,k,MOMX_E);
  Real J_y = arr(i,j,k,MOMY_E);
  Real J_z = arr(i,j,k,MOMZ_E);
  Real E_e = arr(i,j,k,ENER_E);
  Real B_x = arr(i,j,k,BX);
  Real B_y = arr(i,j,k,BY);
  Real B_z = arr(i,j,k,BZ);
  Real E_x = arr(i,j,k,EX);
  Real E_y = arr(i,j,k,EY);
  Real E_z = arr(i,j,k,EZ);
  
  // electron variables
  Vector<Real> u_e = get_electron_var(get_data_zone(arr,i,j,k,0,10));
  Real rho_e = u_e[0];
  Real momX_e = u_e[1];
  Real momY_e = u_e[2];
  Real momZ_e = u_e[3];

  // quasineutrality n_i=n_e, rho_e=rho_i*(m_e/m_i) and m_i=1
  Real n_e = rho;
  Real m_e = 1.0/m;
  Real q_i = r_i;
  
  // define primitive variables
  Real v_x = momX/rho;
  Real v_y = momY/rho;
  Real v_z = momZ/rho;
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;

  // get kinetic energies
  Real v = get_magnitude(v_x, v_y, v_z);
  Real v_e = get_magnitude(v_x_e, v_y_e, v_z_e);
  Real kin = 0.5*rho*v*v;
  Real kin_e = 0.5*rho_e*v_e*v_e;
  Real J = get_magnitude(J_x, J_y, J_z);
  Real kin_J = 0.5*m_e/(q_i*q_i*n_e)*J*J;

  // matrix for momentum, current and electric field source terms update
  MatrixXd matrix = MatrixXd::Constant(9,9,0.0);
  matrix(0,1+3) = B_z/l_r, matrix(0,2+3) = -B_y/l_r;
  matrix(1,0+3) = -B_z/l_r, matrix(1,2+3) = B_x/l_r;
  matrix(2,0+3) = B_y/l_r, matrix(2,1+3) = -B_x/l_r;
  matrix(3,4) = r_e*B_z/l_r, matrix(3,5) = -r_e*B_y/l_r, matrix(3,6) = r_e*r_e*rho_e/l_r;
  matrix(4,3) = -r_e*B_z/l_r, matrix(4,5) = r_e*B_x/l_r, matrix(4,7) = r_e*r_e*rho_e/l_r;
  matrix(5,3) = r_e*B_y/l_r, matrix(5,4) = -r_e*B_x/l_r, matrix(5,8) = r_e*r_e*rho_e/l_r;
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // additional entries for the current
  matrix(3,4-3) = r_e*r_e*B_z/(m*l_r), matrix(3,5-3) = -r_e*r_e*B_y/(m*l_r);
  matrix(4,3-3) = -r_e*r_e*B_z/(m*l_r), matrix(4,5-3) = r_e*r_e*B_x/(m*l_r);
  matrix(5,3-3) = r_e*r_e*B_y/(m*l_r), matrix(5,4-3) = -r_e*r_e*B_x/(m*l_r);
  ////////////////////////////////////////////////////////////////////////////////////////////////
  matrix(6,3) = -1.0/(lambda_d*lambda_d*l_r);
  matrix(7,4) = -1.0/(lambda_d*lambda_d*l_r);
  matrix(8,5) = -1.0/(lambda_d*lambda_d*l_r);

  // calculate the matrix used in to invert
  MatrixXd finalMatrix = identityMat - 0.5*dt*matrix;

  // LU decomposition, useful to calculate the inverse
  Eigen::FullPivLU<MatrixXd> lu_matrix(finalMatrix);
  // calculate inverse
  MatrixXd inverse_matrix = lu_matrix.inverse();
  VectorXd source_var(9);
  source_var << momX,momY,momZ,J_x,J_y,J_z,E_x,E_y,E_z;
  VectorXd source_var_new = inverse_matrix*source_var;  

  // flux function
  Vector<Real> functionInt(NUM_STATE,0.0);
  functionInt[0] = rho;
  functionInt[1] = source_var_new(0);
  functionInt[2] = source_var_new(1);
  functionInt[3] = source_var_new(2);

  functionInt[RHO_E] = 0.0;
  functionInt[MOMX_E] = source_var_new(3);
  functionInt[MOMY_E] = source_var_new(4);
  functionInt[MOMZ_E] = source_var_new(5);

  functionInt[BX] = B_x;
  functionInt[BY] = B_y;
  functionInt[BZ] = B_z;
  functionInt[EX] = source_var_new(6);
  functionInt[EY] = source_var_new(7);
  functionInt[EZ] = source_var_new(8);
  
  functionInt[ENER_I] = E;
  functionInt[ENER_E] = E_e;
  //functionInt[ENER_I] = E + 0.5*dt/l_r*(functionInt[EX]*functionInt[MOMX_E] + functionInt[EY]*functionInt[MOMY_E] + functionInt[EZ]*functionInt[MOMZ_E]);
  //functionInt[ENER_E] = E_e + 0.5*dt/l_r*(functionInt[EX]*(functionInt[MOMX_E]-q_i*functionInt[1]) + functionInt[EY]*(functionInt[MOMY_E]-q_i*functionInt[2]) + functionInt[EZ]*(functionInt[MOMZ_E]-q_i*functionInt[3]));

  for (int n=0; n<NUM_STATE; n++)
    arr(i,j,k,n) = 2.0*functionInt[n]-arr(i,j,k,n);

  // get new electron velocities
  v_x = arr(i,j,k,MOMX_I)/rho;
  v_y = arr(i,j,k,MOMY_I)/rho;
  v_z = arr(i,j,k,MOMZ_I)/rho;
  v_x_e = v_x-arr(i,j,k,MOMX_E)/(n_e*q_i);
  v_y_e = v_y-arr(i,j,k,MOMY_E)/(n_e*q_i);
  v_z_e = v_z-arr(i,j,k,MOMZ_E)/(n_e*q_i);

  // get kinetic energies
  Real v_new = get_magnitude(v_x, v_y, v_z);
  Real v_e_new = get_magnitude(v_x_e, v_y_e, v_z_e);
  Real kin_new = 0.5*rho*v_new*v_new;
  Real kin_e_new = 0.5*rho_e*v_e_new*v_e_new;  
  Real J_new = get_magnitude(arr(i,j,k,MOMX_E), arr(i,j,k,MOMY_E), arr(i,j,k,MOMZ_E));
  Real kin_J_mew = 0.5*m_e/(q_i*q_i*n_e)*J_new*J_new;

  // see eq. 5.48 from the Sandia report of Richard M.J.Kramer
  // entitled "A Plasma Modeling Hierarchy and Verification Approach"
  arr(i,j,k,ENER_I) += kin_new-kin + kin_J_mew-kin_J;
  arr(i,j,k,ENER_E) += kin_e_new-kin_e;
 
}
