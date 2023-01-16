#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

// Eigen library
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

// 9x9 identity matrix
MatrixXd identity = MatrixXd::Identity(9, 9);

void CAMReXmp::sourceUpdate(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  MFIter::allowMultipleMFIters(true);
  
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
		  if (geom.Coord()==1)		  
		    {
		      //Real x = geom.ProbLo()[0]+(double(i)+0.5)*dx[0];
		      //const Real x = (double(i)+0.5)*dx[0];
		      // y <- x in cylindrical
		      const Real y = geom.ProbLo()[1]+(double(j)+0.5)*dx[1];
		      //amrex::Print() << x << " ";
		      cylSourceUpdate(arr,i,j,k,dt,y);
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

void CAMReXmp::cylSourceUpdateImplicitRK2(MultiFab& Sborder, MultiFab& Sold, const Real* dx, Real dt)
{
MFIter::allowMultipleMFIters(true);
  
  for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      
      Array4<Real> arr = Sborder.array(mfi);
      Array4<Real> arr_old = Sold.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  const Real y = geom.ProbLo()[1]+(double(j)+0.5)*dx[1];
		  //arr(i,j,k,BX) -= (1.0-tau)*dt*(arr_old(i,j,k,EZ)/y);      
		  arr(i,j,k,EX) += (1.0-tau)*dt*(c*c*arr_old(i,j,k,BZ)/y); 
				  	    
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

void CAMReXmp::cylSourceUpdate(Array4<Real>& arr, int i, int j, int k, Real dt, Real y)
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
  /*
  arr(i,j,k,0) -= dt*(rho_i*v_x_i/x);
  arr(i,j,k,1) -= dt*(rho_i*v_x_i*v_x_i/x);
  arr(i,j,k,2) -= dt*(rho_i*v_x_i*v_y_i/x);
  arr(i,j,k,3) -= dt*(rho_i*v_x_i*v_z_i/x);
  arr(i,j,k,ENER_I) -= dt*((E_i+p_i)*v_x_i/x);
    
  arr(i,j,k,RHO_E) -= dt*(rho_e*v_x_e/x);
  arr(i,j,k,MOMX_E) -= dt*(rho_e*v_x_e*v_x_e/x);
  arr(i,j,k,MOMY_E) -= dt*(rho_e*v_x_e*v_y_e/x);
  arr(i,j,k,MOMZ_E) -= dt*(rho_e*v_x_e*v_z_e/x);
  arr(i,j,k,ENER_E) -= dt*((E_e+p_e)*v_x_e/x);
  
  if (MaxwellMethod=="HYP")
    {
      arr(i,j,k,BZ) -= dt*(E_y/x);      
      arr(i,j,k,EZ) += dt*(c*c*B_y/x);
    }
  */
  // y <- x in cylindrical
  // v_y <- v_x in cylindrical  
  arr(i,j,k,0) -= dt*(rho_i*v_y_i/y);
  arr(i,j,k,1) -= dt*(rho_i*v_y_i*v_x_i/y);
  arr(i,j,k,2) -= dt*(rho_i*v_y_i*v_y_i/y);
  arr(i,j,k,3) -= dt*(rho_i*v_y_i*v_z_i/y);
  arr(i,j,k,ENER_I) -= dt*((E_i+p_i)*v_y_i/y);
  
  arr(i,j,k,RHO_E) -= dt*(rho_e*v_y_e/y);
  arr(i,j,k,MOMX_E) -= dt*(rho_e*v_y_e*v_x_e/y);
  arr(i,j,k,MOMY_E) -= dt*(rho_e*v_y_e*v_y_e/y);
  arr(i,j,k,MOMZ_E) -= dt*(rho_e*v_y_e*v_z_e/y);
  arr(i,j,k,ENER_E) -= dt*((E_e+p_e)*v_y_e/y);
  
  if (MaxwellMethod=="IM")
    {
      arr(i,j,k,BX) -= tau*dt*(E_z/y);      
      arr(i,j,k,EX) += tau*dt*(c*c*B_z/y);
    }
  else
    {
      arr(i,j,k,BX) -= dt*(E_z/y);      
      arr(i,j,k,EX) += dt*(c*c*B_z/y);    
    }
}
