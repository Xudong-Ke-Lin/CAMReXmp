#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#include "utils.H"

Vector<Real> fluidGOLFlux(const Vector<Real>& u_i, int d){
  // define conserved variables
  Real rho_i = u_i[0];
  Real momX_i = u_i[1+d];
  Real momY_i = u_i[1+(1+d)%3];
  Real momZ_i = u_i[1+(2+d)%3];
  Real E_i = u_i[ENER_I];
  Real J_x = u_i[MOMX_E+d];
  Real J_y = u_i[MOMX_E+(1+d)%3];
  Real J_z = u_i[MOMX_E+(2+d)%3];
  Real E_e = u_i[ENER_E];

  // quasineutrality n_i=n_e_rho_i/m_i, rho_e=rho_i*(m_e/m_i) and m_i=1
  Real n_e = rho_i;
  Real rho_e = rho_i/m;
  Real m_e = 1.0/m;
  Real q_i = r_i;
  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});  
  Real v_x_e = v_x_i-J_x/(n_e*q_i);
  Real v_y_e = v_y_i-J_y/(n_e*q_i);
  Real v_z_e = v_z_i-J_z/(n_e*q_i);
  Real p_e = get_pressure({rho_e,rho_e*v_x_e,rho_e*v_y_e,rho_e*v_z_e,E_e});

  // flux function
  Vector<Real> function(u_i.size(),0.0);
  function[0] = rho_i*v_x_i;
  function[1+d] = rho_i*v_x_i*v_x_i + p_i+p_e;
  function[1+(1+d)%3] = rho_i*v_x_i*v_y_i;
  function[1+(2+d)%3] = rho_i*v_x_i*v_z_i;
  function[ENER_I] = (E_i+p_i)*v_x_i;
  function[MOMX_E+d] = 2.0*J_x*v_x_i - J_x*J_x/(n_e*q_i) -p_e*q_i/m_e;
  function[MOMX_E+(1+d)%3] = v_x_i*J_y + v_y_i*J_x - J_x*J_y/(n_e*q_i);
  function[MOMX_E+(2+d)%3] = v_x_i*J_z + v_z_i*J_x - J_x*J_z/(n_e*q_i);  
  function[ENER_E] = (E_e+p_e)*v_x_e;
  return function;
}
Vector<Real> fluidTWOFlux(const Vector<Real>& u_i, int d){
  // define conserved variables
  Real rho_i = u_i[0];
  Real momX_i = u_i[1+d];
  Real momY_i = u_i[1+(1+d)%3];
  Real momZ_i = u_i[1+(2+d)%3];
  Real E_i = u_i[ENER_I];
  Real rho_e = u_i[RHO_E];
  Real momX_e = u_i[MOMX_E+d];
  Real momY_e = u_i[MOMX_E+(1+d)%3];
  Real momZ_e = u_i[MOMX_E+(2+d)%3];  
  Real E_e = u_i[ENER_E];

  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});

  // flux functiono
  Vector<Real> function(u_i.size(),0.0);
  function[0] = rho_i*v_x_i;
  function[1+d] = rho_i*v_x_i*v_x_i + p_i;
  function[1+(1+d)%3] = rho_i*v_x_i*v_y_i;
  function[1+(2+d)%3] = rho_i*v_x_i*v_z_i;
  function[ENER_I] = (E_i+p_i)*v_x_i;
  function[RHO_E] = rho_e*v_x_e;
  function[MOMX_E+d] = rho_e*v_x_e*v_x_e + p_e;
  function[MOMX_E+(1+d)%3] = rho_e*v_x_e*v_y_e;
  function[MOMX_E+(2+d)%3] = rho_e*v_x_e*v_z_e;  
  function[ENER_E] = (E_e+p_e)*v_x_e;
  return function;
}
Real get_max_seed(const Vector<Real>& u_i, int d)
{
  // define conserved variables
  Real rho_i = u_i[0];
  Real momX_i = u_i[1+d];
  Real momY_i = u_i[1+(1+d)%3];
  Real momZ_i = u_i[1+(2+d)%3];
  Real E_i = u_i[ENER_I];
  Real J_x = u_i[MOMX_E+d];
  Real J_y = u_i[MOMX_E+(1+d)%3];
  Real J_z = u_i[MOMX_E+(2+d)%3];
  Real E_e = u_i[ENER_E];

  // quasineutrality n_i=n_e_rho_i/m_i, rho_e=rho_i*(m_e/m_i) and m_i=1
  Real n_e = rho_i;
  Real rho_e = rho_i/m;
  Real m_e = 1.0/m;
  Real q_i = r_i;
  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real v_x_e = v_x_i-J_x/(n_e*q_i);
  Real v_y_e = v_y_i-J_y/(n_e*q_i);
  Real v_z_e = v_z_i-J_z/(n_e*q_i);
  Vector<Real> u_e = {rho_e,rho_e*v_x_e,rho_e*v_y_e,rho_e*v_z_e,E_e};
  Real a = get_speed(u_e);
  return std::abs(v_x_e)+a;
  
  //return std::abs(v_x_e)+get_speed({rho_e,rho_e*v_x_e,rho_e*v_y_e,rho_e*v_z_e,E_e});
  //return get_speed({rho_e,rho_e*v_x_e,rho_e*v_y_e,rho_e*v_z_e,E_e});
  
}
Vector<Real> flux_Rusanov(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, Real dx, Real dt, int d,
			  std::function<Vector<Real> (const Vector<Real>&, int)> flux_function){
  Vector<Real> flux, func_i, func_iPlus1;
  func_i = flux_function(u_i,d);
  func_iPlus1 = flux_function(u_iPlus1,d);
  Real speed_max = std::max(get_max_seed(u_i,d),get_max_seed(u_iPlus1,d));
  for (int i = 0; i<u_i.size(); i++)
    {
      flux.push_back(0.5*(func_iPlus1[i]+func_i[i]
			  - speed_max*(u_iPlus1[i]-u_i[i])));
    }
  
  return flux;

}
Vector<Real> TVD_fluxGOL(const Array4<Real>& arr, const Array4<Real>& slopes,
			 int i, int j, int k, int iOffset, int jOffset, int kOffset,
			 int start, int len, int startSlope, Real dx, Real dt, int d){

  Vector<Real> u_iMinus1, u_i;

  Vector<Real> slopes_iMinus1, slopes_i;

  for (int n = start; n<start+len; n++)
    {
      u_iMinus1.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_i.push_back(arr(i,j,k,n));
    }
  for (int n = startSlope; n<startSlope+len; n++)
    {

      slopes_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n));
      slopes_i.push_back(slopes(i,j,k,n)); 
    }

  Vector<Real> u_iMinus1R,u_iL;
  for (int n = 0; n<len; n++)
    {
      u_iMinus1R.push_back(u_iMinus1[n] + 0.5*slopes_iMinus1[n]);
      u_iL.push_back(u_i[n] - 0.5*slopes_i[n]);
    }

  // flux depending on the speeds, defined in slides or Toro's book
  Vector<Real> flux = flux_FORCE(u_iMinus1R,u_iL,dx,dt,d,fluidGOLFlux);
  //Vector<Real> flux = flux_Rusanov(u_iMinus1R,u_iL,dx,dt,d,fluidGOLFlux);
  //Vector<Real> flux = flux_MUSTA1(u_iMinus1R,u_iL,dx,dt,d,0.8,fluidGOLFlux);
  //Vector<Real> flux = flux_MUSTA1(u_iMinus1R,u_iL,dx,dt,d,0.8);
  return flux;
}

Vector<Real> flux_LF(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, double dt, int d,
		     std::function<Vector<Real> (const Vector<Real>&, int)> flux_function){
    
    Vector<Real> flux, func_i, func_iPlus1;
    func_i = flux_function(u_i,d);
    func_iPlus1 = flux_function(u_iPlus1,d);
    for (int i = 0; i<u_i.size(); i++)
      {
        flux.push_back(0.5*(dx/dt)*(u_i[i]-u_iPlus1[i])
            + 0.5*(func_iPlus1[i]+func_i[i]));       
      }
    return flux;
}
Vector<Real> flux_RI(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, double dt, int d,
		     std::function<Vector<Real> (const Vector<Real>&, int)> flux_function){

    Vector<Real> u_half, func_i, func_iPlus1;
    func_i = flux_function(u_i,d);
    func_iPlus1 = flux_function(u_iPlus1,d);
    for (int i = 0; i<u_i.size(); i++)
      {
        u_half.push_back(0.5*(u_i[i]+u_iPlus1[i])
                    - 0.5*(dt/dx)*(func_iPlus1[i]-func_i[i]));
      }
    return flux_function(u_half,d);
}
Vector<Real> flux_FORCE(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
			double dx, double dt, int d,
			std::function<Vector<Real> (const Vector<Real>&, int)> flux_function){
  Vector<Real> flux_force;
  Vector<Real> flux_lf = flux_LF(u_i, u_iPlus1, dx, dt, d, flux_function);
  Vector<Real> flux_ri = flux_RI(u_i, u_iPlus1, dx, dt, d, flux_function);
  for(int i=0; i<u_i.size(); i++)
    {
      flux_force.push_back(0.5*(flux_lf[i] + flux_ri[i]));
    }
  return flux_force;
}
////////////////////////////////////////////////////////
Real calc_mutimeGOL(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, int d, Real cfl){

  Real speed_i = get_max_seed(u_i, d);
  Real speed_iPlus1 = get_max_seed(u_iPlus1, d);
  
  return cfl/std::max(speed_i, speed_iPlus1);
}
Vector<Real> flux_GFORCE(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
			 Real dx, Real dt, int d, Real cfl,
			 std::function<Vector<Real> (const Vector<Real>&, int)> flux_function){
  Vector<Real> flux_force;
  Vector<Real> flux_lf = flux_LF(u_i, u_iPlus1, dx, dt, d, flux_function);
  Vector<Real> flux_ri = flux_RI(u_i, u_iPlus1, dx, dt, d, flux_function);

  Real omega = 1.0/(1.0+cfl);
  
  for(int i=0; i<u_i.size(); i++)
    {
      flux_force.push_back((1.0-omega)*flux_lf[i] + omega*flux_ri[i]);
    }
  return flux_force;
}
Vector<Real> flux_MUSTA1(Vector<Real>& u_i, Vector<Real>& u_iPlus1, 
			 double dx, double dt, int d, Real cfl,
			 std::function<Vector<Real> (const Vector<Real>&, int)> flux_function){

  // Evaluate fluxes FL, FR on data u_i and u_iPlus1
  Vector<Real> F_i = flux_function(u_i,d);
  Vector<Real> F_iPlus1 = flux_function(u_iPlus1,d);

  Real mutime = calc_mutimeGOL(u_i, u_iPlus1, d, cfl);

  // Compute predictor fluxes
  Vector<Real> P = flux_GFORCE(u_i, u_iPlus1, 1.0, mutime, d, cfl, flux_function);

  for (int n = 0; n<u_i.size(); n++)
    {
      u_i[n] = u_i[n] - mutime*(P[n] - F_i[n]);
      u_iPlus1[n] = u_iPlus1[n] - mutime*(F_iPlus1[n] - P[n]);
    }

    // Re-compute fluxes fluxes FL, FR on data u_i and u_iPlus1
  F_i = flux_function(u_i,d);
  F_iPlus1 = flux_function(u_iPlus1,d);

  mutime = calc_mutimeGOL(u_i, u_iPlus1, d, cfl);

  // Compute corrector fluxes
  return flux_GFORCE(u_i, u_iPlus1, 1.0, mutime, d, cfl, flux_function);
  
}
