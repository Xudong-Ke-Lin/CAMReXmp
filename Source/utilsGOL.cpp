#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#include "utils.H"

Real get_pressure_totalGOL(const Vector<Real>& u_i)
{
  // define conserved variables
  Real rho = u_i[0];
  Real momX = u_i[1];
  Real momY = u_i[2];
  Real momZ = u_i[3];
  Real E = u_i[ENER_I];
  Real J_x = u_i[MOMX_E];
  Real J_y = u_i[MOMX_E];
  Real J_z = u_i[MOMX_E];
  
  // quasineutrality n_i=n_e, rho_e=rho_i*(m_e/m_i) and m_i=1
  Real n_e = rho;
  Real m_e = 1.0/m;
  Real q_i = r_i;
  
  // define primitive variables
  Real v_x = momX/rho;
  Real v_y = momY/rho;
  Real v_z = momZ/rho;

  Real v2 = get_magnitude_squared(v_x,v_y,v_z);
  Real J2 = get_magnitude_squared(J_x,J_y,J_z);

  return (Gamma-1.0)*(E - 0.5*rho*v2 - 0.5*m_e/(q_i*q_i*n_e)*J2);
}
Real get_energy_totalGOL(const Vector<Real>& w_i)
{
  // define conserved variables
  Real rho = w_i[0];
  Real v_x = w_i[1];
  Real v_y = w_i[2];
  Real v_z = w_i[3];
  Real p = w_i[ENER_I];
  Real J_x = w_i[MOMX_E];
  Real J_y = w_i[MOMX_E];
  Real J_z = w_i[MOMX_E];
  
  // quasineutrality n_i=n_e, rho_e=rho_i*(m_e/m_i) and m_i=1
  Real n_e = rho;
  Real m_e = 1.0/m;
  Real q_i = r_i;

  Real v2 = get_magnitude_squared(v_x,v_y,v_z);
  Real J2 = get_magnitude_squared(J_x,J_y,J_z);

  return p/(Gamma-1.0) + 0.5*rho*v2 + 0.5*m_e/(q_i*q_i*n_e)*J2;
}
Real get_speedGOL(Vector<Real> u_i){
  // speed of sound for ideal gas
  Real rho = u_i[0];
  Real p = get_pressure_totalGOL(u_i);
  return std::sqrt((Gamma*p)/rho);
}
Vector<Real> get_electron_var(const Vector<Real>& u_i){

  // define conserved variables
  Real rho = u_i[0];
  Real momX = u_i[1];
  Real momY = u_i[2];
  Real momZ = u_i[3];
  Real E = u_i[ENER_I];
  Real J_x = u_i[MOMX_E];
  Real J_y = u_i[MOMY_E];
  Real J_z = u_i[MOMZ_E];  
  Real E_e = u_i[ENER_E];
  
  // quasineutrality n_i=n_e, rho_e=rho_i*(m_e/m_i) and m_i=1
  Real n_e = rho;
  Real rho_e = rho/m;
  Real q_i = r_i;
  
  // define primitive variables
  Real v_x = momX/rho;
  Real v_y = momY/rho;
  Real v_z = momZ/rho;

  Real v_x_e = v_x-J_x/(n_e*q_i);
  Real v_y_e = v_y-J_y/(n_e*q_i);
  Real v_z_e = v_z-J_z/(n_e*q_i);
  
  Vector<Real> elec_var(5,0.0);
  elec_var[0] = rho_e;
  elec_var[1] = rho_e*v_x_e;
  elec_var[2] = rho_e*v_y_e;
  elec_var[3] = rho_e*v_z_e;
  elec_var[4] = E_e;

  return elec_var;
  
}
Vector<Real> fluidGOLFlux(const Vector<Real>& u_i, int d){
  // define conserved variables
  Real rho = u_i[0];
  Real momX = u_i[1+d];
  Real momY = u_i[1+(1+d)%3];
  Real momZ = u_i[1+(2+d)%3];
  Real E = u_i[ENER_I];
  Real J_x = u_i[MOMX_E+d];
  Real J_y = u_i[MOMX_E+(1+d)%3];
  Real J_z = u_i[MOMX_E+(2+d)%3];

  // electron variables
  Vector<Real> u_e = get_electron_var(u_i);
  Real rho_e = u_e[0];
  Real momX_e = u_e[1+d];
  Real momY_e = u_e[1+(1+d)%3];
  Real momZ_e = u_e[1+(2+d)%3];
  Real E_e = u_e[4];
  
  // ion variable
  Real E_i = E - E_e;

  // quasineutrality n_i=n_e, rho_e=rho_i*(m_e/m_i) and m_i=1
  Real n_e = rho;
  Real m_e = 1.0/m;
  Real q_i = r_i;
  
  // define primitive variables
  Real v_x = momX/rho;
  Real v_y = momY/rho;
  Real v_z = momZ/rho;
  Real p = get_pressure_totalGOL(u_i);  
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure(u_e);
  //Real p_i = p - p_e;
  
  // flux function
  Vector<Real> function(u_i.size(),0.0);
  function[0] = rho*v_x;
  function[1+d] = rho*v_x*v_x + p;
  //function[1+d] = rho*v_x*v_x + p_i + p_e;
  function[1+(1+d)%3] = rho*v_x*v_y;
  function[1+(2+d)%3] = rho*v_x*v_z;
  function[ENER_I] = (E+p)*v_x + (E_e+p_e)*(-J_x/(n_e*q_i));
  //function[ENER_I] = (E_i+p_i)*v_x + (E_e+p_e)*(v_x-J_x/(n_e*q_i));
  function[MOMX_E+d] = 2.0*J_x*v_x - J_x*J_x/(n_e*q_i) -p_e*q_i/m_e;
  function[MOMX_E+(1+d)%3] = v_x*J_y + v_y*J_x - J_x*J_y/(n_e*q_i);
  function[MOMX_E+(2+d)%3] = v_x*J_z + v_z*J_x - J_x*J_z/(n_e*q_i);  
  function[ENER_E] = (E_e+p_e)*v_x_e;
  return function;
}
// flux function for the advection part when IMEX is used for the electron pressure
/*Vector<Real> fluidGOLFlux(const Vector<Real>& u_i, int d){
  // define conserved variables
  Real rho = u_i[0];
  Real momX = u_i[1+d];
  Real momY = u_i[1+(1+d)%3];
  Real momZ = u_i[1+(2+d)%3];
  Real E = u_i[ENER_I];
  Real J_x = u_i[MOMX_E+d];
  Real J_y = u_i[MOMX_E+(1+d)%3];
  Real J_z = u_i[MOMX_E+(2+d)%3];

  // electron variables
  Vector<Real> u_e = get_electron_var(u_i);
  Real rho_e = u_e[0];
  Real momX_e = u_e[1+d];
  Real momY_e = u_e[1+(1+d)%3];
  Real momZ_e = u_e[1+(2+d)%3];
  Real E_e = u_e[4];
  
  // ion variable
  Real E_i = E - E_e;

  // quasineutrality n_i=n_e, rho_e=rho_i*(m_e/m_i) and m_i=1
  Real n_e = rho;
  Real m_e = 1.0/m;
  Real q_i = r_i;
  
  // define primitive variables
  Real v_x = momX/rho;
  Real v_y = momY/rho;
  Real v_z = momZ/rho;
  Real p = get_pressure_totalGOL(u_i);  
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure(u_e);
  Real p_i = p - p_e;

  // define kinetic energy
  Real rho_k_e = 0.5*rho_e*get_magnitude_squared(v_x_e,v_y_e,v_z_e);
  
  // flux function
  Vector<Real> function(u_i.size(),0.0);
  function[0] = rho*v_x;
  function[1+d] = rho*v_x*v_x + p_i;
  function[1+(1+d)%3] = rho*v_x*v_y;
  function[1+(2+d)%3] = rho*v_x*v_z;
  function[ENER_I] = (E_i+p_i)*v_x + rho_k_e*v_x_e;
  function[MOMX_E+d] = 2.0*J_x*v_x - J_x*J_x/(n_e*q_i);
  function[MOMX_E+(1+d)%3] = v_x*J_y + v_y*J_x - J_x*J_y/(n_e*q_i);
  function[MOMX_E+(2+d)%3] = v_x*J_z + v_z*J_x - J_x*J_z/(n_e*q_i);  
  function[ENER_E] = rho_k_e*v_x_e;
  return function;
}
*/
Real get_max_speed(const Vector<Real>& u_i, int d)
{
  // define conserved variables
  Real rho = u_i[0];
  Real momX = u_i[1+d];
  Real momY = u_i[1+(1+d)%3];
  Real momZ = u_i[1+(2+d)%3];
  Real E = u_i[ENER_I];
  Real J_x = u_i[MOMX_E+d];
  Real J_y = u_i[MOMX_E+(1+d)%3];
  Real J_z = u_i[MOMX_E+(2+d)%3];
  Real E_e = u_i[ENER_E];

  // electron variables
  Vector<Real> u_e = get_electron_var(u_i);
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
  Real p = get_pressure_totalGOL(u_i);  
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure(u_e);

  Real c = get_speedGOL(u_i);
  Real c_e = get_speed(u_e);

  
  return std::max(std::abs(v_x)+c,std::abs(v_x_e)+c_e);
  // IMEX case
  //Real c_i = get_speed({rho,momX,momY,momZ,E-E_e});
  //return std::max(std::abs(v_x)+c_i,2.0*std::abs(v_x_e));
  
}
Vector<Real> flux_Rusanov(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, int d,
			  std::function<Vector<Real> (const Vector<Real>&, int)> flux_function){
  Vector<Real> flux, func_i, func_iPlus1;
  func_i = flux_function(u_i,d);
  func_iPlus1 = flux_function(u_iPlus1,d);
  Real speed_max = std::max(get_max_speed(u_i,d),get_max_speed(u_iPlus1,d));

  for (int i = 0; i<u_i.size(); i++)
    {
      flux.push_back(0.5*(func_iPlus1[i]+func_i[i]
			  - speed_max*(u_iPlus1[i]-u_i[i])));
    }
  
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
// for Rusanov
Vector<Real> monotone_fluxGOL(const Array4<Real>& arr,
			      int i, int j, int k, int iOffset, int jOffset, int kOffset,
			      int start, int len, int d){

  Vector<Real> u_iMinus1, u_i;

  for (int n = start; n<start+len; n++)
    {
      u_iMinus1.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_i.push_back(arr(i,j,k,n));      
    }
  Vector<Real> flux = flux_Rusanov(u_iMinus1,u_i,d,fluidGOLFlux);

  return flux;
}
// for FORCE
Vector<Real> monotone_fluxGOL(const Array4<Real>& arr,
			      int i, int j, int k, int iOffset, int jOffset, int kOffset,
			      int start, int len, int d, Real dx, Real dt){

  Vector<Real> u_iMinus1, u_i;

  for (int n = start; n<start+len; n++)
    {
      u_iMinus1.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_i.push_back(arr(i,j,k,n));      
    }
  Vector<Real> flux = flux_FORCE(u_iMinus1,u_i,dx,dt,d,fluidGOLFlux);
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
  //Vector<Real> flux = flux_FORCE(u_iMinus1R,u_iL,dx,dt,d,fluidGOLFlux);
  Vector<Real> flux = flux_Rusanov(u_iMinus1R,u_iL,d,fluidGOLFlux);
  return flux;
}
Vector<Real> SLIC_TVD_fluxGOL(const Array4<Real>& arr, const Array4<Real>& slopes,
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

  Vector<Real> u_iMinus1L,u_iMinus1R,u_iL,u_iR;
  for (int n = 0; n<len; n++)
    {
      u_iMinus1L.push_back(u_iMinus1[n] - 0.5*slopes_iMinus1[n]);
      u_iMinus1R.push_back(u_iMinus1[n] + 0.5*slopes_iMinus1[n]);
      u_iL.push_back(u_i[n] - 0.5*slopes_i[n]);
      u_iR.push_back(u_i[n] + 0.5*slopes_i[n]);
    }

  // Reiamnn problem left state
  Vector<Real> u_iMinus1_nPlusHalf_R = half_update_R(u_iMinus1L, u_iMinus1R, dx, dt, d, fluidGOLFlux);
  // Reiamnn problem right state
  Vector<Real> u_i_nPlusHalf_L = half_update_L(u_iL, u_iR, dx, dt, d, fluidGOLFlux);
  
  Vector<Real> flux = flux_FORCE(u_iMinus1_nPlusHalf_R,u_i_nPlusHalf_L,dx,dt,d,fluidGOLFlux);

  return flux;
}
