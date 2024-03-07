#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#include "utils.H"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real NUM_STATE = CAMReXmp::NUM_STATE;
Real NUM_STATE_FLUID = CAMReXmp::NUM_STATE_FLUID;
Real NUM_STATE_MAXWELL = CAMReXmp::NUM_STATE_MAXWELL;

//int nSlopes = CAMReXmp::nSlopes;

// read settings file                                                                                
std::string test;

Real Gamma = 0.0;
Real c_h = 0.0;

// charge-mass ratios
Real r_i = 0.0, r_e = 0.0;
// ion-electron mass ratio
Real m = 0.0;
// normalized ion Larmor radius
Real l_r = 0.0;
// normalized speed of light
Real c = 0.0;
// normalized Debye length
Real lambda_d = 0.0;
// resistivity
Real eta = 0.0;

// parameters of the divergence cleaning
Real cb = 1.0;
Real ce = 1.0;

// functions to compute the magnitudes
Real get_magnitude_squared(Real x, Real y, Real z){
    return x*x + y*y + z*z;
}
Real get_magnitude(Real x, Real y, Real z){
    return std::sqrt(get_magnitude_squared(x, y, z));
}
// calculate magnitude velocity from converved variables
Real get_v(const Vector<Real>& u_i){
    Real v_x = u_i[1]/u_i[0];
    Real v_y = u_i[2]/u_i[0];
    Real v_z = u_i[3]/u_i[0];

    return get_magnitude(v_x, v_y, v_z);
}
// calculate squared magnitude magnetic field from converved variables
Real get_B_squared(const Vector<Real>& u_i){
    Real B_x = u_i[BX];
    Real B_y = u_i[BY];
    Real B_z = u_i[BZ];

    return get_magnitude_squared(B_x, B_y, B_z);
}
// calculate magnitude magnetic field from converved variables
Real get_B(const Vector<Real>& u_i){
    return std::sqrt(get_B_squared(u_i));
}
Real get_energy(const Vector<Real>& w_i)
{
  Real rho = w_i[0];
  Real v_x = w_i[1];
  Real v_y = w_i[2];
  Real v_z = w_i[3];
  Real v = get_magnitude(v_x, v_y, v_z);
  Real p = w_i[ENER_I];

  return p/(Gamma-1.0) + 0.5*rho*v*v;
}
Real get_specific_energy(const Vector<Real>& u_i)
{
  Real rho = u_i[0];
  Real v = get_v(u_i);
  Real E = u_i[ENER_I];

  return E - 0.5*rho*v*v;
}
Real get_pressure(const Vector<Real>& u_i)
{
  Real rho = u_i[0];
  Real v = get_v(u_i);
  Real E = u_i[ENER_I];
  return (Gamma-1.0)*(E - 0.5*rho*v*v);
}
Real get_pressure_total(const Vector<Real>& u_i)
{
  Real rho = u_i[0];
  Real v = get_v(u_i);
  Real E = u_i[ENER_I];
  Real B_squared = get_B_squared(u_i);
  return (Gamma-1.0)*(E - 0.5*rho*v*v - 0.5*B_squared) + 0.5*B_squared;
}
Real get_speed(Vector<Real> u_i){
  // speed of sound for ideal gas
  Real rho = u_i[0];
  Real p = get_pressure(u_i);
  return std::sqrt((Gamma*p)/rho);
}
// Note, B_i is the magnetic field in the i direction
// e.g. for flux update in x direction, use B_x
Real get_speed_a(Real rho, Real B_i){
    // Alfvén speed for MHD
    return std::abs(B_i)/std::sqrt(rho);
}
Real get_speed_f(const Vector<Real>& u_i, int bi){
    Real rho = u_i[0];
    Real B_i = u_i[bi];
    Real B_squared = get_magnitude_squared(u_i[BX], u_i[BY], u_i[BZ]);
    // speed of sound
    Real c_s = get_speed(u_i);    
    Real c_squared = c_s*c_s + B_squared/rho;
    // fast magneto-acoustic speed for MHD
    return std::sqrt(
        0.5*(c_squared + std::sqrt(c_squared*c_squared - 4.0*c_s*c_s*B_i*B_i/rho))
    );
}
Vector<Real> func(const Vector<Real>& u_i, int d){
  // define conserved variables
  Real rho_i = u_i[0];
  Real momX_i = u_i[1+d];
  Real momY_i = u_i[2-d];
  Real momZ_i = u_i[MOMZ_I];
  Real E_i = u_i[ENER_I];
  Real rho_e = u_i[RHO_E];
  Real momX_e = u_i[MOMX_E+d];
  Real momY_e = u_i[MOMY_E-d];
  Real momZ_e = u_i[MOMZ_E];
  Real E_e = u_i[ENER_E];
  Real B_x = u_i[BX+d];
  Real B_y = u_i[BY-d];
  Real B_z = u_i[BZ];
  Real E_x = u_i[EX+d];
  Real E_y = u_i[EY-d];
  Real E_z = u_i[EZ];

  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});
  
  // flux function
  Vector<Real> function(NUM_STATE,0.0);
  function[0] = rho_i*v_x_i;
  function[1+d] = rho_i*v_x_i*v_x_i + p_i;
  function[2-d] = rho_i*v_x_i*v_y_i;
  function[3] = rho_i*v_x_i*v_z_i;
  function[ENER_I] = (E_i+p_i)*v_x_i;
  function[RHO_E] = rho_e*v_x_e;
  function[MOMX_E+d] = rho_e*v_x_e*v_x_e + p_e;
  function[MOMY_E-d] = rho_e*v_x_e*v_y_e;
  function[MOMZ_E] = rho_e*v_x_e*v_z_e;
  function[ENER_E] = (E_e+p_e)*v_x_e;
  function[BX+d] = 0.0;
  function[BY-d] = std::pow(-1.0,d+1)*E_z;
  function[BZ] = std::pow(-1.0,d)*E_y;
  function[EX+d] = 0.0;
  function[EY-d] = -std::pow(-1.0,d+1)*c*c*B_z;
  function[EZ] = -std::pow(-1.0,d)*c*c*B_y;
  return function;
}
Vector<Real> fluidFlux(const Vector<Real>& u_i, int d){
  // define conserved variables
  Real rho_i = u_i[0];
  Real momX_i = u_i[1+d];
  Real momY_i = u_i[1+(1+d)%3];
  Real momZ_i = u_i[1+(2+d)%3];
  Real E_i = u_i[ENER_I];

  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});

  // flux function
  Vector<Real> function(u_i.size(),0.0);
  function[0] = rho_i*v_x_i;
  function[1+d] = rho_i*v_x_i*v_x_i + p_i;
  function[1+(1+d)%3] = rho_i*v_x_i*v_y_i;
  function[1+(2+d)%3] = rho_i*v_x_i*v_z_i;
  function[ENER_I] = (E_i+p_i)*v_x_i;
  return function;
}
Vector<Real> MaxwellFlux(const Vector<Real>& u_i, int d){
  // define conserved variables
  Real B_x = u_i[BX_LOCAL+d];
  Real B_y = u_i[BX_LOCAL+(1+d)%3];
  Real B_z = u_i[BX_LOCAL+(2+d)%3];
  Real E_x = u_i[EX_LOCAL+d];
  Real E_y = u_i[EX_LOCAL+(1+d)%3];
  Real E_z = u_i[EX_LOCAL+(2+d)%3];
#if (AMREX_SPACEDIM >= 2)
  Real psi_b = u_i[DIVB_LOCAL];
  Real psi_e = u_i[DIVE_LOCAL];
#endif

  // flux function
  Vector<Real> function(u_i.size(),0.0);

  function[BX_LOCAL+(1+d)%3] = -E_z;
  function[BX_LOCAL+(2+d)%3] = E_y;
  function[EX_LOCAL+(1+d)%3] = c*c*B_z;
  function[EX_LOCAL+(2+d)%3] = -c*c*B_y;

#if (AMREX_SPACEDIM == 1)
  function[BX_LOCAL+d] = 0.0;
  function[EX_LOCAL+d] = 0.0;
#else
  function[BX_LOCAL+d] = cb*psi_b;
  function[EX_LOCAL+d] = ce*c*c*psi_e;
  function[DIVB_LOCAL] = cb*c*c*B_x;
  function[DIVE_LOCAL] = ce*E_x;
#endif
  
  return function;
}

Real v_exact_dir(Real B_squared, Real B, Real B_three_halves, Real EdotB, Real v_dotB, 
		 Real EcrossB_i, Real v_crossB_i, Real Bcross_EcrossB_i, Real Bcross_v_crossB_i,
		 Real B_i, Real t, Real factor){
  return EcrossB_i/B_squared + B_i*(v_dotB + factor*t*EdotB/l_r)/B_squared +
    std::cos(factor*B*t/l_r)*(Bcross_v_crossB_i - EcrossB_i)/B_squared +
    std::sin(factor*B*t/l_r)*(B_squared*v_crossB_i + Bcross_EcrossB_i)/B_three_halves;
}
Vector<Real> cross_product(Vector<Real> a, Vector<Real> b){
  return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
}
Real dot_product(Vector<Real> a, Vector<Real> b){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
Real dotProduct(Vector<Real> a, Vector<Real> b){
  Real res = 0.0;
  for (int n=0; n<a.size(); n++)
    res += a[n]*b[n];

  return res;
}
// select the method for slope limiter
// epsilon is the slope limiter
// r is the slope ratio
Real get_epsilon(Real r){
  // MIMBEE
  /*if (r <= 0){
    return 0.0;
  } else if (r>0 && r<=1.0){
    return r;
  } else{
    Real epsilon_R = 2.0/(1.0+r);
    return std::min(1.0, epsilon_R);
    }*/
  // Van-Leer
  if (r <= 0){
    return 0.0;
  } else{
    Real epsilon_R = 2.0/(1.0+r);
    Real epsilon_L = 2.0*r/(1.0+r);
    return std::min(epsilon_L, epsilon_R);
  }
  // Superbee
  /*if (r <= 0){
    return 0.0;
  } else if (r>0 && r<=0.5){
    return 2.0*r;
  } else if (r>0.5 && r<=1.0){
    return 1.0;
  } else{
    Real epsilon_R = 2.0/(1.0+r);
    return std::min(std::min(r, epsilon_R),2.0);
    }*/
}
// slope ratio for slope limiter defined in Toro's book
Real get_r(const Real& q_iMinus1, const Real& q_i, const Real& q_iPlus1){
  if (q_iPlus1==q_i){
    return 0.0;
  }
  return (q_i-q_iMinus1)/(q_iPlus1-q_i);
}

// measure of the slope in a linear reconstruction
Vector<Real> get_delta_i(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1){
  Vector<Real> delta_iMinusHalf;
  Vector<Real> delta_iPlusHalf;
  Vector<Real> delta_i;
  for (int i = 0; i<u_i.size(); i++)
    {
      delta_iMinusHalf.push_back(u_i[i] - u_iMinus1[i]);
      delta_iPlusHalf.push_back(u_iPlus1[i] - u_i[i]);
      delta_i.push_back(0.5*(delta_iPlusHalf[i]+delta_iMinusHalf[i]));
    }
  return delta_i;
}
// half time step evolution for left state in a linear reconstruction
// define this function for left and right state
Vector<Real> half_step_evolution_L(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1,
				   Vector<int> limiting_idx, Real dx, Real dt, int d,
				   std::function<Vector<Real> (const Vector<Real>&, int)> flux_function){
  
    // slope measure
    Vector<Real> delta_i = get_delta_i(u_iMinus1, u_i, u_iPlus1);

    // cell boundary extrapolated values in a linear reconstruction
    Vector<Real> u_iL;
    Vector<Real> u_iR;
    for (int i = 0; i<u_i.size(); i++)
      {
	// slope ratio 
	Real ri = get_r(u_iMinus1[limiting_idx[i]], u_i[limiting_idx[i]], u_iPlus1[limiting_idx[i]]);
	
	// slope limiter
	Real epsilon_i = get_epsilon(ri);

        u_iL.push_back(u_i[i] - 0.5*epsilon_i*delta_i[i]);
        u_iR.push_back(u_i[i] + 0.5*epsilon_i*delta_i[i]);
      }

    Vector<Real> func_iR = flux_function(u_iR, d);
    Vector<Real> func_iL = flux_function(u_iL, d);
    Vector<Real> diff;
    Vector<Real> half_step_evolution;
    for (int i = 0; i<u_i.size(); i++)
      {
	diff.push_back(func_iR[i]-func_iL[i]);
	half_step_evolution.push_back(u_iL[i] - 0.5*(dt/dx)*diff[i]);
      }
    return half_step_evolution;
}
// half time step evolution for right state in a linear reconstruction
Vector<Real> half_step_evolution_R(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1,
				   Vector<int> limiting_idx, Real dx, Real dt, int d,
				   std::function<Vector<Real> (const Vector<Real>&, int)> flux_function){
    // slope measure
    Vector<Real> delta_i = get_delta_i(u_iMinus1, u_i, u_iPlus1);

    // cell boundary extrapolated values in a linear reconstruction
    Vector<Real> u_iL;
    Vector<Real> u_iR;
    for (int i = 0; i<u_i.size(); i++)
      {
	// slope ratio
	Real ri = get_r(u_iMinus1[limiting_idx[i]], u_i[limiting_idx[i]], u_iPlus1[limiting_idx[i]]);
	
	// slope limiter
	Real epsilon_i = get_epsilon(ri);

        u_iL.push_back(u_i[i] - 0.5*epsilon_i*delta_i[i]);
        u_iR.push_back(u_i[i] + 0.5*epsilon_i*delta_i[i]);
      }

    Vector<Real> func_iR = flux_function(u_iR, d);
    Vector<Real> func_iL = flux_function(u_iL, d);
    Vector<Real> diff;
    Vector<Real> half_step_evolution;
    for (int i = 0; i<u_i.size(); i++)
      {
	diff.push_back(func_iR[i]-func_iL[i]);
	half_step_evolution.push_back(u_iR[i] - 0.5*(dt/dx)*diff[i]);
      }
    return half_step_evolution;
}

Vector<Real> flux_LF(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, double dt, int d){
    
    Vector<Real> flux, func_i, func_iPlus1;
    func_i = func(u_i,d);
    func_iPlus1 = func(u_iPlus1,d);
    for (int i = 0; i<u_i.size(); i++)
      {
        flux.push_back(0.5*(dx/dt)*(u_i[i]-u_iPlus1[i])
            + 0.5*(func_iPlus1[i]+func_i[i]));       
      }
    return flux;
}
Vector<Real> flux_RI(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, double dt, int d){

    Vector<Real> u_half, func_i, func_iPlus1;
    func_i = func(u_i,d);
    func_iPlus1 = func(u_iPlus1,d);
    for (int i = 0; i<u_i.size(); i++)
      {
        u_half.push_back(0.5*(u_i[i]+u_iPlus1[i])
                    - 0.5*(dt/dx)*(func_iPlus1[i]-func_i[i]));
      }
    return func(u_half,d);
}
Vector<Real> flux_FORCE(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
			double dx, double dt, int d){
  Vector<Real> flux_force;
  Vector<Real> flux_lf = flux_LF(u_i, u_iPlus1, dx, dt, d);
  Vector<Real> flux_ri = flux_RI(u_i, u_iPlus1, dx, dt, d);
  for(int i=0; i<u_i.size(); i++)
    {
      flux_force.push_back(0.5*(flux_lf[i] + flux_ri[i]));
    }
  return flux_force;
}
Vector<Real> flux_SLIC(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset,
                       Real dx, Real dt, int d){
  Vector<Real> u_iMinus1;
  Vector<Real> u_i;
  Vector<Real> u_iPlus1;
  Vector<Real> u_iPlus2;

  for (int n = 0; n<u_i.size(); n++)
    {
      u_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      u_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_iPlus1.push_back(arr(i,j,k,n));
      u_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));
    }

  // Reiamnn problem left state                                                  
  Vector<Real> u_i_nPlusHalf_R = half_step_evolution_R(u_iMinus1, u_i,
                                                       u_iPlus1, dx, dt, d);
  // Reiamnn problem right state                                           
  Vector<Real> u_iPlus1_nPlusHalf_L = half_step_evolution_L(u_i, u_iPlus1,
                                                            u_iPlus2, dx, dt, d);
  return flux_FORCE(u_i_nPlusHalf_R, u_iPlus1_nPlusHalf_L, dx, dt, d);
  
}
Vector<Real> flux_LF(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, Vector<Real> dt, int d){
    
    Vector<Real> flux, func_i, func_iPlus1;
    func_i = func(u_i,d);
    func_iPlus1 = func(u_iPlus1,d);
    for (int i = 0; i<u_i.size(); i++)
      {
        flux.push_back(0.5*(dx/dt[i])*(u_i[i]-u_iPlus1[i])
            + 0.5*(func_iPlus1[i]+func_i[i]));       
      }
    return flux;
}
Vector<Real> flux_RI(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, Vector<Real> dt, int d){

    Vector<Real> u_half, func_i, func_iPlus1;
    func_i = func(u_i,d);
    func_iPlus1 = func(u_iPlus1,d);
    for (int i = 0; i<u_i.size(); i++)
      {
        u_half.push_back(0.5*(u_i[i]+u_iPlus1[i])
                    - 0.5*(dt[i]/dx)*(func_iPlus1[i]-func_i[i]));
      }
    return func(u_half,d);
}
Vector<Real> flux_GFORCE(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
			 Real dx, Vector<Real> dt, int d, Real cfl){
  Vector<Real> flux_force;
  Vector<Real> flux_lf = flux_LF(u_i, u_iPlus1, dx, dt, d);
  Vector<Real> flux_ri = flux_RI(u_i, u_iPlus1, dx, dt, d);

  Real omega = 1.0/(1.0+cfl);
  
  for(int i=0; i<u_i.size(); i++)
    {
      flux_force.push_back((1.0-omega)*flux_lf[i] + omega*flux_ri[i]);
    }
  return flux_force;
}
Vector<Real> flux_MUSTA1(Vector<Real>& u_i, Vector<Real>& u_iPlus1, 
			 double dx, double dt, int d, Real cfl){

  // Evaluate fluxes FL, FR on data u_i and u_iPlus1
  Vector<Real> F_i = func(u_i,d);
  Vector<Real> F_iPlus1 = func(u_iPlus1,d);

  Real mutime = calc_mutime(u_i, u_iPlus1, d, cfl);
  Real mutime_elec = calc_mutime_elec(u_i, u_iPlus1, d, cfl);
  Real mutime_Maxwell = calc_mutime_Maxwell(u_i, u_iPlus1, d, cfl);
  Vector<Real> mutime_array = {mutime,mutime,mutime,mutime,mutime,
			       mutime_elec,mutime_elec,mutime_elec,mutime_elec,mutime_elec,
			       mutime_Maxwell,mutime_Maxwell,mutime_Maxwell,mutime_Maxwell,mutime_Maxwell,mutime_Maxwell};

  // Compute predictor fluxes
  Vector<Real> P = flux_GFORCE(u_i, u_iPlus1, 1.0, mutime_array, d, cfl);

  for (int n = 0; n<u_i.size(); n++)
    {
      u_i[n] = u_i[n] - mutime_array[n]*(P[n] - F_i[n]);
      u_iPlus1[n] = u_iPlus1[n] - mutime_array[n]*(F_iPlus1[n] - P[n]);
    }

    // Re-compute fluxes fluxes FL, FR on data u_i and u_iPlus1
  F_i = func(u_i,d);
  F_iPlus1 = func(u_iPlus1,d);

  mutime = calc_mutime(u_i, u_iPlus1, d, cfl);
  mutime_elec = calc_mutime_elec(u_i, u_iPlus1, d, cfl);
  mutime_array = {mutime,mutime,mutime,mutime,mutime,
		  mutime_elec,mutime_elec,mutime_elec,mutime_elec,mutime_elec,
		  mutime_Maxwell,mutime_Maxwell,mutime_Maxwell,mutime_Maxwell,mutime_Maxwell,mutime_Maxwell};

  // Compute corrector fluxes
  return flux_GFORCE(u_i, u_iPlus1, 1.0, mutime_array, d, cfl);
  
}
Vector<Real> flux_SLIC(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset,
                       Real dx, Real dt, int d, Real cfl){
  Vector<Real> u_iMinus1;
  Vector<Real> u_i;
  Vector<Real> u_iPlus1;
  Vector<Real> u_iPlus2;

  for (int n = 0; n<NUM_STATE; n++)
    {
      u_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      u_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_iPlus1.push_back(arr(i,j,k,n));
      u_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));
    }

  // Reiamnn problem left state                                                  
  Vector<Real> u_i_nPlusHalf_R = half_step_evolution_R(u_iMinus1, u_i,
                                                       u_iPlus1, dx, dt, d);
  // Reiamnn problem right state                                           
  Vector<Real> u_iPlus1_nPlusHalf_L = half_step_evolution_L(u_i, u_iPlus1,
                                                            u_iPlus2, dx, dt, d);
  return flux_MUSTA1(u_i_nPlusHalf_R, u_iPlus1_nPlusHalf_L, dx, dt, d, cfl);
  
}
Real calc_mutime(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, int d, Real cfl){
  Real v_x_i = u_i[1+d]/u_i[0];
  Real c_sound_i = get_speed(u_i);

  Real v_x_iPlus1 = u_iPlus1[1+d]/u_iPlus1[0];
  Real c_sound_iPlus1 = get_speed(u_iPlus1);

  return cfl/std::max(std::abs(v_x_i)+c_sound_i, std::abs(v_x_iPlus1)+c_sound_iPlus1);
}
Real calc_mutime_elec(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, int d, Real cfl){
  Real v_x_i = u_i[MOMX_E+d]/u_i[RHO_E];
  Real c_sound_i = get_speed({u_i[RHO_E],u_i[MOMX_E],u_i[MOMY_E],u_i[MOMZ_E],u_i[ENER_E]});

  Real v_x_iPlus1 = u_iPlus1[MOMX_E+d]/u_iPlus1[RHO_E];
  Real c_sound_iPlus1 = get_speed({u_iPlus1[RHO_E],u_iPlus1[MOMX_E],u_iPlus1[MOMY_E],u_iPlus1[MOMZ_E],u_iPlus1[ENER_E]});

  return cfl/std::max(std::abs(v_x_i)+c_sound_i, std::abs(v_x_iPlus1)+c_sound_iPlus1);
}
Real calc_mutime_Maxwell(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, int d, Real cfl){
  return cfl/c;
}

// S* in Riemann problem for 2D MHD
Real get_S_star(const Vector<Real>& u_L, const Vector<Real>& u_R,
                Real S_L, Real S_R, int d){
  // define variables for left and right states
  Real rho_L = u_L[0], rho_R = u_R[0];
  
  Real v_x_L = u_L[1+d]/rho_L;
  Real v_x_R = u_R[1+d]/rho_R;

  Real p_L = get_pressure(u_L), p_R = get_pressure(u_R);
    
  return (p_R-p_L+rho_L*v_x_L*(S_L-v_x_L)-rho_R*v_x_R*(S_R-v_x_R))/(rho_L*(S_L-v_x_L)-rho_R*(S_R-v_x_R));
}
// HLLC left or right states
Vector<Real> get_u_HLLC(const Vector<Real>& u_K, Real S_K, Real S_star, int d){
  // unpack the input data
  Real rho = u_K[0], E = u_K[ENER_I];
  Real v_x = u_K[1+d]/rho, v_y = u_K[1+(1+d)%3]/rho, v_z = u_K[1+(2+d)%3]/rho;
  Real p = get_pressure(u_K);

  // pre-factor in HLLC states
  Real constants = (S_K-v_x)/(S_K-S_star);
    
  // star states for conserved variables
  Vector<Real> u_HLLC(5,0.0);
  u_HLLC[0] = constants*rho;
  u_HLLC[1+d] = constants*rho*S_star;
  u_HLLC[1+(1+d)%3] = constants*rho*v_y;
  u_HLLC[1+(2+d)%3] = constants*rho*v_z;
  u_HLLC[ENER_I] = constants*rho*(E/rho + (S_star-v_x)*(S_star+p/(rho*(S_K-v_x))));

  return u_HLLC;
}
// S_R (and S_L) definition
Real get_S_K(const Vector<Real>& u_L, const Vector<Real>& u_R, int d){
  Real v_x_L = u_L[1+d]/u_L[0], v_x_R = u_R[1+d]/u_R[0];

  Real c_L = get_speed(u_L);
  Real c_R = get_speed(u_R);

  // choose fastest speeds
  return std::max(std::abs(v_x_L)+c_L, std::abs(v_x_R)+c_R);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

Real computeDerivative(Real u_iMinus1, Real u_iPlus1, Real dx){
  return (u_iPlus1-u_iMinus1)/(2.0*dx);
}
Real computeSecondDerivative(Real u_iMinus1, Real u_i, Real u_iPlus1, Real dx){
  return (u_iMinus1-2*u_i+u_iPlus1)/(dx*dx);
}

Vector<Real> half_step_evolution_L(Vector<Real> u_iMinus1, Vector<Real> u_i,
                            Vector<Real> u_iPlus1, Real dx, Real dt, int d){
    
    // slope measure
    Vector<Real> delta_i = get_delta_i(u_iMinus1, u_i, u_iPlus1);
    // slope ratio 
    // index ENER corresponds to energy (limiting on the energy)
    Real ri = get_r(u_iMinus1[ENER_I], u_i[ENER_I], u_iPlus1[ENER_I]);
    Real re = get_r(u_iMinus1[ENER_E], u_i[ENER_E], u_iPlus1[ENER_E]);
    // slope limiter
    Real epsilon_i = get_epsilon(ri);
    Real epsilon_e = get_epsilon(re);
    // cell boundary extrapolated values in a linear reconstruction
    Vector<Real> u_iL;
    Vector<Real> u_iR;
    for (int i = 0; i<=ENER_I; i++)
      {
        u_iL.push_back(u_i[i] - 0.5*epsilon_i*delta_i[i]);
        u_iR.push_back(u_i[i] + 0.5*epsilon_i*delta_i[i]);
      }
    for (int i = RHO_E; i<=ENER_E; i++)
      {
        u_iL.push_back(u_i[i] - 0.5*epsilon_e*delta_i[i]);
        u_iR.push_back(u_i[i] + 0.5*epsilon_e*delta_i[i]);
      }
    for (int i = BX; i<NUM_STATE; i++)
      {
	Real r = get_r(u_iMinus1[i], u_i[i], u_iPlus1[i]);
	Real epsilon = get_epsilon(r);
	u_iL.push_back(u_i[i] - 0.5*epsilon*delta_i[i]);
	u_iR.push_back(u_i[i] + 0.5*epsilon*delta_i[i]);
      }
    Vector<Real> func_iR = func(u_iR, d);
    Vector<Real> func_iL = func(u_iL, d);
    Vector<Real> diff;
    Vector<Real> half_step_evolution;
    for (int i = 0; i<NUM_STATE; i++)
      {
	diff.push_back(func_iR[i]-func_iL[i]);
	half_step_evolution.push_back(u_iL[i] - 0.5*(dt/dx)*diff[i]);
      }
    return half_step_evolution;
}
// half time step evolution for right state in a linear reconstruction
Vector<Real> half_step_evolution_R(Vector<Real> u_iMinus1, Vector<Real> u_i,
				   Vector<Real> u_iPlus1, Real dx, Real dt, int d){
    // slope measure
    Vector<Real> delta_i = get_delta_i(u_iMinus1, u_i, u_iPlus1);
    // slope ratio 
    // index ENER corresponds to energy (limiting on the energy)
    Real ri = get_r(u_iMinus1[ENER_I], u_i[ENER_I], u_iPlus1[ENER_I]);
    Real re = get_r(u_iMinus1[ENER_E], u_i[ENER_E], u_iPlus1[ENER_E]);
    // slope limiter                                                                      
    Real epsilon_i = get_epsilon(ri);
    Real epsilon_e = get_epsilon(re);
    // cell boundary extrapolated values in a linear reconstruction                         
    Vector<Real> u_iL;
    Vector<Real> u_iR;
    for (int i = 0; i<=ENER_I; i++)
      {
        u_iL.push_back(u_i[i] - 0.5*epsilon_i*delta_i[i]);
        u_iR.push_back(u_i[i] + 0.5*epsilon_i*delta_i[i]);
      }
    for (int i = RHO_E; i<=ENER_E; i++)
      {
        u_iL.push_back(u_i[i] - 0.5*epsilon_e*delta_i[i]);
        u_iR.push_back(u_i[i] + 0.5*epsilon_e*delta_i[i]);
      }
    for (int i = BX; i<NUM_STATE; i++)
      {
        Real r = get_r(u_iMinus1[i], u_i[i], u_iPlus1[i]);
        Real epsilon = get_epsilon(r);
        u_iL.push_back(u_i[i] - 0.5*epsilon*delta_i[i]);
        u_iR.push_back(u_i[i] + 0.5*epsilon*delta_i[i]);
      }
    Vector<Real> func_iR = func(u_iR, d);
    Vector<Real> func_iL = func(u_iL, d);
    Vector<Real> diff;
    Vector<Real> half_step_evolution;
    for (int i = 0; i<NUM_STATE; i++)
      {
	diff.push_back(func_iR[i]-func_iL[i]);
	half_step_evolution.push_back(u_iR[i] - 0.5*(dt/dx)*diff[i]);
      }
    return half_step_evolution;
}

Real get_delta_i(Real u_iMinus1, Real u_i, Real u_iPlus1){
  Real delta_iMinusHalf, delta_iPlusHalf, delta_i;
  delta_iMinusHalf = u_i - u_iMinus1;
  delta_iPlusHalf = u_iPlus1 - u_i;
  delta_i = 0.5*(delta_iPlusHalf+delta_iMinusHalf);
  return delta_i;
}
Real TVD_slope(Real u_iMinus1, Real u_i, Real u_iPlus1)
{
  // slope measure
  Real delta_i = get_delta_i(u_iMinus1, u_i, u_iPlus1);
  // slope ratio
  Real ri = get_r(u_iMinus1,u_i,u_iPlus1);
  // slope limiter
  Real epsilon_i = get_epsilon(ri);

  return epsilon_i*delta_i;
}
// second order linear reconstruction on the left
Vector<Real> TVD2_reconstruction_L(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1,
				   Vector<int> limiting_idx){
  // slope measure
    Vector<Real> delta_i = get_delta_i(u_iMinus1, u_i, u_iPlus1);

    // cell boundary extrapolated values in a linear reconstruction
    Vector<Real> u_iL;
    for (int i = 0; i<u_i.size() ; i++)    
      {
	// slope ratio 
	Real ri = get_r(u_iMinus1[limiting_idx[i]], u_i[limiting_idx[i]], u_iPlus1[limiting_idx[i]]);
	
	// slope limiter
	Real epsilon_i = get_epsilon(ri);

        u_iL.push_back(u_i[i] - 0.5*epsilon_i*delta_i[i]);
      }

    return u_iL;
}
// second order linear reconstruction on the right
Vector<Real> TVD2_reconstruction_R(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1,
				   Vector<int> limiting_idx){
  // slope measure
    Vector<Real> delta_i = get_delta_i(u_iMinus1, u_i, u_iPlus1);

    // cell boundary extrapolated values in a linear reconstruction
    Vector<Real> u_iR;
    for (int i = 0; i<u_i.size() ; i++)    
      {
	// slope ratio 
	Real ri = get_r(u_iMinus1[limiting_idx[i]], u_i[limiting_idx[i]], u_iPlus1[limiting_idx[i]]);
	
	// slope limiter
	Real epsilon_i = get_epsilon(ri);

        u_iR.push_back(u_i[i] + 0.5*epsilon_i*delta_i[i]);
      }

    return u_iR;
}
// halft time step update evolution on the left
Vector<Real> half_update_L(Vector<Real> u_iL, Vector<Real> u_iR, Real dx, Real dt, int d,
			   std::function<Vector<Real> (const Vector<Real>&, int)> flux_function)
{
  Vector<Real> func_iR = flux_function(u_iR, d);
  Vector<Real> func_iL = flux_function(u_iL, d);
  Vector<Real> diff;
  Vector<Real> half_step_evolution;
  for (int i = 0; i<u_iL.size(); i++)
      {
	diff.push_back(func_iR[i]-func_iL[i]);
	half_step_evolution.push_back(u_iL[i] - 0.5*(dt/dx)*diff[i]);
      }
  return half_step_evolution;
}
// halft time step update evolution on the right
Vector<Real> half_update_R(Vector<Real> u_iL, Vector<Real> u_iR, Real dx, Real dt, int d,
			   std::function<Vector<Real> (const Vector<Real>&, int)> flux_function)
{
  Vector<Real> func_iR = flux_function(u_iR, d);
  Vector<Real> func_iL = flux_function(u_iL, d);
  Vector<Real> diff;
  Vector<Real> half_step_evolution;
  for (int i = 0; i<u_iL.size(); i++)
      {
	diff.push_back(func_iR[i]-func_iL[i]);
	half_step_evolution.push_back(u_iR[i] - 0.5*(dt/dx)*diff[i]);
      }
  return half_step_evolution;
}

Vector<Real> MUSCL_Hancock_TVD_flux(const Array4<Real>& arr, const Array4<Real>& slopes,
				    int i, int j, int k, int iOffset, int jOffset, int kOffset,
				    int start, int len,
				    Real dx, Real dt, int d,
				    std::function<Vector<Real> (const Vector<Real>&, int)> flux_function,
				    std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver){

  Vector<Real> u_iMinus1, u_i;

  Vector<Real> slopes_iMinus1, slopes_i;

  int offset;
  if (d==0)
    offset = 0;
  else if (d==1)
    offset = NUM_STATE_FLUID;
  
  for (int n = start; n<start+len; n++)
    {
      //ui_iMinus2.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      u_iMinus1.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_i.push_back(arr(i,j,k,n));
      //ui_iPlus1.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));

      slopes_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+offset));
      slopes_i.push_back(slopes(i,j,k,n+offset)); 
    }

  // Slope limiting variable index
  //Vector<int> limiting_idx(NUM_STATE_FLUID/2, NUM_STATE_FLUID/2-1);
  
  // Cell boundary extrapolated values at the left and the right for the ion
  //Vector<Real> slopes_iMinus1_i = TVD_slopes(ui_iMinus2, ui_iMinus1, ui_i, limiting_idx);
  //Vector<Real> slopes_i_i = TVD_slopes(ui_iMinus1, ui_i, ui_iPlus1, limiting_idx);
  Vector<Real> u_iMinus1L,u_iMinus1R,u_iL,u_iR;
  for (int n = 0; n<len; n++)
    {
      u_iMinus1L.push_back(u_iMinus1[n] - 0.5*slopes_iMinus1[n]);
      u_iMinus1R.push_back(u_iMinus1[n] + 0.5*slopes_iMinus1[n]);
      u_iL.push_back(u_i[n] - 0.5*slopes_i[n]);
      u_iR.push_back(u_i[n] + 0.5*slopes_i[n]);
    }

  // Reiamnn problem left state
  Vector<Real> u_iMinus1_nPlusHalf_R = half_update_R(u_iMinus1L, u_iMinus1R, dx, dt, d, flux_function);
  // Reiamnn problem right state
  Vector<Real> u_i_nPlusHalf_L = half_update_L(u_iL, u_iR, dx, dt, d, flux_function);
  
  Vector<Real> flux = solver(u_iMinus1_nPlusHalf_R,u_i_nPlusHalf_L,d);

  return flux;
}

Vector<Real> EM_quadraticFunc(const Array4<Real>& Bc, const Array4<Real>& Ec,  int i, int j, int k, Real x, Real y, Real z, const Real* dx)
{
  // For second order, there are 18 coeff.
  // i.e. a0,ax,ay,az,axx,...,b0,bx,...,c0,cx,...,czz  
  int a0=0,ax=1,ay=2,az=3,axx=4,ayy=5,azz=6,axy=7,ayz=8,axz=9,axxx=10,axxy=11,axxz=12,axyy=13,axzz=14,axyz=15;
  int b0=16,bx=17,by=18,bz=19,bxx=20,byy=21,bzz=22,bxy=23,byz=24,bxz=25,byyy=26,bxyy=27,byyz=28,bxxy=29,byzz=30,bxyz=31;
  int c0=32,cx=33,cy=34,cz=35,cxx=36,cyy=37,czz=38,cxy=39,cyz=40,cxz=41,czzz=42,cxzz=43,cyzz=44,cxxz=45,cyyz=46,cxyz=47;
  
  Vector<Real> EM(NUM_STATE_MAXWELL,0.0);

  Real xdx = x/dx[0], ydy = y/dx[1];
     
  EM[BX_LOCAL] = Bc(i,j,k,a0) + Bc(i,j,k,ax)*xdx
    + Bc(i,j,k,ay)*ydy // + az*
    + Bc(i,j,k,axx)*(xdx*xdx - 1.0/12.0)
    + Bc(i,j,k,ayy)*(ydy*ydy - 1.0/12.0) 
    + Bc(i,j,k,axy)*xdx*ydy
    + Bc(i,j,k,axxx)*(xdx*xdx*xdx - 3.0/20.0*xdx)
    + Bc(i,j,k,axxy)*(xdx*xdx - 1.0/12.0)*ydy
    + Bc(i,j,k,axyy)*(ydy*ydy - 1.0/12.0)*xdx;
  EM[BY_LOCAL] = Bc(i,j,k,b0) + Bc(i,j,k,bx)*xdx
    + Bc(i,j,k,by)*ydy // + bz*
    + Bc(i,j,k,bxx)*(xdx*xdx - 1.0/12.0)
    + Bc(i,j,k,byy)*(ydy*ydy - 1.0/12.0)
    + Bc(i,j,k,bxy)*xdx*ydy
    + Bc(i,j,k,byyy)*(ydy*ydy*ydy - 3.0/20.0*ydy)
    + Bc(i,j,k,bxyy)*(ydy*ydy - 1.0/12.0)*xdx
    + Bc(i,j,k,bxxy)*(xdx*xdx - 1.0/12.0)*ydy;
  EM[BZ_LOCAL] = Bc(i,j,k,c0) + Bc(i,j,k,cx)*xdx
    + Bc(i,j,k,cy)*ydy
    + Bc(i,j,k,cxx)*(xdx*xdx - 1.0/12.0)
    + Bc(i,j,k,cyy)*(ydy*ydy - 1.0/12.0)
    + Bc(i,j,k,cxy)*xdx*ydy;

  EM[EX_LOCAL] = Ec(i,j,k,a0) + Ec(i,j,k,ax)*xdx
    + Ec(i,j,k,ay)*ydy // + az*
    + Ec(i,j,k,axx)*(xdx*xdx - 1.0/12.0)
    + Ec(i,j,k,ayy)*(ydy*ydy - 1.0/12.0) 
    + Ec(i,j,k,axy)*xdx*ydy
    + Ec(i,j,k,axxx)*(xdx*xdx*xdx - 3.0/20.0*xdx)
    + Ec(i,j,k,axxy)*(xdx*xdx - 1.0/12.0)*ydy
    + Ec(i,j,k,axyy)*(ydy*ydy - 1.0/12.0)*xdx;
  EM[EY_LOCAL] = Ec(i,j,k,b0) + Ec(i,j,k,bx)*xdx
    + Ec(i,j,k,by)*ydy // + bz*
    + Ec(i,j,k,bxx)*(xdx*xdx - 1.0/12.0)
    + Ec(i,j,k,byy)*(ydy*ydy - 1.0/12.0)
    + Ec(i,j,k,bxy)*xdx*ydy
    + Ec(i,j,k,byyy)*(ydy*ydy*ydy - 3.0/20.0*ydy)
    + Ec(i,j,k,bxyy)*(ydy*ydy - 1.0/12.0)*xdx
    + Ec(i,j,k,bxxy)*(xdx*xdx - 1.0/12.0)*ydy;
  EM[EZ_LOCAL] = Ec(i,j,k,c0) + Ec(i,j,k,cx)*xdx
    + Ec(i,j,k,cy)*ydy
    + Ec(i,j,k,cxx)*(xdx*xdx - 1.0/12.0)
    + Ec(i,j,k,cyy)*(ydy*ydy - 1.0/12.0)
    + Ec(i,j,k,cxy)*xdx*ydy;
  
  return EM;
  
}
Vector<Real> EM_linearFunc(const Array4<Real>& Bc, const Array4<Real>& Ec,  int i, int j, int k, Real x, Real y, Real z, const Real* dx)
{
  // For second order, there are 27 coeff.
  // i.e. a0,ax,ay,az,axx,...,b0,bx,...,c0,cx,...,czz  
  int a0=0,ax=1,ay=2,az=3,axx=4,axy=5,axz=6;
  int b0=7,bx=8,by=9,bz=10,byy=11,bxy=12,byz=13;
  int c0=14,cx=15,cy=16,cz=17,czz=18,cxz=19,cyz=20;
  
  Vector<Real> EM(NUM_STATE_MAXWELL,0.0);

  Real xdx = x/dx[0], ydy = y/dx[1];
  
  EM[BX_LOCAL] = Bc(i,j,k,a0) + Bc(i,j,k,ax)*xdx
    + Bc(i,j,k,ay)*ydy // + az*
    + Bc(i,j,k,axx)*(xdx*xdx - 1.0/12.0)
    + Bc(i,j,k,axy)*xdx*ydy; // + axz*
  EM[BY_LOCAL] = Bc(i,j,k,b0) + Bc(i,j,k,bx)*xdx
    + Bc(i,j,k,by)*ydy // + bz*
    + Bc(i,j,k,byy)*(ydy*ydy - 1.0/12.0)
    + Bc(i,j,k,bxy)*xdx*ydy; // + bxz*
  EM[BZ_LOCAL] = Bc(i,j,k,c0) + Bc(i,j,k,cx)*xdx
    + Bc(i,j,k,cy)*ydy; // + cz* + cxz* + cyz* + czz		  

  EM[EX_LOCAL] = Ec(i,j,k,a0) + Ec(i,j,k,ax)*xdx
    + Ec(i,j,k,ay)*ydy // + az*
    + Ec(i,j,k,axx)*(xdx*xdx - 1.0/12.0)
    + Ec(i,j,k,axy)*xdx*ydy; // + axz*
  EM[EY_LOCAL] = Ec(i,j,k,b0) + Ec(i,j,k,bx)*xdx
    + Ec(i,j,k,by)*ydy // + bz*
    + Ec(i,j,k,byy)*(ydy*ydy - 1.0/12.0)
    + Ec(i,j,k,bxy)*xdx*ydy; // + bxz*
  EM[EZ_LOCAL] = Ec(i,j,k,c0) + Ec(i,j,k,cx)*xdx
    + Ec(i,j,k,cy)*ydy; // + cz* + cxz* + cyz* + czz		  

  return EM;
  
}
Vector<Real> HLLC(Vector<Real> uL, Vector<Real> uR, int d)
{
  // speeds for approximate Riemann problem
  Real S_R = get_S_K(uL, uR, d);
  Real S_L = -S_R;
  Real S_star = get_S_star(uL, uR, S_L, S_R, d);  

  // flux depending on the speeds, defined in slides or Toro's book
  Vector<Real> flux(uL.size(),0.0);
  // HLLC
  if (S_L>=0){
    Vector<Real> function = fluidFlux(uL, d);
    for (int n = 0; n<=ENER_I; n++)
      {       
        flux[n] = function[n];
      }
  } else if (S_L<=0 && S_star>=0){
    Vector<Real> u_HLLC_L= get_u_HLLC(uL, S_L, S_star, d);
    Vector<Real> function = fluidFlux(uL, d);
    Vector<Real> diff;
    for (int n = 0; n<=ENER_I; n++)
      {
        diff.push_back(S_L*(u_HLLC_L[n]-uL[n]));
	flux[n] = function[n]+diff[n];
      }
  } else if (S_star<=0 && S_R>=0){
    Vector<Real> u_HLLC_R = get_u_HLLC(uR, S_R, S_star, d);
    Vector<Real> function = fluidFlux(uR, d);
    Vector<Real> diff;
    for (int n = 0; n<=ENER_I; n++)
      {
        diff.push_back(S_R*(u_HLLC_R[n]-uR[n]));
	flux[n] = function[n]+diff[n];
      }
  } else {    
    Vector<Real> function = fluidFlux(uR, d);
    for (int n = 0; n<=ENER_I; n++)
      {
        flux[n] = function[n];
      }
  }

  return flux;
}
Vector<Real> RusanovEuler(Vector<Real> uL, Vector<Real> uR, int d)
{
  Vector<Real> flux, func_i, func_iPlus1;
  func_i = fluidFlux(uL,d);
  func_iPlus1 = fluidFlux(uR,d);
  Real speed_max = std::max(std::abs(uL[1+d]/uL[0])+get_speed(uL),
			    std::abs(uR[1+d]/uR[0])+get_speed(uR));

  for (int i = 0; i<uL.size(); i++)
    {
      flux.push_back(0.5*(func_iPlus1[i]+func_i[i]
			  - speed_max*(uR[i]-uL[i])));
    }
  
  return flux;

}
Vector<Real> RankineHugoniot(Vector<Real> uL, Vector<Real> uR, int d)
{
  Vector<Real> flux(NUM_STATE_MAXWELL,0.0);
  
  Real B_y_star = 0.5*(uR[BX_LOCAL+(1+d)%3]+uL[BX_LOCAL+(1+d)%3])
    +0.5/c*(uR[EX_LOCAL+(2+d)%3]-uL[EX_LOCAL+(2+d)%3]);
  Real B_z_star = 0.5*(uR[BX_LOCAL+(2+d)%3]+uL[BX_LOCAL+(2+d)%3])
    -0.5/c*(uR[EX_LOCAL+(1+d)%3]-uL[EX_LOCAL+(1+d)%3]);
  Real E_y_star = 0.5*(uR[EX_LOCAL+(1+d)%3]+uL[EX_LOCAL+(1+d)%3])
    -0.5*c*(uR[BX_LOCAL+(2+d)%3]-uL[BX_LOCAL+(2+d)%3]);
  Real E_z_star = 0.5*(uR[EX_LOCAL+(2+d)%3]+uL[EX_LOCAL+(2+d)%3])
    +0.5*c*(uR[BX_LOCAL+(1+d)%3]-uL[BX_LOCAL+(1+d)%3]);  

#if (AMREX_SPACEDIM != 1)  
  Real B_x_star = 0.5*(uR[BX_LOCAL+d]+uL[BX_LOCAL+d])-0.5/c*(uR[DIVB_LOCAL]-uL[DIVB_LOCAL]);
  Real E_x_star = 0.5*(uR[EX_LOCAL+d]+uL[EX_LOCAL+d])-0.5*c*(uR[DIVE_LOCAL]-uL[DIVE_LOCAL]);
  Real psi_b_star = 0.5*(uR[DIVB_LOCAL]+uL[DIVB_LOCAL])-0.5*c*(uR[BX_LOCAL+d]-uL[BX_LOCAL+d]);
  Real psi_e_star = 0.5*(uR[DIVE_LOCAL]+uL[DIVE_LOCAL])-0.5/c*(uR[EX_LOCAL+d]-uL[EX_LOCAL+d]);
#endif
  
  // EM HLLC states
  flux[BX_LOCAL+(1+d)%3] = -E_z_star;
  flux[BX_LOCAL+(2+d)%3] = E_y_star;
  flux[EX_LOCAL+(1+d)%3] = c*c*B_z_star;
  flux[EX_LOCAL+(2+d)%3] =  -c*c*B_y_star;

#if (AMREX_SPACEDIM == 1)
  flux[BX_LOCAL+d] = 0.0;
  flux[EX_LOCAL+d] = 0.0;  
#else
  flux[BX_LOCAL+d] = cb*psi_b_star;
  flux[EX_LOCAL+d] = ce*c*c*psi_e_star;
  flux[DIVB_LOCAL] = cb*c*c*B_x_star;
  flux[DIVE_LOCAL] = ce*E_x_star; 
#endif
  
  return flux;
}
Vector<Real> MUSCL_Hancock_WENO_flux(const Array4<Real>& arr, const Array4<Real>& slopes,
				     int i, int j, int k, int iOffset, int jOffset, int kOffset,
				     int start, int len,
				     Real dx, Real dt, int d,
				     std::function<Vector<Real> (const Vector<Real>&, int)> flux_function,
				     std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver){
  Vector<Real> u_iMinus1;
  Vector<Real> u_i;

  Vector<Real> ux_iMinus1,uxx_iMinus1,ux_i,uxx_i;
  
#if (AMREX_SPACEDIM >= 2)
  Vector<Real> uy_iMinus1,uyy_iMinus1,uy_i,uyy_i;
  Vector<Real> uxy_iMinus1,uxy_i;
#endif

  int xOffset,yOffset;
  if (d==0)
    xOffset = 0, yOffset = 2*NUM_STATE_FLUID;
  else if (d==1)
    xOffset = 2*NUM_STATE_FLUID, yOffset = 0;
  
  for (int n = start; n<start+len; n++)
    {
      u_iMinus1.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_i.push_back(arr(i,j,k,n));

      // slopes in x-direction
      ux_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+xOffset));
      uxx_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+xOffset+NUM_STATE_FLUID));
      ux_i.push_back(slopes(i,j,k,n+xOffset));
      uxx_i.push_back(slopes(i,j,k,n+xOffset+NUM_STATE_FLUID));

#if (AMREX_SPACEDIM >= 2)
      uy_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+yOffset));
      uyy_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+yOffset+NUM_STATE_FLUID));
      uy_i.push_back(slopes(i,j,k,n+yOffset));
      uyy_i.push_back(slopes(i,j,k,n+yOffset+NUM_STATE_FLUID));

      uxy_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+4*NUM_STATE_FLUID));
      uxy_i.push_back(slopes(i,j,k,n+4*NUM_STATE_FLUID));
#endif      
    }

  
  Real x = 0.0, Lx, Lxx;
  // Cell boundary extrapolated values at the left and the right
  Vector<Real> u_iMinus1L,u_iMinus1R,u_iL,u_iR;

#if (AMREX_SPACEDIM >= 2)
  Real y = 0.0, Ly, Lyy;
  Vector<Real> u_iMinus1L2,u_iMinus1R2,u_iL2,u_iR2;
#endif
  
  for (int n = 0; n<u_i.size(); n++)
    {
      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1L.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			   + uxx_iMinus1[n]*Lxx);
      u_iL.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx);

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			   + uxx_iMinus1[n]*Lxx);
      u_iR.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx);

#if (AMREX_SPACEDIM >= 2)
      y = 1.0/(2.0*std::sqrt(3.0));
      Ly = y, Lyy = y*y - 1.0/12.0;

      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1L[n] +=  uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly;
      u_iL[n] +=  uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly;

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R[n] +=  uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly;
      u_iR[n] +=  uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly;

      // second quadrature point
      y = -1.0/(2.0*std::sqrt(3.0));
      Ly = y, Lyy = y*y - 1.0/12.0;     
      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1L2.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			     + uxx_iMinus1[n]*Lxx
			     + uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly);
      u_iL2.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx
		      + uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly);

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R2.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			    + uxx_iMinus1[n]*Lxx
			    + uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly);
      u_iR2.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx
		      + uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly);
      
#endif
    }

  // Reiamnn problem left state
  Vector<Real> u_iMinus1_nPlusHalf_R = half_update_R(u_iMinus1L, u_iMinus1R, dx, dt, d, flux_function);
  // Reiamnn problem right state
  Vector<Real> u_i_nPlusHalf_L = half_update_L(u_iL, u_iR, dx, dt, d, flux_function);

  Vector<Real> flux = solver(u_iMinus1_nPlusHalf_R,u_i_nPlusHalf_L,d);

#if (AMREX_SPACEDIM >= 2)
  Vector<Real> u_iMinus1_nPlusHalf_R2 = half_update_R(u_iMinus1L2, u_iMinus1R2, dx, dt, d, flux_function);
  Vector<Real> u_i_nPlusHalf_L2 = half_update_L(u_iL2, u_iR2, dx, dt, d, flux_function);

  Vector<Real> flux2 = solver(u_iMinus1_nPlusHalf_R2,u_i_nPlusHalf_L2,d);
  for (int n = 0; n<u_i.size(); n++)
    {
      flux[n] += flux2[n];
      flux[n] *= 0.5;
    }
#endif  
  
  return flux;
}
Vector<Real> flux_LLF(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, int d){
    
    Vector<Real> flux, func_i, func_iPlus1;
    func_i = fluidFlux(u_i,d);
    func_iPlus1 = fluidFlux(u_iPlus1,d);
    Real speed_i, speed_iPlus1, speed_max;
    speed_i = get_speed(u_i);
    speed_iPlus1 = get_speed(u_iPlus1);
    speed_max = std::max(speed_i,speed_iPlus1);

    for (int i = 0; i<u_i.size(); i++)
      {
        flux.push_back(0.5*(func_iPlus1[i]+func_i[i]
			    - speed_max*(u_iPlus1[i]-u_i[i])));
      }
    return flux;
}
Vector<Real> monotone_flux(const Array4<Real>& arr,
 			   int i, int j, int k, int iOffset, int jOffset, int kOffset,
 			   int start, int len, int d,
 			   std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver){

  Vector<Real> u_iMinus1, u_i;

  for (int n = start; n<start+len; n++)
    {
      u_iMinus1.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_i.push_back(arr(i,j,k,n));
    }

  // flux depending on the speeds, defined in slides or Toro's book
  Vector<Real> flux = solver(u_iMinus1,u_i,d);

  return flux;
}
Vector<Real> TVD_flux(const Array4<Real>& arr, const Array4<Real>& slopes,
		      int i, int j, int k, int iOffset, int jOffset, int kOffset,
		      int start, int len, int startSlope, int d,
		      std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver){

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
  Vector<Real> flux = solver(u_iMinus1R,u_iL,d);

  return flux;
}
Vector<Real> WENO_flux(const Array4<Real>& arr, const Array4<Real>& slopes,
		       int i, int j, int k, int iOffset, int jOffset, int kOffset,
		       int start, int len, int d,
		       std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver){
  Vector<Real> u_iMinus1;
  Vector<Real> u_i;
  
  Vector<Real> ux_iMinus1,uxx_iMinus1,ux_i,uxx_i;
  
#if (AMREX_SPACEDIM >= 2)
  Vector<Real> uy_iMinus1,uyy_iMinus1,uy_i,uyy_i;
  Vector<Real> uxy_iMinus1,uxy_i;
#endif

  int xOffset,yOffset;
  if (d==0)
    xOffset = 0, yOffset = 2*NUM_STATE_FLUID;
  else if (d==1)
    xOffset = 2*NUM_STATE_FLUID, yOffset = 0;
  
  for (int n = start; n<start+len; n++)
    {
      u_iMinus1.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_i.push_back(arr(i,j,k,n));      
      
      // slopes in x-direction
      ux_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+xOffset));
      uxx_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+xOffset+NUM_STATE_FLUID));
      ux_i.push_back(slopes(i,j,k,n+xOffset));
      uxx_i.push_back(slopes(i,j,k,n+xOffset+NUM_STATE_FLUID));

#if (AMREX_SPACEDIM >= 2)
      uy_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+yOffset));
      uyy_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+yOffset+NUM_STATE_FLUID));
      uy_i.push_back(slopes(i,j,k,n+yOffset));
      uyy_i.push_back(slopes(i,j,k,n+yOffset+NUM_STATE_FLUID));

      uxy_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+4*NUM_STATE_FLUID));
      uxy_i.push_back(slopes(i,j,k,n+4*NUM_STATE_FLUID));
#endif      
    }

  Real x = 0.0, Lx, Lxx;
  // Cell boundary extrapolated values at the left and the right
  Vector<Real> u_iMinus1R,u_iL;
#if (AMREX_SPACEDIM >= 2)
  Real y = 0.0, Ly, Lyy;
  Vector<Real> u_iMinus1R2,u_iL2;
  //Vector<Real> u_iMinus1R3,u_iL3;
#endif
  
  for (int n = 0; n<u_i.size(); n++)
    {
      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iL.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx);
  
      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			   + uxx_iMinus1[n]*Lxx);

#if (AMREX_SPACEDIM >= 2)
      y = 1.0/(2.0*std::sqrt(3.0));
      Ly = y, Lyy = y*y - 1.0/12.0;

      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iL[n] +=  uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly;

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R[n] +=  uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly;

      // second quadrature point
      y = -1.0/(2.0*std::sqrt(3.0));
      Ly = y, Lyy = y*y - 1.0/12.0;     
      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iL2.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx
		      + uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly);

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R2.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			    + uxx_iMinus1[n]*Lxx
			    + uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly); 
      /*
      // Simpson rule
      y = -0.5;
      Ly = y, Lyy = y*y - 1.0/12.0;

      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iL[n] +=  uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly;

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R[n] +=  uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly;

      // second quadrature point
      y = 0.0;
      Ly = y, Lyy = y*y - 1.0/12.0;     
      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iL2.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx
		      + uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly);

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R2.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			    + uxx_iMinus1[n]*Lxx
			    + uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly);

      // third quadrature point
      y = 0.5;
      Ly = y, Lyy = y*y - 1.0/12.0;     
      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iL3.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx
		      + uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly);

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R3.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			    + uxx_iMinus1[n]*Lxx
			    + uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly);
      */
     
#endif
    }

  Vector<Real> flux = solver(u_iMinus1R,u_iL,d);
  //Vector<Real> flux = flux_LLF(u_iMinus1R,u_iL,d);
  
#if (AMREX_SPACEDIM >= 2)  
  Vector<Real> flux2 = solver(u_iMinus1R2,u_iL2,d);
  //Vector<Real> flux3 = solver(u_iMinus1R3,u_iL3,d);
  for (int n = 0; n<u_i.size(); n++)
    {
      flux[n] += flux2[n];
      flux[n] *= 0.5;
      //flux[n] = flux[n]/6.0 + 2.0*flux2[n]/3.0 + flux3[n]/6.0;
    }
#endif  

  /*if (i==160 && j==60 && start==5)
    std::cout << "After reconstruction " << get_specific_energy(u_iMinus1R) << " " << get_specific_energy(u_iL) << " " << get_specific_energy(u_iMinus1R2) << " " << get_specific_energy(u_iL2) << std::endl;
  */
  return flux;
}
Vector<Real> WENO_flux_flat(const Array4<Real>& arr, const Array4<Real>& slopes, const Array4<Real>& tau,
			    int i, int j, int k, int iOffset, int jOffset, int kOffset,
			    int start, int len, int d,
			    std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver){
  Vector<Real> u_iMinus1;
  Vector<Real> u_i;

  Vector<Real> ux_iMinus1,uxx_iMinus1,ux_i,uxx_i;
  
#if (AMREX_SPACEDIM >= 2)
  Vector<Real> uy_iMinus1,uyy_iMinus1,uy_i,uyy_i;
  Vector<Real> uxy_iMinus1,uxy_i;
#endif

  int xOffset,yOffset;
  if (d==0)
    xOffset = 0, yOffset = 2*NUM_STATE_FLUID;
  else if (d==1)
    xOffset = 2*NUM_STATE_FLUID, yOffset = 0;
  
  for (int n = start; n<start+len; n++)
    {
      u_iMinus1.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_i.push_back(arr(i,j,k,n));

      // slopes in x-direction
      ux_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+xOffset));
      uxx_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+xOffset+NUM_STATE_FLUID));
      ux_i.push_back(slopes(i,j,k,n+xOffset));
      uxx_i.push_back(slopes(i,j,k,n+xOffset+NUM_STATE_FLUID));

#if (AMREX_SPACEDIM >= 2)
      uy_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+yOffset));
      uyy_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+yOffset+NUM_STATE_FLUID));
      uy_i.push_back(slopes(i,j,k,n+yOffset));
      uyy_i.push_back(slopes(i,j,k,n+yOffset+NUM_STATE_FLUID));

      uxy_iMinus1.push_back(slopes(i-iOffset,j-jOffset,k-kOffset,n+4*NUM_STATE_FLUID));
      uxy_i.push_back(slopes(i,j,k,n+4*NUM_STATE_FLUID));
#endif      
    }

  
  Real x = 0.0, Lx, Lxx;
  // Cell boundary extrapolated values at the left and the right
  Vector<Real> u_iMinus1R,u_iL;

#if (AMREX_SPACEDIM >= 2)
  Real y = 0.0, Ly, Lyy;
  Vector<Real> u_iMinus1R2,u_iL2;
#endif
  
  for (int n = 0; n<u_i.size(); n++)
    {
      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iL.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx);
  
      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			   + uxx_iMinus1[n]*Lxx);
      
#if (AMREX_SPACEDIM >= 2)
      y = 1.0/(2.0*std::sqrt(3.0));
      Ly = y, Lyy = y*y - 1.0/12.0;

      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iL[n] +=  uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly;

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R[n] +=  uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly;

      // second quadrature point
      y = -1.0/(2.0*std::sqrt(3.0));
      Ly = y, Lyy = y*y - 1.0/12.0;     
      x = -0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iL2.push_back(u_i[n] + ux_i[n]*Lx + uxx_i[n]*Lxx
		      + uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly);

      x = 0.5;
      Lx = x, Lxx = x*x - 1.0/12.0;
      u_iMinus1R2.push_back(u_iMinus1[n] + ux_iMinus1[n]*Lx
			    + uxx_iMinus1[n]*Lxx
			    + uy_iMinus1[n]*Ly + uyy_iMinus1[n]*Lyy + uxy_iMinus1[n]*Lx*Ly);
      
#endif

      for (int tau_n=0; tau_n<2; tau_n++)
	{
	  if (tau(i,j,k,tau_n)<1.0)
	    {
	      //std::cout << "Within function iL " << i << " " << tau_n << " " << tau(i,j,k,tau_n) << " " << u_iL[n] << " ";
	      u_iL[n] = (1.0-tau(i,j,k,tau_n))*u_i[n] + tau(i,j,k,tau_n)*u_iL[n];
#if (AMREX_SPACEDIM >= 2)
	      u_iL2[n] = (1.0-tau(i,j,k,tau_n))*u_i[n] + tau(i,j,k,tau_n)*u_iL2[n];
#endif
	      //std::cout <<  u_iL[n] << std::endl; 
	    }
	  if (tau(i-1,j,k,tau_n)<1.0)
	    {
	      //std::cout << "Within function iMinus1R " << i << " " << tau_n << " " << tau(i-1,j,k,tau_n) << " " << u_iMinus1R[n] << " ";
	      u_iMinus1R[n] = (1.0-tau(i-1,j,k,tau_n))*u_iMinus1[n] + tau(i-1,j,k,tau_n)*u_iMinus1R[n];
#if (AMREX_SPACEDIM >= 2)
	      u_iMinus1R2[n] = (1.0-tau(i-1,j,k,tau_n))*u_iMinus1[n] + tau(i-1,j,k,tau_n)*u_iMinus1R2[n];
#endif
	      //std::cout <<  u_iMinus1R[n] << std::endl;
	    }
	}            
    }
  /*
  if (tau(i,j,k,1)<1.0)
    {
      std::cout << "Within function iL pressure " << i << " " << get_pressure(u_iL) << std::endl;
    }
  if (tau(i-1,j,k,1)<1.0)
    {
      std::cout << "Within function iMinus1R pressure " << i << " " << get_pressure(u_iMinus1R) << std::endl;
    }
  */

  Vector<Real> flux = solver(u_iMinus1R,u_iL,d);
  //Vector<Real> flux = flux_LLF(u_iMinus1R,u_iL,d);

#if (AMREX_SPACEDIM >= 2)

  Vector<Real> flux2 = solver(u_iMinus1R2,u_iL2,d);
  for (int n = 0; n<u_i.size(); n++)
    {
      flux[n] += flux2[n];
      flux[n] *= 0.5;
    }
#endif  
  
  return flux;
}
Vector<Real> get_data_zone(const Array4<Real>& arr, int i, int j, int k, int start, int length)
{
  Vector<Real> data;
  for (int n=start; n<start+length; n++)
    data.push_back(arr(i,j,k,n));
  return data;
}
Vector<Real> get_data_stencil(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset, int n)
{
  Vector<Real> data;
  for (int dataiOffset = -iOffset; dataiOffset <= iOffset; dataiOffset++)
    {
      for (int datajOffset = -jOffset; datajOffset <= jOffset; datajOffset++)
	{
	  for (int datakOffset = -kOffset; datakOffset <= kOffset; datakOffset++)
	    {
	      data.push_back(arr(i+dataiOffset,
				 j+datajOffset,
				 k+datakOffset,n));
	    }
	}
    }
  return data;
}
Real TVD_slope(Vector<Real> u, Vector<Real> limiter)
{
  Real u_iMinus1 = u[0];
  Real u_i = u[1];
  Real u_iPlus1 = u[2];

  Real lim_iMinus1 = limiter[0];
  Real lim_i = limiter[1];
  Real lim_iPlus1 = limiter[2];

  // slope measure
  Real delta_i = get_delta_i(u_iMinus1, u_i, u_iPlus1);
  // slope ratio
  Real ri = get_r(lim_iMinus1,lim_i,lim_iPlus1);
  // slope limiter
  Real epsilon_i = get_epsilon(ri);

  return epsilon_i*delta_i;
}

std::array<Real, 2> WENO3_slope(Vector<Real> u)
{
  Real u_iMinus2 = u[0];
  Real u_iMinus1 = u[1];
  Real u_i = u[2];
  Real u_iPlus1 = u[3];
  Real u_iPlus2 = u[4];

  int L=0,C=1,R=2;
  
  Vector<Real> ux(3, 0.0);
  Vector<Real> uxx(3, 0.0);
  ux[L] = -2.0*u_iMinus1 + 0.5*u_iMinus2 + 1.5*u_i;
  uxx[L] = 0.5*u_iMinus2 - u_iMinus1 + 0.5*u_i;
  ux[C] = 0.5*(u_iPlus1-u_iMinus1);
  uxx[C] = 0.5*u_iMinus1 - u_i + 0.5*u_iPlus1;
  ux[R] = -1.5*u_i + 2.0*u_iPlus1 - 0.5*u_iPlus2;
  uxx[R] = 0.5*u_i  - u_iPlus1 + 0.5*u_iPlus2;

  // smothness indicator
  Vector<Real> IS(ux.size(), 0.0);
  for (int n=0; n<ux.size(); n++)
    IS[n] = ux[n]*ux[n] + 13.0/3.0*uxx[n]*uxx[n];

  Real e = 1e-12;
  // linear weights
  Vector<Real> d(ux.size(), 0.0);
  d[L] = 1.0, d[C] = 100.0, d[R] = 1.0;
  //d[0] = 0.1, d[1] = 0.6, d[2] = 0.3;
  
  Vector<Real> alpha(ux.size(), 0.0);
  for (int n=0; n<ux.size(); n++)
    alpha[n] = d[n]/(std::pow(e+IS[n],4));
  /*
  Real alpha1_i = d1/((e+IS1)*(e+IS1));
  Real alpha2_i = d2/((e+IS2)*(e+IS2));
  Real alpha3_i = d3/((e+IS3)*(e+IS3));
  */
  Real alpha_sum = 0.0;
  for (auto& n : alpha)
    alpha_sum += n;
  //Real alpha_sum = alpha1_i + alpha2_i + alpha3_i;

  // non-linear weights
  Vector<Real> omega(ux.size(), 0.0);
  for (int n=0; n<ux.size(); n++)
    omega[n] = alpha[n]/alpha_sum;
  /*
  Real omega1_i = alpha1_i/alpha_i;
  Real omega2_i = alpha2_i/alpha_i;
  Real omega3_i = alpha3_i/alpha_i;
  */
  /*for (int n=0; n<5; n++)
    std::cout << u[n] << " ";
    std::cout << omega[L] << " " << omega[C] << " " << omega[R] << " " << omega[0]+omega[1]+omega[2] << std::endl;*/
  // Slopes
  //Real ux = omega1_i*ux1 + omega2_i*ux2 + omega3_i*ux3;
  //Real uxx = omega1_i*uxx1 + omega2_i*uxx2 + omega3_i*uxx3;
  Real ux_final = 0.0, uxx_final = 0.0;
  for (int n=0; n<ux.size(); n++)
    {
      ux_final += omega[n]*ux[n];
      uxx_final += omega[n]*uxx[n];
    }
  std::array<Real, 2> slopes = {ux_final, uxx_final};

  return slopes;
}
Real WENO3_slopeCross(Vector<Real> u, Vector<Real> slopes)
{
  Real u_iMinus1jMinus1 = u[0];
  Real u_iMinus1jPlus1 = u[2];
  Real u_ij = u[4];
  Real u_iPlus1jMinus1 = u[6];
  Real u_iPlus1jPlus1 = u[8];

  Real ux = slopes[0];
  Real uxx = slopes[1];
  Real uy = slopes[2];
  Real uyy = slopes[3];

  Vector<Real> uxy(4, 0.0);
  uxy[0] = u_iPlus1jPlus1 - u_ij - ux - uy - uxx - uyy;
  uxy[1] = -u_iPlus1jMinus1 + u_ij + ux - uy + uxx + uyy;
  uxy[2] = -u_iMinus1jPlus1 + u_ij - ux + uy + uxx + uyy;
  uxy[3] = u_iMinus1jMinus1 - u_ij + ux + uy - uxx - uyy;
  
  /*
  Real uxy1 = u_iPlus1jPlus1 - u_ij - ux - uy - uxx - uyy;
  Real uxy2 = -u_iPlus1jMinus1 + u_ij + ux - uy + uxx + uyy;
  Real uxy3 = -u_iMinus1jMinus1 + u_ij - ux + uy + uxx + uyy;
  Real uxy4 = u_iMinus1jMinus1 - u_ij + ux + uy - uxx - uyy;
  */
  // smothness indicator
  Vector<Real> IS(uxy.size(), 0.0);
  for (int n=0; n<uxy.size(); n++)
    IS[n] = 4.0*uxx*uxx + 4.0*uyy*uyy + uxy[n]*uxy[n];
  /*
  Real IS1 = 4.0*uxx*uxx + 4.0*uyy*uyy + uxy1*uxy1;
  Real IS2 = 4.0*uxx*uxx + 4.0*uyy*uyy + uxy2*uxy2;
  Real IS3 = 4.0*uxx*uxx + 4.0*uyy*uyy + uxy3*uxy3;
  Real IS4 = 4.0*uxx*uxx + 4.0*uyy*uyy + uxy4*uxy4;
  */

  Real d = 0.25;
  Real e = 1e-12;

  Vector<Real> alpha(uxy.size(), 0.0);
  for (int n=0; n<uxy.size(); n++)
    alpha[n] = d/(std::pow(e+IS[n],4));
  /*
  Real alpha1_i = d/((e+IS1)*(e+IS1));
  Real alpha2_i = d/((e+IS2)*(e+IS2));
  Real alpha3_i = d/((e+IS3)*(e+IS3));
  Real alpha4_i = d/((e+IS4)*(e+IS4));
  */
  //Real alpha_i = alpha1_i + alpha2_i + alpha3_i + alpha4_i;
  Real alpha_sum = 0.0;
  for (auto& n : alpha)
    alpha_sum += n;

  // non-linear weights
  Vector<Real> omega(uxy.size(), 0.0);
  for (int n=0; n<uxy.size(); n++)
    omega[n] = alpha[n]/alpha_sum;
  /*
  Real omega1_i = alpha1_i/alpha_i;
  Real omega2_i = alpha2_i/alpha_i;
  Real omega3_i = alpha3_i/alpha_i;
  Real omega4_i = alpha4_i/alpha_i;
  */
  // Cell boundary extrapolated values at the left and the right
  //Real uxy = omega1_i*uxy1 + omega2_i*uxy2 + omega3_i*uxy3 + omega4_i*uxy4;
  Real uxy_final = 0.0;
  for (int n=0; n<uxy.size(); n++)
    {
      uxy_final += omega[n]*uxy[n];
    }

  return uxy_final;
} 
Vector<Real> get_charge_scaled(Vector<Real> u_i, Vector<Real> u_e)
{
  Vector<Real> charge_scaled;
  for (int n=0; n<u_i.size(); n++)
    charge_scaled.push_back((r_i*u_i[n] + r_e*u_e[n])/(lambda_d*lambda_d*l_r));
  return charge_scaled;
}
void WENOcharacteristic(const Array4<Real>& arr, const Array4<Real>& slopes,
			int i, int j, int k, int iOffset, int jOffset, int kOffset,
			int start, int len, int d){
  // local characteristic reconstruction		      
  Vector<Real> u_i = get_data_zone(arr,i,j,k,start,start+len);
  // define conserved variables
  Real rho_i = u_i[0];
  Real momX_i = u_i[1];
  Real momY_i = u_i[2];
  Real momZ_i = u_i[3];
  Real E_i = u_i[ENER_I];
		  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real v_squared = get_magnitude_squared(v_x_i,v_y_i,v_z_i);
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  // sound speed
  Real c_i = get_speed(u_i);
  // enthalpy
  Real H_i = (E_i+p_i)/rho_i;
  // useful variables
  Real b1 = (Gamma-1.0)/(c_i*c_i);
  Real b2 = 0.5*v_squared*b1;

  Vector<Vector<Real>> left, right;
  if (d==0){
    left.push_back({0.5*(b2+v_x_i/c_i),-0.5*(b1*v_x_i+1.0/c_i),-0.5*b1*v_y_i,-0.5*b1*v_z_i,0.5*b1});
    left.push_back({1.0-b2,b1*v_x_i,b1*v_y_i,b1*v_z_i,-b1});
    left.push_back({0.5*(b2-v_x_i/c_i),-0.5*(b1*v_x_i-1.0/c_i),-0.5*b1*v_y_i,-0.5*b1*v_z_i,0.5*b1});
    left.push_back({-v_y_i,0.0,1.0,0.0,0.0});
    left.push_back({-v_z_i,0.0,0.0,1.0,0.0});
  } else if (d==1){
    left.push_back({0.5*(b2+v_y_i/c_i),-0.5*b1*v_x_i,-0.5*(b1*v_y_i+1.0/c_i),-0.5*b1*v_z_i,0.5*b1});
    left.push_back({1.0-b2,b1*v_x_i,b1*v_y_i,b1*v_z_i,-b1});
    left.push_back({0.5*(b2-v_y_i/c_i),-0.5*b1*v_x_i,-0.5*(b1*v_y_i-1.0/c_i),-0.5*b1*v_z_i,0.5*b1});
    left.push_back({v_x_i,-1.0,0.0,0.0,0.0});
    left.push_back({v_z_i,0.0,0.0,-1.0,0.0});			
  }
  
  // Conserved variables covering the stencils
  Vector<Real> u_iMinus2 = get_data_zone(arr,i-2*iOffset,j-2*jOffset,k-2*kOffset,start,start+len);
  Vector<Real> u_iMinus1 = get_data_zone(arr,i-iOffset,j-jOffset,k-kOffset,start,start+len);
  Vector<Real> u_iPlus1 = get_data_zone(arr,i+iOffset,j+jOffset,k+kOffset,start,start+len);
  Vector<Real> u_iPlus2 = get_data_zone(arr,i+2*iOffset,j+2*jOffset,k+2*kOffset,start,start+len);
  // Characteristic variables covering the stencils
  Vector<Real> w_iMinus2,w_iMinus1,w_i,w_iPlus1,w_iPlus2;  
  for (int n=0; n<left.size(); n++)
    {
      w_iMinus2.push_back(dotProduct(left[n],u_iMinus2));
      w_iMinus1.push_back(dotProduct(left[n],u_iMinus1));
      w_i.push_back(dotProduct(left[n],u_i));
      w_iPlus1.push_back(dotProduct(left[n],u_iPlus1));
      w_iPlus2.push_back(dotProduct(left[n],u_iPlus2));
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
  Vector<Real> w5 = {w_iMinus2[4],w_iMinus1[4],w_i[4],
		     w_iPlus1[4],w_iPlus2[4]};
		  
  std::array<Real, 2> slopesw1 = WENO3_slope(w1);		      
  std::array<Real, 2> slopesw2 = WENO3_slope(w2);
  std::array<Real, 2> slopesw3 = WENO3_slope(w3);
  std::array<Real, 2> slopesw4 = WENO3_slope(w4);
  std::array<Real, 2> slopesw5 = WENO3_slope(w5);

  if (d==0){
    right.push_back({1.0,1.0,1.0,0.0,0.0});
    right.push_back({v_x_i-c_i,v_x_i,v_x_i+c_i,0.0,0.0});
    right.push_back({v_y_i,v_y_i,v_y_i,1.0,0.0});
    right.push_back({v_z_i,v_z_i,v_z_i,0.0,1.0});
    right.push_back({H_i-v_x_i*c_i,0.5*v_squared,H_i+v_x_i*c_i,v_y_i,v_z_i});
  } else if (d==1){
    right.push_back({1.0,1.0,1.0,0.0,0.0});
    right.push_back({v_x_i,v_x_i,v_x_i,-1.0,0.0});
    right.push_back({v_y_i-c_i,v_y_i,v_y_i+c_i,0.0,0.0});
    right.push_back({v_z_i,v_z_i,v_z_i,0.0,-1.0});
    right.push_back({H_i-v_y_i*c_i,0.5*v_squared,H_i+v_y_i*c_i,-v_x_i,-v_z_i});
  }

  for (int n = 0; n<right.size(); n++)
    {
      slopes(i,j,k,n+2*NUM_STATE_FLUID*d+start) = dotProduct(right[n],{slopesw1[0],slopesw2[0],slopesw3[0],slopesw4[0],slopesw5[0]});
      slopes(i,j,k,n+NUM_STATE_FLUID+2*NUM_STATE_FLUID*d+start) = dotProduct(right[n],{slopesw1[1],slopesw2[1],slopesw3[1],slopesw4[1],slopesw5[1]});
    }
}
