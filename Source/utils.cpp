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
    Real E = u_i[ENER_I];
    Real B_i = u_i[bi];
    Real B_squared = get_magnitude_squared(u_i[BX], u_i[BY], u_i[BZ]);
    // speed of sound
    Real c_s = get_speed(u_i);
    // Alfvén speed
    Real c_a = get_speed_a(rho, B_i);
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
  Real B_x = u_i[BX+d-NUM_STATE_FLUID];
  Real B_y = u_i[BX+(1+d)%3-NUM_STATE_FLUID];
  Real B_z = u_i[BX+(2+d)%3-NUM_STATE_FLUID];
  Real E_x = u_i[EX+d-NUM_STATE_FLUID];
  Real E_y = u_i[EX+(1+d)%3-NUM_STATE_FLUID];
  Real E_z = u_i[EX+(2+d)%3-NUM_STATE_FLUID];

  // flux function
  Vector<Real> function(u_i.size(),0.0);
  function[BX+d-NUM_STATE_FLUID] = 0.0;
  function[BX+(1+d)%3-NUM_STATE_FLUID] = -E_z;
  function[BX+(2+d)%3-NUM_STATE_FLUID] = E_y;
  function[EX+d-NUM_STATE_FLUID] = 0.0;
  function[EX+(1+d)%3-NUM_STATE_FLUID] = c*c*B_z;
  function[EX+(2+d)%3-NUM_STATE_FLUID] = -c*c*B_y;
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
// select the method for slope limiter
// epsilon is the slope limiter
// r is the slope ratio
Real get_epsilon(Real r){
  // MIMBEE
  /* if (r <= 0){
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
  Real rho_L = u_L[0], E_L = u_L[ENER_I], rho_R = u_R[0], E_R = u_R[ENER_I];
  
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
Vector<Real> fluid_flux_HLLC(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset,
			     Real dx, Real dt, int d){
  Vector<Real> ui_iMinus1, ue_iMinus1;
  Vector<Real> ui_i, ue_i;
  Vector<Real> ui_iPlus1, ue_iPlus1;
  Vector<Real> ui_iPlus2, ue_iPlus2;

  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
    {
      ui_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      ui_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      ui_iPlus1.push_back(arr(i,j,k,n));
      ui_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));

      ue_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n+NUM_STATE_FLUID/2));
      ue_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n+NUM_STATE_FLUID/2));
      ue_iPlus1.push_back(arr(i,j,k,n+NUM_STATE_FLUID/2));
      ue_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n+NUM_STATE_FLUID/2));
    }

  // Slope limiting variable index
  Vector<int> limiting_idx(NUM_STATE_FLUID/2, NUM_STATE_FLUID/2-1);

  // Reiamnn problem left state for the ion variables
  Vector<Real> u_i_nPlusHalf_R_i = half_step_evolution_R(ui_iMinus1, ui_i, ui_iPlus1,
							 limiting_idx, dx, dt, d, &fluidFlux); 
  // Reiamnn problem right state for the ion variables
  Vector<Real> u_iPlus1_nPlusHalf_L_i = half_step_evolution_L(ui_i, ui_iPlus1, ui_iPlus2,
							      limiting_idx, dx, dt, d, &fluidFlux);

  // Reiamnn problem left state for the electron variables
  Vector<Real> u_i_nPlusHalf_R_e = half_step_evolution_R(ue_iMinus1, ue_i, ue_iPlus1,
							 limiting_idx, dx, dt, d, &fluidFlux); 
  // Reiamnn problem right state for the electron variables
  Vector<Real> u_iPlus1_nPlusHalf_L_e = half_step_evolution_L(ue_i, ue_iPlus1, ue_iPlus2,
							      limiting_idx, dx, dt, d, &fluidFlux);
  
  // speeds for approximate Riemann problem
  Real S_R_i = get_S_K(u_i_nPlusHalf_R_i, u_iPlus1_nPlusHalf_L_i, d);
  Real S_L_i = -S_R_i;
  Real S_star_i = get_S_star(u_i_nPlusHalf_R_i, u_iPlus1_nPlusHalf_L_i, S_L_i, S_R_i, d);  

  Real S_R_e = get_S_K(u_i_nPlusHalf_R_e, u_iPlus1_nPlusHalf_L_e, d);
  Real S_L_e = -S_R_e;
  Real S_star_e = get_S_star(u_i_nPlusHalf_R_e, u_iPlus1_nPlusHalf_L_e, S_L_e, S_R_e, d);

  // flux depending on the speeds, defined in slides or Toro's book
  Vector<Real> flux(NUM_STATE_FLUID,0.0);
  // Ion HLLC
  if (S_L_i>=0){
    Vector<Real> function = fluidFlux(u_i_nPlusHalf_R_i, d);
    for (int i = 0; i<=ENER_I; i++)
      {       
        flux[i] = function[i];
      }
  } else if (S_L_i<=0 && S_star_i>=0){
    Vector<Real> u_HLLC_L= get_u_HLLC(u_i_nPlusHalf_R_i, S_L_i, S_star_i, d);
    Vector<Real> function = fluidFlux(u_i_nPlusHalf_R_i, d);
    Vector<Real> diff;
    for (int i = 0; i<=ENER_I; i++)
      {
        diff.push_back(S_L_i*(u_HLLC_L[i]-u_i_nPlusHalf_R_i[i]));
	flux[i] = function[i]+diff[i];
      }
  } else if (S_star_i<=0 && S_R_i>=0){
    Vector<Real> u_HLLC_R = get_u_HLLC(u_iPlus1_nPlusHalf_L_i, S_R_i, S_star_i, d);
    Vector<Real> function = fluidFlux(u_iPlus1_nPlusHalf_L_i, d);
    Vector<Real> diff;
    for (int i = 0; i<=ENER_I; i++)
      {
        diff.push_back(S_R_i*(u_HLLC_R[i]-u_iPlus1_nPlusHalf_L_i[i]));
	flux[i] = function[i]+diff[i];
      }
  } else {    
    Vector<Real> function = fluidFlux(u_iPlus1_nPlusHalf_L_i, d);
    for (int i = 0; i<=ENER_I; i++)
      {
        flux[i] = function[i];
      }
  }
  // Electron HLLC
  if (S_L_e>=0){
    Vector<Real> function = fluidFlux(u_i_nPlusHalf_R_e, d);
    for (int i = RHO_E; i<=ENER_E; i++)
      {       
        flux[i] = function[i-5];
      }
  } else if (S_L_e<=0 && S_star_e>=0){
    Vector<Real> u_HLLC_L= get_u_HLLC(u_i_nPlusHalf_R_e, S_L_e, S_star_e, d);
    Vector<Real> function = fluidFlux(u_i_nPlusHalf_R_e, d);
    Vector<Real> diff;
    for (int i = RHO_E; i<=ENER_E; i++)
      {
        diff.push_back(S_L_e*(u_HLLC_L[i-5]-u_i_nPlusHalf_R_e[i-5]));
	flux[i] = function[i-5]+diff[i-5];
      }
  } else if (S_star_e<=0 && S_R_e>=0){
    Vector<Real> u_HLLC_R = get_u_HLLC(u_iPlus1_nPlusHalf_L_e, S_R_e, S_star_e, d);
    Vector<Real> function = fluidFlux(u_iPlus1_nPlusHalf_L_e, d);
    Vector<Real> diff;
    for (int i = RHO_E; i<=ENER_E; i++)
      {
        diff.push_back(S_R_e*(u_HLLC_R[i-5]-u_iPlus1_nPlusHalf_L_e[i-5]));
	flux[i] = function[i-5]+diff[i-5];
      }
  } else {    
    Vector<Real> function = fluidFlux(u_iPlus1_nPlusHalf_L_e, d);
    for (int i = RHO_E; i<=ENER_E; i++)
      {
        flux[i] = function[i-5];
      }
  }

  return flux;
}
Vector<Real> Maxwell_flux_Godunov(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset, 
		       Real dx, Real dt, int d){

  Vector<Real> u_iMinus1;
  Vector<Real> u_i;
  Vector<Real> u_iPlus1;
  Vector<Real> u_iPlus2;

  // Slope limiting variable index
  Vector<int> limiting_idx;

  for (int n = NUM_STATE_FLUID; n<NUM_STATE; n++)
    {
      u_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      u_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_iPlus1.push_back(arr(i,j,k,n));
      u_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));

      // each variable limits itself
      limiting_idx.push_back(n-NUM_STATE_FLUID);
    }
  
  if (CAMReXmp::RKOrder==2)
    {
      /*// Reiamnn problem left state
      Vector<Real> u_i_nPlusHalf_R = half_step_evolution_R(u_iMinus1, u_i, u_iPlus1,
						       limiting_idx, dx, dt, d, &MaxwellFlux); 
      // Reiamnn problem right state
      Vector<Real> u_iPlus1_nPlusHalf_L = half_step_evolution_L(u_i, u_iPlus1, u_iPlus2,
							    limiting_idx, dx, dt, d, &MaxwellFlux);
      */
      // Cell boundary extrapolated values at the left and the right for the ion
      //Vector<Real> u_iL = TVD2_reconstruction_L(u_iMinus1, u_i, u_iPlus1, limiting_idx);
      Vector<Real> u_iR = TVD2_reconstruction_R(u_iMinus1, u_i, u_iPlus1, limiting_idx);
      Vector<Real> u_iPlus1L = TVD2_reconstruction_L(u_i, u_iPlus1, u_iPlus2, limiting_idx);
      //Vector<Real> u_iPlus1R = TVD2_reconstruction_R(u_i, u_iPlus1, u_iPlus2, limiting_idx);

      // Second order in space, not in time
      u_i = u_iR;
      u_iPlus1 = u_iPlus1L;
      
      // Reiamnn problem left state
      //Vector<Real> u_i_nPlusHalf_R = half_update_R(u_iL, u_iR, dx, dt, d, &MaxwellFlux); 
      // Reiamnn problem right state
      //Vector<Real> u_iPlus1_nPlusHalf_L = half_update_L(u_iPlus1L, u_iPlus1R, dx, dt, d, &MaxwellFlux);

      //u_i = u_i_nPlusHalf_R;
      //u_iPlus1 = u_iPlus1_nPlusHalf_L;
    }
  //std::cout << u_i[2] << " " << u_i_nPlusHalf_R[2] << std::endl;
  // godunov flux
  Vector<Real> flux(NUM_STATE_MAXWELL,0.0);
  
  Real B_y_star = 0.5*(u_iPlus1[BX_LOCAL+(1+d)%3]+u_i[BX_LOCAL+(1+d)%3])
    +0.5/c*(u_iPlus1[EX_LOCAL+(2+d)%3]-u_i[EX_LOCAL+(2+d)%3]);
  Real B_z_star = 0.5*(u_iPlus1[BX_LOCAL+(2+d)%3]+u_i[BX_LOCAL+(2+d)%3])
    -0.5/c*(u_iPlus1[EX_LOCAL+(1+d)%3]-u_i[EX_LOCAL+(1+d)%3]);
  Real E_y_star = 0.5*(u_iPlus1[EX_LOCAL+(1+d)%3]+u_i[EX_LOCAL+(1+d)%3])
    -0.5*c*(u_iPlus1[BX_LOCAL+(2+d)%3]-u_i[BX_LOCAL+(2+d)%3]);
  Real E_z_star = 0.5*(u_iPlus1[EX_LOCAL+(2+d)%3]+u_i[EX_LOCAL+(2+d)%3])
    +0.5*c*(u_iPlus1[BX_LOCAL+(1+d)%3]-u_i[BX_LOCAL+(1+d)%3]);  
  
  // EM HLLC states
  flux[BX_LOCAL+d] = 0.0;
  flux[BX_LOCAL+(1+d)%3] = -E_z_star;
  flux[BX_LOCAL+(2+d)%3] = E_y_star;
  flux[EX_LOCAL+d] = 0.0;
  flux[EX_LOCAL+(1+d)%3] = c*c*B_z_star;
  flux[EX_LOCAL+(2+d)%3] =  -c*c*B_y_star;

  return flux;
}
Vector<Real> Maxwell_flux_HLLE(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset, 
		       Real dx, Real dt, int d){

  Vector<Real> u_iMinus1;
  Vector<Real> u_i;
  Vector<Real> u_iPlus1;
  Vector<Real> u_iPlus2;

  // Slope limiting variable index
  Vector<int> limiting_idx;

  for (int n = NUM_STATE_FLUID; n<NUM_STATE; n++)
    {
      u_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      u_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_iPlus1.push_back(arr(i,j,k,n));
      u_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));

      // each variable limits itself
      limiting_idx.push_back(n-NUM_STATE_FLUID);
    }
  
  // Reiamnn problem left state
  Vector<Real> u_i_nPlusHalf_R = half_step_evolution_R(u_iMinus1, u_i, u_iPlus1,
						       limiting_idx, dx, dt, d, &MaxwellFlux); 
  // Reiamnn problem right state
  Vector<Real> u_iPlus1_nPlusHalf_L = half_step_evolution_L(u_i, u_iPlus1, u_iPlus2,
							    limiting_idx, dx, dt, d, &MaxwellFlux);
  if (CAMReXmp::RKOrder==2)
    {
      u_i = u_i_nPlusHalf_R; 
      u_iPlus1 = u_iPlus1_nPlusHalf_L;
    }

  Real S_L = -c;
  Real S_R = c;
  
  // HLLE flux
  Vector<Real> flux(NUM_STATE_MAXWELL,0.0);
  if (S_L>=0){
    Vector<Real> function = MaxwellFlux(u_i, d);
    for (int i = 0; i<NUM_STATE_MAXWELL; i++)
      {       
        flux[i] = function[i];
      }
  } else if (S_L<=0 && S_R>=0){
    //Vector<Real> u_HLLC_L= get_u_HLLC(u_i_nPlusHalf_R_i, S_L_i, S_star_i, d);
    Vector<Real> functionL = MaxwellFlux(u_i, d);
    Vector<Real> functionR = MaxwellFlux(u_iPlus1, d);
    //Vector<Real> diff;
    for (int i = 0; i<NUM_STATE_MAXWELL; i++)
      {
        //diff.push_back(S_L_i*(u_HLLC_L[i]-u_i_nPlusHalf_R_i[i]));
	flux[i] = 1.0/(S_R-S_L)*(S_R*functionL[i] - S_L*functionR[i] + S_R*S_L*(u_iPlus1[i]-u_i[i]));
      }
  } else {    
    Vector<Real> function = MaxwellFlux(u_iPlus1, d);
    for (int i = 0; i<NUM_STATE_MAXWELL; i++)
      {
        flux[i] = function[i];
      }
  }

  return flux;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

Real computeDerivative(Real u_iMinus1, Real u_iPlus1, Real dx){
  return (u_iPlus1-u_iMinus1)/(2.0*dx);
}
Real computeSecondDerivative(Real u_iMinus1, Real u_i, Real u_iPlus1, Real dx){
  return (u_iMinus1-2*u_i+u_iPlus1)/(dx*dx);
}
Real computeCurrentFluxWithPressure(const Vector<Real>& u, int d){

 return r_i*u[MOMX_I+d]*u[MOMX_I+d]/u[RHO_I] +
    r_e*u[MOMX_E+d]*u[MOMX_E+d]/u[RHO_E] +
    r_i*get_pressure({u[RHO_I],u[MOMX_I],u[MOMY_I],u[MOMZ_I],u[ENER_I]}) +
    r_e*get_pressure({u[RHO_E],u[MOMX_E],u[MOMY_E],u[MOMZ_E],u[ENER_E]});
}
Real computeCurrentFluxWithoutPressure(const Vector<Real>& u, int d, int comp){

 return r_i*u[MOMX_I+d]*u[MOMX_I+(comp+d)%3]/u[RHO_I] +
    r_e*u[MOMX_E+d]*u[MOMX_E+(comp+d)%3]/u[RHO_E];
}
Vector<Real> currentUpdate(const Array4<const Real>& arr, int i, int j, int k, const Real* dx, Real dt){
  Vector<Real> u_ij;
  Vector<Vector<Real> > u_iPlus1j, u_iMinus1j;
  u_iPlus1j.resize(SpaceDim, Vector<Real>(NUM_STATE));
  u_iMinus1j.resize(SpaceDim, Vector<Real>(NUM_STATE));
  for (int n = 0; n<NUM_STATE; n++)
    {
      u_ij.push_back(arr(i,j,k,n));      
    }
  for (int d = 0; d < amrex::SpaceDim ; d++)
    {
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
      
      for (int n = 0; n<NUM_STATE; n++)
	{
	  u_iPlus1j[d][n] = arr(i+iOffset,j+jOffset,k+kOffset,n);
	  u_iMinus1j[d][n] = arr(i-iOffset,j-jOffset,k-kOffset,n);
	}
    }
  // define conserved variables                                                                           
  Real rho_i = u_ij[0];
  Real momX_i = u_ij[1];
  Real momY_i = u_ij[2];
  Real momZ_i = u_ij[MOMZ_I];
  Real E_i = u_ij[ENER_I];
  Real rho_e = u_ij[RHO_E];
  Real momX_e = u_ij[MOMX_E];
  Real momY_e = u_ij[MOMY_E];
  Real momZ_e = u_ij[MOMZ_E];
  Real E_e = u_ij[ENER_E];
  Real B_x = u_ij[BX];
  Real B_y = u_ij[BY];
  Real B_z = u_ij[BZ];
  Real E_x = u_ij[EX];
  Real E_y = u_ij[EY];
  Real E_z = u_ij[EZ];
  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});

  Real Jx = r_i*rho_i*v_x_i + r_e*rho_e*v_x_e;
  Real Jy = r_i*rho_i*v_y_i + r_e*rho_e*v_y_e;
  Real Jz = r_i*rho_i*v_z_i + r_e*rho_e*v_z_e;
  Real charge_scaled = r_i*r_i*rho_i + r_e*r_e*rho_e;
  /*Vector<Real> current = {r_i*rho_i*v_x_i + r_e*rho_e*v_x_e,
			  r_i*rho_i*v_y_i + r_e*rho_e*v_y_e,
			  r_i*rho_i*v_z_i + r_e*rho_e*v_z_e};*/
  Vector<Real> current_scaled = {r_i*r_i*rho_i*v_x_i + r_e*r_e*rho_e*v_x_e,
				 r_i*r_i*rho_i*v_y_i + r_e*r_e*rho_e*v_y_e,
				 r_i*r_i*rho_i*v_z_i + r_e*r_e*rho_e*v_z_e};
  //Vector<Real> currentUpdated = cross_product(current_scaled,{B_x, B_y, B_z});
  //currentUpdated[0] += 
  Vector<Real> currentUpdated = {Jx + dt*(charge_scaled*E_x + current_scaled[1]*B_z - current_scaled[2]*B_y)/l_r,
				 Jy + dt*(charge_scaled*E_y + current_scaled[2]*B_x - current_scaled[0]*B_z)/l_r,
				 Jz + dt*(charge_scaled*E_z + current_scaled[0]*B_y - current_scaled[1]*B_x)/l_r};
  
  for (int d = 0; d < amrex::SpaceDim ; d++)
    {      
      Real PxxMinus1 = computeCurrentFluxWithPressure(u_iMinus1j[d],d);
      Real PxxPlus1 = computeCurrentFluxWithPressure(u_iPlus1j[d],d); 

      Real PyxMinus1 = computeCurrentFluxWithoutPressure(u_iMinus1j[d],d,1);
      Real PyxPlus1 = computeCurrentFluxWithoutPressure(u_iPlus1j[d],d,1); 

      Real PzxMinus1 = computeCurrentFluxWithoutPressure(u_iMinus1j[d],d,2);
      Real PzxPlus1 = computeCurrentFluxWithoutPressure(u_iPlus1j[d],d,2);

      Real dxPxx = computeDerivative(PxxMinus1, PxxPlus1, dx[d]);
      Real dxPyx = computeDerivative(PyxMinus1, PyxPlus1, dx[d]);
      Real dxPzx = computeDerivative(PzxMinus1, PzxPlus1, dx[d]);

      currentUpdated[d] -= dt*dxPxx;
      currentUpdated[(1+d)%3] -= dt*dxPyx;
      currentUpdated[(2+d)%3] -= dt*dxPzx;

      // cylindrical source terms
      if (AMReX::top()->getDefaultGeometry()->Coord()==1 && d==1)
	{
	  // y <- x in cylindrical
	  const Real y = AMReX::top()->getDefaultGeometry()->ProbLo()[1]+(double(j)+0.5)*dx[1];
	  Real Pyx = computeCurrentFluxWithoutPressure(u_ij,d,0);
	  Real Pyy = computeCurrentFluxWithoutPressure(u_ij,d,1);
	  Real Pyz = computeCurrentFluxWithoutPressure(u_ij,d,2);
	  
	  currentUpdated[d] -= dt*Pyx/std::abs(y);
	  currentUpdated[(1+d)%3] -= dt*Pyy/std::abs(y);
	  currentUpdated[(2+d)%3] -= dt*Pyz/std::abs(y); 
	}
    }

  return currentUpdated;
}
/*
Real computeJx_nPlusHalf(const Array4<const Real>& arr, int i, int j, int k,
			Real dx, Real dy, Real dt){
  Vector<Real> u_ij, u_iPlus1j, u_iMinus1j, u_ijPlus1, u_ijMinus1;
  for (int n = 0; n<16; n++)
    {
      u_ij.push_back(arr(i,j,k,n));      
      u_iPlus1j.push_back(arr(i+1,j,k,n));
      u_iMinus1j.push_back(arr(i-1,j,k,n));
      // 2D
      u_ijPlus1.push_back(arr(i,j+1,k,n));
      u_ijMinus1.push_back(arr(i,j-1,k,n));
    }
  
  // define conserved variables                                                                                               
  Real rho_i = u_ij[0];
  Real momX_i = u_ij[1];
  Real momY_i = u_ij[2];
  Real momZ_i = u_ij[MOMZ_I];
  Real E_i = u_ij[ENER_I];
  Real rho_e = u_ij[RHO_E];
  Real momX_e = u_ij[MOMX_E];
  Real momY_e = u_ij[MOMY_E];
  Real momZ_e = u_ij[MOMZ_E];
  Real E_e = u_ij[ENER_E];
  Real B_x = u_ij[BX];
  Real B_y = u_ij[BY];
  Real B_z = u_ij[BZ];
  Real E_x = u_ij[EX];
  Real E_y = u_ij[EY];
  Real E_z = u_ij[EZ];
  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});

  Real PxxMinus1 = r_i*u_iMinus1j[MOMX_I]*u_iMinus1j[MOMX_I]/u_iMinus1j[RHO_I] +
    r_e*u_iMinus1j[MOMX_E]*u_iMinus1j[MOMX_E]/u_iMinus1j[RHO_E] +
    r_i*get_pressure({u_iMinus1j[RHO_I],u_iMinus1j[MOMX_I],u_iMinus1j[MOMY_I],u_iMinus1j[MOMZ_I],u_iMinus1j[ENER_I]}) +
    r_e*get_pressure({u_iMinus1j[RHO_E],u_iMinus1j[MOMX_E],u_iMinus1j[MOMY_E],u_iMinus1j[MOMZ_E],u_iMinus1j[ENER_E]});
  Real PxxPlus1 = r_i*u_iPlus1j[MOMX_I]*u_iPlus1j[MOMX_I]/u_iPlus1j[RHO_I] +
    r_e*u_iPlus1j[MOMX_E]*u_iPlus1j[MOMX_E]/u_iPlus1j[RHO_E] +
    r_i*get_pressure({u_iPlus1j[RHO_I],u_iPlus1j[MOMX_I],u_iPlus1j[MOMY_I],u_iPlus1j[MOMZ_I],u_iPlus1j[ENER_I]}) +
    r_e*get_pressure({u_iPlus1j[RHO_E],u_iPlus1j[MOMX_E],u_iPlus1j[MOMY_E],u_iPlus1j[MOMZ_E],u_iPlus1j[ENER_E]});

  /////////////////////////////////////////////////////////////////////////  
  // 2D
  Real PxyMinus1 = r_i*u_ijMinus1[MOMX_I]*u_ijMinus1[MOMY_I]/u_ijMinus1[RHO_I] +
    r_e*u_ijMinus1[MOMX_E]*u_ijMinus1[MOMY_E]/u_ijMinus1[RHO_E];
  Real PxyPlus1 = r_i*u_ijPlus1[MOMX_I]*u_ijPlus1[MOMY_I]/u_ijPlus1[RHO_I] +
    r_e*u_ijPlus1[MOMX_E]*u_ijPlus1[MOMY_E]/u_ijPlus1[RHO_E];
  /////////////////////////////////////////////////////////////////////////
  
  Real Jx = r_i*rho_i*v_x_i + r_e*rho_e*v_x_e;
  Real dxPxx = computeDerivative(PxxMinus1, PxxPlus1, dx);
  /////////////////////////////////////////////////////////////////////////  
  // 2D
  Real dyPxy = computeDerivative(PxyMinus1, PxyPlus1, dy);
  /////////////////////////////////////////////////////////////////////////
  Real charge_scaled = r_i*r_i*rho_i + r_e*r_e*rho_e;
  Vector<Real> current_scaled = {r_i*r_i*rho_i*v_x_i + r_e*r_e*rho_e*v_x_e,
				 r_i*r_i*rho_i*v_y_i + r_e*r_e*rho_e*v_y_e,
				 r_i*r_i*rho_i*v_z_i + r_e*r_e*rho_e*v_z_e};

  // 2D
  return Jx - 0.5*dt*(dxPxx+dyPxy) + 0.5*dt*(charge_scaled*E_x + current_scaled[1]*B_z - current_scaled[2]*B_y)/l_r;
  // 1D
  //return Jx - 0.5*dt*dxPxx + 0.5*dt*(charge_scaled*E_x + current_scaled[1]*B_z - current_scaled[2]*B_y)/l_r;   
}
Real computeJy_nPlusHalf(const Array4<const Real>& arr, int i, int j, int k,
			Real dx, Real dy, Real dt){
  Vector<Real> u_ij, u_iPlus1j, u_iMinus1j, u_ijPlus1, u_ijMinus1;
  for (int n = 0; n<16; n++)
    {
      u_ij.push_back(arr(i,j,k,n));      
      u_iPlus1j.push_back(arr(i+1,j,k,n));
      u_iMinus1j.push_back(arr(i-1,j,k,n));
      // 2D
      u_ijPlus1.push_back(arr(i,j+1,k,n));
      u_ijMinus1.push_back(arr(i,j-1,k,n));
    }
  
  // define conserved variables                                                                                               
  Real rho_i = u_ij[0];
  Real momX_i = u_ij[1];
  Real momY_i = u_ij[2];
  Real momZ_i = u_ij[MOMZ_I];
  Real E_i = u_ij[ENER_I];
  Real rho_e = u_ij[RHO_E];
  Real momX_e = u_ij[MOMX_E];
  Real momY_e = u_ij[MOMY_E];
  Real momZ_e = u_ij[MOMZ_E];
  Real E_e = u_ij[ENER_E];
  Real B_x = u_ij[BX];
  Real B_y = u_ij[BY];
  Real B_z = u_ij[BZ];
  Real E_x = u_ij[EX];
  Real E_y = u_ij[EY];
  Real E_z = u_ij[EZ];
  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});

  Real PyxMinus1 = r_i*u_iMinus1j[MOMX_I]*u_iMinus1j[MOMY_I]/u_iMinus1j[RHO_I] +
    r_e*u_iMinus1j[MOMX_E]*u_iMinus1j[MOMY_E]/u_iMinus1j[RHO_E];
  Real PyxPlus1 = r_i*u_iPlus1j[MOMX_I]*u_iPlus1j[MOMY_I]/u_iPlus1j[RHO_I] +
    r_e*u_iPlus1j[MOMX_E]*u_iPlus1j[MOMY_E]/u_iPlus1j[RHO_E];
  /////////////////////////////////////////////////////////////////////////
  // 2D
  Real PyyMinus1 = r_i*u_ijMinus1[MOMY_I]*u_ijMinus1[MOMY_I]/u_ijMinus1[RHO_I] +
    r_e*u_ijMinus1[MOMY_E]*u_ijMinus1[MOMY_E]/u_ijMinus1[RHO_E] +
    r_i*get_pressure({u_ijMinus1[RHO_I],u_ijMinus1[MOMX_I],u_ijMinus1[MOMY_I],u_ijMinus1[MOMZ_I],u_ijMinus1[ENER_I]}) +
    r_e*get_pressure({u_ijMinus1[RHO_E],u_ijMinus1[MOMX_E],u_ijMinus1[MOMY_E],u_ijMinus1[MOMZ_E],u_ijMinus1[ENER_E]});
  Real PyyPlus1 = r_i*u_ijPlus1[MOMY_I]*u_ijPlus1[MOMY_I]/u_ijPlus1[RHO_I] +
    r_e*u_ijPlus1[MOMY_E]*u_ijPlus1[MOMY_E]/u_ijPlus1[RHO_E] +
    r_i*get_pressure({u_ijPlus1[RHO_I],u_ijPlus1[MOMX_I],u_ijPlus1[MOMY_I],u_ijPlus1[MOMZ_I],u_ijPlus1[ENER_I]}) +
    r_e*get_pressure({u_ijPlus1[RHO_E],u_ijPlus1[MOMX_E],u_ijPlus1[MOMY_E],u_ijPlus1[MOMZ_E],u_ijPlus1[ENER_E]});
  /////////////////////////////////////////////////////////////////////////

  Real Jy = r_i*rho_i*v_y_i + r_e*rho_e*v_y_e;
  Real dxPyx = computeDerivative(PyxMinus1, PyxPlus1, dx);
  /////////////////////////////////////////////////////////////////////////
  // 2D
  Real dyPyy = computeDerivative(PyyMinus1, PyyPlus1, dy);
  /////////////////////////////////////////////////////////////////////////  
  Real charge_scaled = r_i*r_i*rho_i + r_e*r_e*rho_e;
  Vector<Real> current_scaled = {r_i*r_i*rho_i*v_x_i + r_e*r_e*rho_e*v_x_e,
				 r_i*r_i*rho_i*v_y_i + r_e*r_e*rho_e*v_y_e,
				 r_i*r_i*rho_i*v_z_i + r_e*r_e*rho_e*v_z_e};
  // 2D
  return Jy - 0.5*dt*(dxPyx+dyPyy) + 0.5*dt*(charge_scaled*E_y + current_scaled[2]*B_z - current_scaled[0]*B_z)/l_r;
  // 1D
  //return Jy - 0.5*dt*dxPyx + 0.5*dt*(charge_scaled*E_y + current_scaled[2]*B_x - current_scaled[0]*B_z)/l_r;  
}
Real computeJz_nPlusHalf(const Array4<const Real>& arr, int i, int j, int k,
			Real dx, Real dy, Real dt){
  Vector<Real> u_ij, u_iPlus1j, u_iMinus1j, u_ijPlus1, u_ijMinus1;
  for (int n = 0; n<16; n++)
    {
      u_ij.push_back(arr(i,j,k,n));      
      u_iPlus1j.push_back(arr(i+1,j,k,n));
      u_iMinus1j.push_back(arr(i-1,j,k,n));
      // 2D
      u_ijPlus1.push_back(arr(i,j+1,k,n));
      u_ijMinus1.push_back(arr(i,j-1,k,n));
    }
  
  // define conserved variables                                                                                               
  Real rho_i = u_ij[0];
  Real momX_i = u_ij[1];
  Real momY_i = u_ij[2];
  Real momZ_i = u_ij[MOMZ_I];
  Real E_i = u_ij[ENER_I];
  Real rho_e = u_ij[RHO_E];
  Real momX_e = u_ij[MOMX_E];
  Real momY_e = u_ij[MOMY_E];
  Real momZ_e = u_ij[MOMZ_E];
  Real E_e = u_ij[ENER_E];
  Real B_x = u_ij[BX];
  Real B_y = u_ij[BY];
  Real B_z = u_ij[BZ];
  Real E_x = u_ij[EX];
  Real E_y = u_ij[EY];
  Real E_z = u_ij[EZ];
  
  // define primitive variables
  Real v_x_i = momX_i/rho_i;
  Real v_y_i = momY_i/rho_i;
  Real v_z_i = momZ_i/rho_i;
  Real p_i = get_pressure({rho_i,momX_i,momY_i,momZ_i,E_i});
  Real v_x_e = momX_e/rho_e;
  Real v_y_e = momY_e/rho_e;
  Real v_z_e = momZ_e/rho_e;
  Real p_e = get_pressure({rho_e,momX_e,momY_e,momZ_e,E_e});

  Real PzxMinus1 = r_i*u_iMinus1j[MOMX_I]*u_iMinus1j[MOMZ_I]/u_iMinus1j[RHO_I] +
    r_e*u_iMinus1j[MOMX_E]*u_iMinus1j[MOMZ_E]/u_iMinus1j[RHO_E];
  Real PzxPlus1 = r_i*u_iPlus1j[MOMX_I]*u_iPlus1j[MOMZ_I]/u_iPlus1j[RHO_I] +
    r_e*u_iPlus1j[MOMX_E]*u_iPlus1j[MOMZ_E]/u_iPlus1j[RHO_E];
  /////////////////////////////////////////////////////////////////////////
  // 2D
  Real PzyMinus1 = r_i*u_ijMinus1[MOMY_I]*u_ijMinus1[MOMZ_I]/u_ijMinus1[RHO_I] +
    r_e*u_ijMinus1[MOMY_E]*u_ijMinus1[MOMZ_E]/u_ijMinus1[RHO_E];
  Real PzyPlus1 = r_i*u_ijPlus1[MOMY_I]*u_ijPlus1[MOMZ_I]/u_ijPlus1[RHO_I] +
    r_e*u_ijPlus1[MOMY_E]*u_ijPlus1[MOMZ_E]/u_ijPlus1[RHO_E];
  /////////////////////////////////////////////////////////////////////////
  
  Real Jz = r_i*rho_i*v_z_i + r_e*rho_e*v_z_e;
  Real dxPzx = computeDerivative(PzxMinus1, PzxPlus1, dx);
  /////////////////////////////////////////////////////////////////////////
  // 2D
  Real dyPzy = computeDerivative(PzyMinus1, PzyPlus1, dy);
  /////////////////////////////////////////////////////////////////////////  
  Real charge_scaled = r_i*r_i*rho_i + r_e*r_e*rho_e;
  Vector<Real> current_scaled = {r_i*r_i*rho_i*v_x_i + r_e*r_e*rho_e*v_x_e,
				 r_i*r_i*rho_i*v_y_i + r_e*r_e*rho_e*v_y_e,
				 r_i*r_i*rho_i*v_z_i + r_e*r_e*rho_e*v_z_e};

  // 2D
  return Jz - 0.5*dt*(dxPzx+dyPzy) + 0.5*dt*(charge_scaled*E_z + current_scaled[0]*B_y - current_scaled[1]*B_x)/l_r;
  // 1D
  //return Jz - 0.5*dt*dxPzx + 0.5*dt*(charge_scaled*E_z + current_scaled[0]*B_y - current_scaled[1]*B_x)/l_r;
}
*/

Vector<Real> half_step_evolution_L(Vector<Real> u_iMinus1, Vector<Real> u_i,
                            Vector<Real> u_iPlus1, Real dx, Real dt, int d){
    
    // slope measure
    Vector<Real> delta_i = get_delta_i(u_iMinus1, u_i, u_iPlus1);
    // slope ratio 
    // index ENER corresponds to energy (limiting on the energy)
    Real r_i = get_r(u_iMinus1[ENER_I], u_i[ENER_I], u_iPlus1[ENER_I]);
    Real r_e = get_r(u_iMinus1[ENER_E], u_i[ENER_E], u_iPlus1[ENER_E]);
    // slope limiter
    Real epsilon_i = get_epsilon(r_i);
    Real epsilon_e = get_epsilon(r_e);
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
    Real r_i = get_r(u_iMinus1[ENER_I], u_i[ENER_I], u_iPlus1[ENER_I]);
    Real r_e = get_r(u_iMinus1[ENER_E], u_i[ENER_E], u_iPlus1[ENER_E]);
    // slope limiter                                                                      
    Real epsilon_i = get_epsilon(r_i);
    Real epsilon_e = get_epsilon(r_e);
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

Vector<Real> flux_HLLC(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset, 
		       Real dx, Real dt, int d){
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

    // define vectors for electron conserved variables
  Vector<Real> u_i_nPlusHalf_R_e = {u_i_nPlusHalf_R[RHO_E],u_i_nPlusHalf_R[MOMX_E],u_i_nPlusHalf_R[MOMY_E],
				    u_i_nPlusHalf_R[MOMZ_E],u_i_nPlusHalf_R[ENER_E]};
  Vector<Real> u_iPlus1_nPlusHalf_L_e = {u_iPlus1_nPlusHalf_L[RHO_E],u_iPlus1_nPlusHalf_L[MOMX_E],
					 u_iPlus1_nPlusHalf_L[MOMY_E],u_iPlus1_nPlusHalf_L[MOMZ_E],
					 u_iPlus1_nPlusHalf_L[ENER_E]};

  //amrex::Print() << u_i_nPlusHalf_R[0] << " " << u_iPlus1_nPlusHalf_L[0] << " " << u_i_nPlusHalf_R_e[0] << " " << u_iPlus1_nPlusHalf_L_e[0] << std::endl;
  
  // speeds for approximate Riemann problem
  Real S_R_i = get_S_K(u_i_nPlusHalf_R, u_iPlus1_nPlusHalf_L, d);
  Real S_L_i = -S_R_i;
  Real S_star_i = get_S_star(u_i_nPlusHalf_R, u_iPlus1_nPlusHalf_L, S_L_i, S_R_i, d);  

  Real S_R_e = get_S_K(u_i_nPlusHalf_R_e, u_iPlus1_nPlusHalf_L_e, d);
  Real S_L_e = -S_R_e;
  Real S_star_e = get_S_star(u_i_nPlusHalf_R_e, u_iPlus1_nPlusHalf_L_e, S_L_e, S_R_e, d);

  // flux depending on the speeds, defined in slides or Toro's book
  Vector<Real> flux(NUM_STATE,0.0);
  // Ion HLLC
  if (S_L_i>=0){
    Vector<Real> function = func(u_i_nPlusHalf_R, d);
    for (int i = 0; i<=ENER_I; i++)
      {       
        flux[i] = function[i];
	//	std::cout << function[i] << " ";
      }
  } else if (S_L_i<=0 && S_star_i>=0){
    Vector<Real> u_HLLC_L= get_u_HLLC(u_i_nPlusHalf_R, S_L_i, S_star_i, d);
    Vector<Real> function = func(u_i_nPlusHalf_R, d);
    Vector<Real> diff;
    for (int i = 0; i<=ENER_I; i++)
      {
        diff.push_back(S_L_i*(u_HLLC_L[i]-u_i_nPlusHalf_R[i]));
	flux[i] = function[i]+diff[i];
	//	std::cout << function[i]+diff[i] << " ";
      }
  } else if (S_star_i<=0 && S_R_i>=0){
    Vector<Real> u_HLLC_R = get_u_HLLC(u_iPlus1_nPlusHalf_L, S_R_i, S_star_i, d);
    Vector<Real> function = func(u_iPlus1_nPlusHalf_L, d);
    Vector<Real> diff;
    for (int i = 0; i<=ENER_I; i++)
      {
        diff.push_back(S_R_i*(u_HLLC_R[i]-u_iPlus1_nPlusHalf_L[i]));
	flux[i] = function[i]+diff[i];
	//	std::cout << function[i]+diff[i] << " ";
      }
  } else {    
    Vector<Real> function = func(u_iPlus1_nPlusHalf_L, d);
    for (int i = 0; i<=ENER_I; i++)
      {
        flux[i] = function[i];
	//	std::cout << function[i] << " ";
      }
  }
  // Electron HLLC
  if (S_L_e>=0){
    Vector<Real> function = func(u_i_nPlusHalf_R, d);
    for (int i = RHO_E; i<=ENER_E; i++)
      {       
        flux[i] = function[i];
	//	std::cout << function[i] << " ";
      }
  } else if (S_L_e<=0 && S_star_e>=0){
    Vector<Real> u_HLLC_L= get_u_HLLC(u_i_nPlusHalf_R_e, S_L_e, S_star_e, d);
    Vector<Real> function = func(u_i_nPlusHalf_R, d);
    Vector<Real> diff;
    for (int i = RHO_E; i<=ENER_E; i++)
      {
        diff.push_back(S_L_e*(u_HLLC_L[i-5]-u_i_nPlusHalf_R[i]));
	flux[i] = function[i]+diff[i-5];
	//	std::cout << function[i]+diff[i-5] << " ";
      }
  } else if (S_star_e<=0 && S_R_e>=0){
    Vector<Real> u_HLLC_R = get_u_HLLC(u_iPlus1_nPlusHalf_L_e, S_R_e, S_star_e, d);
    Vector<Real> function = func(u_iPlus1_nPlusHalf_L, d);
    Vector<Real> diff;
    for (int i = RHO_E; i<=ENER_E; i++)
      {
        diff.push_back(S_R_e*(u_HLLC_R[i-5]-u_iPlus1_nPlusHalf_L[i]));
	flux[i] = function[i]+diff[i-5];
	//	std::cout << function[i]+diff[i-5] << " ";
      }
  } else {    
    Vector<Real> function = func(u_iPlus1_nPlusHalf_L, d);
    for (int i = RHO_E; i<=ENER_E; i++)
      {
        flux[i] = function[i];
	//	std::cout << function[i] << " ";
      }
  }

  Real B_y_star = 0.5*(u_iPlus1[BX+(1+d)%3]+u_i[BX+(1+d)%3])+0.5/c*(u_iPlus1[EX+(2+d)%3]-u_i[EX+(2+d)%3]);
  Real B_z_star = 0.5*(u_iPlus1[BX+(2+d)%3]+u_i[BX+(2+d)%3])-0.5/c*(u_iPlus1[EX+(1+d)%3]-u_i[EX+(1+d)%3]);
  Real E_y_star = 0.5*(u_iPlus1[EX+(1+d)%3]+u_i[EX+(1+d)%3])-0.5*c*(u_iPlus1[BX+(2+d)%3]-u_i[BX+(2+d)%3]);
  Real E_z_star = 0.5*(u_iPlus1[EX+(2+d)%3]+u_i[EX+(2+d)%3])+0.5*c*(u_iPlus1[BX+(1+d)%3]-u_i[BX+(1+d)%3]);
  // EM HLLC states
  flux[BX+d] = 0.0;
  flux[BX+(1+d)%3] = -E_z_star;
  flux[BX+(2+d)%3] = E_y_star;
  flux[EX+d] = 0.0;
  flux[EX+(1+d)%3] = c*c*B_z_star;
  flux[EX+(2+d)%3] =  -c*c*B_y_star;
  
  return flux;
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

Vector<Real> MUSCL_Hancock_HLLC_flux(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset,
				     Real dx, Real dt, int d){
  Vector<Real> ui_iMinus1, ue_iMinus1;
  Vector<Real> ui_i, ue_i;
  Vector<Real> ui_iPlus1, ue_iPlus1;
  Vector<Real> ui_iPlus2, ue_iPlus2;

  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
    {
      ui_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      ui_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      ui_iPlus1.push_back(arr(i,j,k,n));
      ui_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));
      
      ue_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n+NUM_STATE_FLUID/2));
      ue_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n+NUM_STATE_FLUID/2));
      ue_iPlus1.push_back(arr(i,j,k,n+NUM_STATE_FLUID/2));
      ue_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n+NUM_STATE_FLUID/2));
    }

  // Slope limiting variable index
  Vector<int> limiting_idx(NUM_STATE_FLUID/2, NUM_STATE_FLUID/2-1);
  
  // Cell boundary extrapolated values at the left and the right for the ion
  Vector<Real> u_iL_i = TVD2_reconstruction_L(ui_iMinus1, ui_i, ui_iPlus1, limiting_idx);
  Vector<Real> u_iR_i = TVD2_reconstruction_R(ui_iMinus1, ui_i, ui_iPlus1, limiting_idx);
  Vector<Real> u_iPlus1L_i = TVD2_reconstruction_L(ui_i, ui_iPlus1, ui_iPlus2, limiting_idx);
  Vector<Real> u_iPlus1R_i = TVD2_reconstruction_R(ui_i, ui_iPlus1, ui_iPlus2, limiting_idx);

  // Cell boundary extrapolated values at the left and the right for the electron
  Vector<Real> u_iL_e = TVD2_reconstruction_L(ue_iMinus1, ue_i, ue_iPlus1, limiting_idx);
  Vector<Real> u_iR_e = TVD2_reconstruction_R(ue_iMinus1, ue_i, ue_iPlus1, limiting_idx);
  Vector<Real> u_iPlus1L_e = TVD2_reconstruction_L(ue_i, ue_iPlus1, ue_iPlus2, limiting_idx);
  Vector<Real> u_iPlus1R_e = TVD2_reconstruction_R(ue_i, ue_iPlus1, ue_iPlus2, limiting_idx);
  
  // Reiamnn problem left state for the ion variables 
  Vector<Real> u_i_nPlusHalf_R_i = half_update_R(u_iL_i, u_iR_i, dx, dt, d, &fluidFlux); 
  // Reiamnn problem right state for the ion variables
  Vector<Real> u_iPlus1_nPlusHalf_L_i = half_update_L(u_iPlus1L_i, u_iPlus1R_i, dx, dt, d, &fluidFlux);

  // Reiamnn problem left state for the electron variables
  Vector<Real> u_i_nPlusHalf_R_e = half_update_R(u_iL_e, u_iR_e, dx, dt, d, &fluidFlux); 
  // Reiamnn problem right state for the electron variables
  Vector<Real> u_iPlus1_nPlusHalf_L_e = half_update_L(u_iPlus1L_e, u_iPlus1R_e, dx, dt, d, &fluidFlux);
  /*
  if ((i==256) && (j>252)){
    std::cout << "Reconstruction: " << i << " " << j << " ";
    for(int n=0;n<5;n++)
      std::cout << u_iL_i[n] << " " << u_iR_i[n] << " " << u_iPlus1L_i[n] << " " << u_iPlus1R_i[n] << " " << u_i_nPlusHalf_R_i[n] << " " << u_iPlus1_nPlusHalf_L_i[n] << " " << ui_iMinus1[n] << " " << ui_i[n] << " " <<  ui_iPlus1[n] << " " <<  ui_iPlus2[n] << std::endl;
    Vector<Real> func_iR = fluidFlux(u_iL_i, d);
    Vector<Real> func_iL = fluidFlux(u_iR_i, d);
    Vector<Real> diff;
    Vector<Real> half_step_evolution;
    for (int n = 0; n<u_iR_i.size(); n++)
      {
	diff.push_back(func_iR[n]-func_iL[n]);
	half_step_evolution.push_back(u_iR_i[n] - 0.5*(dt/dx)*diff[n]);
	std::cout << "Half evolution " << u_iR_i[n] << " " << func_iR[n] << " " << func_iL[n] << " " << diff[n] << std::endl;
      }
    
  }
  */
  // 1st order in space
  /*
  Vector<Real> u_i_nPlusHalf_R_i = ui_i;
  Vector<Real> u_iPlus1_nPlusHalf_L_i = ui_iPlus1;
  Vector<Real> u_i_nPlusHalf_R_e = ue_i;
  Vector<Real> u_iPlus1_nPlusHalf_L_e = ue_iPlus1;
  */
  //amrex::Print() << u_i_nPlusHalf_R_i[0] << " " << u_iPlus1_nPlusHalf_L_i[0] << " " << u_i_nPlusHalf_R_e[0] << " " << u_iPlus1_nPlusHalf_L_e[0] << std::endl;
  
  // speeds for approximate Riemann problem
  Real S_R_i = get_S_K(u_i_nPlusHalf_R_i, u_iPlus1_nPlusHalf_L_i, d);
  Real S_L_i = -S_R_i;
  Real S_star_i = get_S_star(u_i_nPlusHalf_R_i, u_iPlus1_nPlusHalf_L_i, S_L_i, S_R_i, d);  

  Real S_R_e = get_S_K(u_i_nPlusHalf_R_e, u_iPlus1_nPlusHalf_L_e, d);
  Real S_L_e = -S_R_e;
  Real S_star_e = get_S_star(u_i_nPlusHalf_R_e, u_iPlus1_nPlusHalf_L_e, S_L_e, S_R_e, d);
  
  // flux depending on the speeds, defined in slides or Toro's book
  Vector<Real> flux(NUM_STATE_FLUID,0.0);
  // Ion HLLC
  if (S_L_i>=0){
    Vector<Real> function = fluidFlux(u_i_nPlusHalf_R_i, d);
    for (int n = 0; n<=ENER_I; n++)
      {       
        flux[n] = function[n];
	//if ((i==256) && (j>253))
	//std::cout << "Fun1: " << i << " " << j << " " << flux[n] << " " << function[n] << " " << u_i_nPlusHalf_R_i[n] << std::endl;
      }
  } else if (S_L_i<=0 && S_star_i>=0){
    Vector<Real> u_HLLC_L= get_u_HLLC(u_i_nPlusHalf_R_i, S_L_i, S_star_i, d);
    Vector<Real> function = fluidFlux(u_i_nPlusHalf_R_i, d);
    Vector<Real> diff;
    for (int n = 0; n<=ENER_I; n++)
      {
        diff.push_back(S_L_i*(u_HLLC_L[n]-u_i_nPlusHalf_R_i[n]));
	flux[n] = function[n]+diff[n];
	//if ((i==256) && (j>252)){
	//std::cout << "Fun2: " << i << " " << j << " " << flux[n] << " " << function[n] << " " << u_i_nPlusHalf_R_i[n] << " " << u_HLLC_L[n] << " " << diff[n] << std::endl;
	  //std::cout << S_L_i << " " << S_star_i << " " << u_i_nPlusHalf_R_i[2]/u_i_nPlusHalf_R_i[0] << std::endl;
	  //std::cout << get_pressure(u_i_nPlusHalf_R_i) << " " << get_pressure( u_iPlus1_nPlusHalf_L_i) << " " << S_R_i << " " << u_iPlus1_nPlusHalf_L_i[2]/u_iPlus1_nPlusHalf_L_i[0] << std::endl;
	//}
      }
  } else if (S_star_i<=0 && S_R_i>=0){
    Vector<Real> u_HLLC_R = get_u_HLLC(u_iPlus1_nPlusHalf_L_i, S_R_i, S_star_i, d);
    Vector<Real> function = fluidFlux(u_iPlus1_nPlusHalf_L_i, d);
    Vector<Real> diff;
    for (int n = 0; n<=ENER_I; n++)
      {
        diff.push_back(S_R_i*(u_HLLC_R[n]-u_iPlus1_nPlusHalf_L_i[n]));
	flux[n] = function[n]+diff[n];
	//if ((i==32) && (j>28))
	//std::cout << "Fun3: " << i << " " << j << " " << flux[n] << " " << function[n] << " " << u_iPlus1_nPlusHalf_L_i[n] << " " << u_HLLC_R[n] << std::endl;
      }
  } else {    
    Vector<Real> function = fluidFlux(u_iPlus1_nPlusHalf_L_i, d);
    for (int n = 0; n<=ENER_I; n++)
      {
        flux[n] = function[n];
	//if ((i==256) && (j>253))
	//std::cout << "Fun4: " << i << " " << j << " " << flux[n] << " " << function[n] << " " << u_iPlus1_nPlusHalf_L_i[n] << std::endl;

      }
  }
  // Electron HLLC
  if (S_L_e>=0){
    Vector<Real> function = fluidFlux(u_i_nPlusHalf_R_e, d);
    for (int n = RHO_E; n<=ENER_E; n++)
      {       
        flux[n] = function[n-5];
      }
  } else if (S_L_e<=0 && S_star_e>=0){
    Vector<Real> u_HLLC_L= get_u_HLLC(u_i_nPlusHalf_R_e, S_L_e, S_star_e, d);
    Vector<Real> function = fluidFlux(u_i_nPlusHalf_R_e, d);
    Vector<Real> diff;
    for (int n = RHO_E; n<=ENER_E; n++)
      {
        diff.push_back(S_L_e*(u_HLLC_L[n-5]-u_i_nPlusHalf_R_e[n-5]));
	flux[n] = function[n-5]+diff[n-5];
      }
  } else if (S_star_e<=0 && S_R_e>=0){
    Vector<Real> u_HLLC_R = get_u_HLLC(u_iPlus1_nPlusHalf_L_e, S_R_e, S_star_e, d);
    Vector<Real> function = fluidFlux(u_iPlus1_nPlusHalf_L_e, d);
    Vector<Real> diff;
    for (int n = RHO_E; n<=ENER_E; n++)
      {
        diff.push_back(S_R_e*(u_HLLC_R[n-5]-u_iPlus1_nPlusHalf_L_e[n-5]));
	flux[n] = function[n-5]+diff[n-5];
      }
  } else {    
    Vector<Real> function = fluidFlux(u_iPlus1_nPlusHalf_L_e, d);
    for (int n = RHO_E; n<=ENER_E; n++)
      {
        flux[n] = function[n-5];
      }
  }

  return flux;
}
std::array<Vector<Real>, 2> MUSCL_Hancock_Godunov_Maxwell_and_transverse(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset, 
							  Real dx, Real dt, int d){

  Vector<Real> u_iMinus1;
  Vector<Real> u_i;
  Vector<Real> u_iPlus1;
  Vector<Real> u_iPlus2;

  // Slope limiting variable index
  Vector<int> limiting_idx;

  for (int n = NUM_STATE_FLUID; n<NUM_STATE; n++)
    {
      u_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      u_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_iPlus1.push_back(arr(i,j,k,n));
      u_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));

      // each variable limits itself
      limiting_idx.push_back(n-NUM_STATE_FLUID);
    }
  
  // Cell boundary extrapolated values at the left and the right for the ion
  Vector<Real> u_iL = TVD2_reconstruction_L(u_iMinus1, u_i, u_iPlus1, limiting_idx);
  Vector<Real> u_iR = TVD2_reconstruction_R(u_iMinus1, u_i, u_iPlus1, limiting_idx);
  Vector<Real> u_iPlus1L = TVD2_reconstruction_L(u_i, u_iPlus1, u_iPlus2, limiting_idx);
  Vector<Real> u_iPlus1R = TVD2_reconstruction_R(u_i, u_iPlus1, u_iPlus2, limiting_idx);

  // Transverse components
  Vector<Real> trans_comp(NUM_STATE_MAXWELL,0.0);

  for (int n=0; n<NUM_STATE_MAXWELL; n++)
    trans_comp[n] = 0.5*(u_iR[n]+u_iPlus1L[n]);

  // Reiamnn problem left state
  Vector<Real> u_i_nPlusHalf_R = half_update_R(u_iL, u_iR, dx, dt, d, &MaxwellFlux); 
  // Reiamnn problem right state
  Vector<Real> u_iPlus1_nPlusHalf_L = half_update_L(u_iPlus1L, u_iPlus1R, dx, dt, d, &MaxwellFlux);

  // Second order in space, not in time
  // Not fully MUSCL-Hancock
  //u_i = u_iR;
  //u_iPlus1 = u_iPlus1L;
  
  if (CAMReXmp::RKOrder==2)
    {
      u_i = u_i_nPlusHalf_R;
      u_iPlus1 = u_iPlus1_nPlusHalf_L;
    }

  Vector<Real> flux(NUM_STATE_MAXWELL,0.0);
  
  Real B_y_star = 0.5*(u_iPlus1[BX_LOCAL+(1+d)%3]+u_i[BX_LOCAL+(1+d)%3])
    +0.5/c*(u_iPlus1[EX_LOCAL+(2+d)%3]-u_i[EX_LOCAL+(2+d)%3]);
  Real B_z_star = 0.5*(u_iPlus1[BX_LOCAL+(2+d)%3]+u_i[BX_LOCAL+(2+d)%3])
    -0.5/c*(u_iPlus1[EX_LOCAL+(1+d)%3]-u_i[EX_LOCAL+(1+d)%3]);
  Real E_y_star = 0.5*(u_iPlus1[EX_LOCAL+(1+d)%3]+u_i[EX_LOCAL+(1+d)%3])
    -0.5*c*(u_iPlus1[BX_LOCAL+(2+d)%3]-u_i[BX_LOCAL+(2+d)%3]);
  Real E_z_star = 0.5*(u_iPlus1[EX_LOCAL+(2+d)%3]+u_i[EX_LOCAL+(2+d)%3])
    +0.5*c*(u_iPlus1[BX_LOCAL+(1+d)%3]-u_i[BX_LOCAL+(1+d)%3]);  
  
  // EM HLLC states
  flux[BX_LOCAL+d] = 0.0;
  flux[BX_LOCAL+(1+d)%3] = -E_z_star;
  flux[BX_LOCAL+(2+d)%3] = E_y_star;
  flux[EX_LOCAL+d] = 0.0;
  flux[EX_LOCAL+(1+d)%3] = c*c*B_z_star;
  flux[EX_LOCAL+(2+d)%3] =  -c*c*B_y_star;

  std::array<Vector<Real>, 2> final = {flux, trans_comp};
  
  return final;
}

Vector<Real> Maxwell_transverse_comp(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset){

  Vector<Real> u_iMinus1;
  Vector<Real> u_i;
  Vector<Real> u_iPlus1;
  Vector<Real> u_iPlus2;

  // Slope limiting variable index
  Vector<int> limiting_idx;

  int offset = (arr.nComp()==NUM_STATE_MAXWELL) ? 0 : NUM_STATE_FLUID;
  
  for (int n = offset; n<offset+NUM_STATE_MAXWELL; n++)
    {
      u_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      u_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_iPlus1.push_back(arr(i,j,k,n));
      u_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));      

      // each variable limits itself
      limiting_idx.push_back(n-offset);
    }
  
  // Cell boundary extrapolated values at the right on current cell i
  Vector<Real> u_iR = TVD2_reconstruction_R(u_iMinus1, u_i, u_iPlus1, limiting_idx);
  // Cell boundary extrapolated values at the left on next cell i+1
  Vector<Real> u_iPlus1L = TVD2_reconstruction_L(u_i, u_iPlus1, u_iPlus2, limiting_idx);
  
  Vector<Real> trans_comp(NUM_STATE_MAXWELL,0.0);

  for (int n=0; n<NUM_STATE_MAXWELL; n++)
    trans_comp[n] = 0.5*(u_iR[n]+u_iPlus1L[n]);
  
  /* First order
  Vector<Real> trans_comp(NUM_STATE_MAXWELL,0.0);
  for (int n=0; n<NUM_STATE_MAXWELL; n++)
    trans_comp[n] = 0.5*(u_i[n]+u_iPlus1[n]);
  */
  
  return trans_comp;
}
Vector<Real> Maxwell_corner(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset, Real dx, Real dt, int d){

  Vector<Real> u_iMinus1;
  Vector<Real> u_i;
  Vector<Real> u_iPlus1;
  Vector<Real> u_iPlus2;

  // Slope limiting variable index
  Vector<int> limiting_idx;

  int offset = (arr.nComp()==NUM_STATE_MAXWELL) ? 0 : NUM_STATE_FLUID;
  
  for (int n = offset; n<offset+NUM_STATE_MAXWELL; n++)
    {
      u_iMinus1.push_back(arr(i-2*iOffset,j-2*jOffset,k-2*kOffset,n));
      u_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      u_iPlus1.push_back(arr(i,j,k,n));
      u_iPlus2.push_back(arr(i+iOffset,j+jOffset,k+kOffset,n));

      // each variable limits itself
      limiting_idx.push_back(n-offset);
    }
  
  // Cell boundary extrapolated values at the right on current cell i
  Vector<Real> u_iR = TVD2_reconstruction_R(u_iMinus1, u_i, u_iPlus1, limiting_idx);
  // Cell boundary extrapolated values at the left on next cell i+1
  Vector<Real> u_iPlus1L = TVD2_reconstruction_L(u_i, u_iPlus1, u_iPlus2, limiting_idx);
  /*
  // Second order in space and time
  // Fully MUSCL-Hancock
  Vector<Real> u_iL = TVD2_reconstruction_L(u_iMinus1, u_i, u_iPlus1, limiting_idx);
  Vector<Real> u_iPlus1R = TVD2_reconstruction_R(u_i, u_iPlus1, u_iPlus2, limiting_idx);
  // Reiamnn problem left state
  Vector<Real> u_i_nPlusHalf_R = half_update_R(u_iL, u_iR, dx, dt, d, &MaxwellFlux); 
  // Reiamnn problem right state
  Vector<Real> u_iPlus1_nPlusHalf_L = half_update_L(u_iPlus1L, u_iPlus1R, dx, dt, d, &MaxwellFlux);  

  u_iR = u_i_nPlusHalf_R;
  u_iPlus1L = u_iPlus1_nPlusHalf_L;
  */
  Vector<Real> corner(NUM_STATE_MAXWELL,0.0);
  corner[BX_LOCAL+(1+d)%3] = 0.5*(u_iR[BX_LOCAL+(1+d)%3]+u_iPlus1L[BX_LOCAL+(1+d)%3]);
  corner[BX_LOCAL+(2+d)%3] = 0.5*(u_iR[BX_LOCAL+(2+d)%3]+u_iPlus1L[BX_LOCAL+(2+d)%3]);
  corner[EX_LOCAL+(1+d)%3] = 0.5*(u_iR[EX_LOCAL+(1+d)%3]+u_iPlus1L[EX_LOCAL+(1+d)%3]);
  corner[EX_LOCAL+(2+d)%3] = 0.5*(u_iR[EX_LOCAL+(2+d)%3]+u_iPlus1L[EX_LOCAL+(2+d)%3]);

  // Assumes 2D code
  corner[BZ_LOCAL] *= 0.5;
  corner[EZ_LOCAL] *= 0.5;
  
  corner[BX_LOCAL+(1+d)%3] += 0.5/c*(u_iPlus1L[EX_LOCAL+(2+d)%3]-u_iR[EX_LOCAL+(2+d)%3]);
  corner[BX_LOCAL+(2+d)%3] -= 0.5/c*(u_iPlus1L[EX_LOCAL+(1+d)%3]-u_iR[EX_LOCAL+(1+d)%3]);  
  corner[EX_LOCAL+(1+d)%3] -= 0.5*c*(u_iPlus1L[BX_LOCAL+(2+d)%3]-u_iR[BX_LOCAL+(2+d)%3]);
  corner[EX_LOCAL+(2+d)%3] += 0.5*c*(u_iPlus1L[BX_LOCAL+(1+d)%3]-u_iR[BX_LOCAL+(1+d)%3]);
  
  /*// First order
  Vector<Real> corner(NUM_STATE_MAXWELL,0.0);
  corner[BX_LOCAL+(1+d)%3] = 0.5*(u_i[BX_LOCAL+(1+d)%3]+u_iPlus1[BX_LOCAL+(1+d)%3]);
  corner[BX_LOCAL+(2+d)%3] = 0.5*(u_i[BX_LOCAL+(2+d)%3]+u_iPlus1[BX_LOCAL+(2+d)%3]);
  corner[EX_LOCAL+(1+d)%3] = 0.5*(u_i[EX_LOCAL+(1+d)%3]+u_iPlus1[EX_LOCAL+(1+d)%3]);
  corner[EX_LOCAL+(2+d)%3] = 0.5*(u_i[EX_LOCAL+(2+d)%3]+u_iPlus1[EX_LOCAL+(2+d)%3]);

  corner[BZ_LOCAL] *= 0.5;
  corner[EZ_LOCAL] *= 0.5;
  
  //corner[BX_LOCAL+d] += 0.0;
  corner[BX_LOCAL+(1+d)%3] += 0.5/c*(u_iPlus1[EX_LOCAL+(2+d)%3]-u_i[EX_LOCAL+(2+d)%3]);
  corner[BX_LOCAL+(2+d)%3] -= 0.5/c*(u_iPlus1[EX_LOCAL+(1+d)%3]-u_i[EX_LOCAL+(1+d)%3]);
  //corner[EX_LOCAL+d] += 0.0;
  corner[EX_LOCAL+(1+d)%3] -= 0.5*c*(u_iPlus1[BX_LOCAL+(2+d)%3]-u_i[BX_LOCAL+(2+d)%3]);
  corner[EX_LOCAL+(2+d)%3] += 0.5*c*(u_iPlus1[BX_LOCAL+(1+d)%3]-u_i[BX_LOCAL+(1+d)%3]);
  */
  
  return corner;
}
Vector<Real> EM_linearFunc(const Array4<Real>& Bc, const Array4<Real>& Ec,  int i, int j, int k, Real x, Real y, Real z, const Real* dx)
{
  // For second order, there are 27 coeff.
  // i.e. a0,ax,ay,az,axx,...,b0,bx,...,c0,cx,...,czz  
  int a0=0,ax=1,ay=2,az=3,axx=4,axy=5,axz=6;
  int b0=7,bx=8,by=9,bz=10,byy=11,bxy=12,byz=13;
  int c0=14,cx=15,cy=16,cz=17,czz=18,cxz=19,cyz=20;
  
  Vector<Real> EM(NUM_STATE_MAXWELL,0.0);
  
  EM[BX_LOCAL] = Bc(i,j,k,a0) + Bc(i,j,k,ax)*(x/dx[0])
    + Bc(i,j,k,ay)*(y/dx[1]) // + az*
    + Bc(i,j,k,axx)*((x/dx[0])*(x/dx[0]) - 1.0/12.0)
    + Bc(i,j,k,axy)*(x/dx[0])*(y/dx[1]); // + axz*
  EM[BY_LOCAL] = Bc(i,j,k,b0) + Bc(i,j,k,bx)*(x/dx[0])
    + Bc(i,j,k,by)*(y/dx[1]) // + bz*
    + Bc(i,j,k,byy)*((y/dx[1])*(y/dx[1]) - 1.0/12.0)
    + Bc(i,j,k,bxy)*(x/dx[0])*(y/dx[1]); // + bxz*
  EM[BZ_LOCAL] = Bc(i,j,k,c0) + Bc(i,j,k,cx)*(x/dx[0])
    + Bc(i,j,k,cy)*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

  EM[EX_LOCAL] = Ec(i,j,k,a0) + Ec(i,j,k,ax)*(x/dx[0])
    + Ec(i,j,k,ay)*(y/dx[1]) // + az*
    + Ec(i,j,k,axx)*((x/dx[0])*(x/dx[0]) - 1.0/12.0)
    + Ec(i,j,k,axy)*(x/dx[0])*(y/dx[1]); // + axz*
  EM[EY_LOCAL] = Ec(i,j,k,b0) + Ec(i,j,k,bx)*(x/dx[0])
    + Ec(i,j,k,by)*(y/dx[1]) // + bz*
    + Ec(i,j,k,byy)*((y/dx[1])*(y/dx[1]) - 1.0/12.0)
    + Ec(i,j,k,bxy)*(x/dx[0])*(y/dx[1]); // + bxz*
  EM[EZ_LOCAL] = Ec(i,j,k,c0) + Ec(i,j,k,cx)*(x/dx[0])
    + Ec(i,j,k,cy)*(y/dx[1]); // + cz* + cxz* + cyz* + czz		  

  return EM;
  
}
