#ifndef _utils_H_
#define _utils_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>

#include <memory>
#include <iostream>

#include "CAMReXmp.H"

// needed for EMwaveTM
#include <complex>


using namespace amrex;

extern std::string test;
extern Real Gamma;
extern Real c_h;
extern Real r_i, r_e;
extern Real m;
extern Real l_r;
extern Real c; 
extern Real lambda_d;
extern Real eta;
extern Real cb;
extern Real ce;

// functions to compute the magnitudes
Real get_magnitude_squared(Real x, Real y, Real z);
Real get_magnitude(Real x, Real y, Real z);
// calculate magnitude velocity from converved variables
Real get_v(const Vector<Real>& u_i);
// calculate squared magnitude magnetic field from converved variables
Real get_B_squared(const Vector<Real>& u_i);
// calculate magnitude magnetic field from converved variables
Real get_B(const Vector<Real>& u_i);

Real get_energy(const Vector<Real>& w_i);
Real get_specific_energy(const Vector<Real>& u_i);
Real get_pressure(const Vector<Real>& u_i);
Real get_pressure_total(const Vector<Real>& u_i);
Real get_speed(Vector<Real> u_i);
Real get_speed_e(Vector<Real> u_i);
// Note, B_i is the magnetic field in the i direction
// e.g. for flux update in x direction, use B_x
Real get_speed_a(Real rho, Real B_i);
Real get_speed_f(const Vector<Real>& u_i, int d);

//Vector<Real> func(const Vector<Real>& u_i, int d);
Vector<Real> fluidFlux(const Vector<Real>& u_i, int d);
Vector<Real> MaxwellFlux(const Vector<Real>& u_i, int d);
Real v_exact_dir(Real B_squared, Real B, Real B_three_halves, Real EdotB, Real v_dotB,
                 Real EcrossB_i, Real v_crossB_i, Real Bcross_EcrossB_i, Real Bcross_v_crossB_i,
                 Real B_i, Real t, Real factor);
Vector<Real> cross_product(Vector<Real>, Vector<Real>);
Real dot_product(Vector<Real>, Vector<Real>);
Real dotProduct(Vector<Real> a, Vector<Real> b);

// select the method for slope limiter
// epsilon is the slope limiter
// r is the slope ratio
Real get_epsilon(Real r);
// slope ratio for slope limiter defined in Toro's book
Real get_r(const Real& q_iMinus1, const Real& q_i, const Real& q_iPlus1);
// measure of the slope in a linear reconstruction
Vector<Real> get_delta_i(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1);
// half time step evolution for left state in a linear reconstruction
// define this function for left and right state
Vector<Real> half_step_evolution_L(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1,
				   Vector<int> limiting_idx, Real dx, Real dt, int d,
				   std::function<Vector<Real> (const Vector<Real>&, int)> flux_function);
// half time step evolution for right state in a linear reconstruction
Vector<Real> half_step_evolution_R(Vector<int> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1,
				   Vector<Real> limiting_idx, Real dx, Real dt, int d,
				   std::function<Vector<Real> (const Vector<Real>&, int)> flux_function);
Vector<Real> flux_LF(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, double dt, int d);
Vector<Real> flux_RI(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, double dt, int d);
Vector<Real> flux_FORCE(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
			double dx, double dt, int d);
Vector<Real> flux_SLIC(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset,
                       Real dx, Real dt, int d);
// S* in Riemann problem
Real get_S_star(const Vector<Real>& u_L, const Vector<Real>& u_R,
                Real S_L, Real S_R, int d);
// HLLC left or right states
Vector<Real> get_u_HLLC(const Vector<Real>& u_K, Real S_K, Real S_star, int d);
// S_R (and S_L) definition
Real get_S_K(const Vector<Real>& u_L, const Vector<Real>& u_R, int d);
Real computeDerivative(Real u_iMinus1, Real u_iPlus1, Real dx);
Real computeSecondDerivative(Real u_iMinus1, Real u_i, Real u_iPlus1, Real dx);
Vector<Real> half_step_evolution_L(Vector<Real> u_iMinus1, Vector<Real> u_i,
                            Vector<Real> u_iPlus1, Real dx, Real dt, int d);
// half time step evolution for right state in a linear reconstruction
Vector<Real> half_step_evolution_R(Vector<Real> u_iMinus1, Vector<Real> u_i,
                            Vector<Real> u_iPlus1, Real dx, Real dt, int d);

Vector<Real> flux_LF(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, Vector<Real> dt, int d);
Vector<Real> flux_RI(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, Vector<Real> dt, int d);
Vector<Real> flux_GFORCE(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
			 Real dx, Vector<Real> dt, int d, Real cfl);
Vector<Real> flux_MUSTA1(Vector<Real>& u_i, Vector<Real>& u_iPlus1, 
			 double dx, double dt, int d, Real cfl);
Vector<Real> flux_SLIC(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset,
                       Real dx, Real dt, int d, Real cfl);
Real calc_mutime(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, int d, Real cfl);
Real calc_mutime_elec(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, int d, Real cfl);
Real calc_mutime_Maxwell(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, int d, Real cfl);


Real get_delta_i(Real u_iMinus1, Real u_i, Real u_iPlus1);
Real TVD_slope(Real u_iMinus1, Real u_i, Real u_iPlus1);
Vector<Real> TVD2_reconstruction_L(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1,
				   Vector<int> limiting_idx);
Vector<Real> TVD2_reconstruction_R(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1,
				   Vector<int> limiting_idx);
Vector<Real> half_update_L(Vector<Real> u_iL, Vector<Real> u_iR, Real dx, Real dt, int d,
			   std::function<Vector<Real> (const Vector<Real>&, int)> flux_function);
Vector<Real> half_update_R(Vector<Real> u_iL, Vector<Real> u_iR, Real dx, Real dt, int d,
			   std::function<Vector<Real> (const Vector<Real>&, int)> flux_function);

Vector<Real> HLLC(Vector<Real> uL, Vector<Real> uR, int d);
Vector<Real> RankineHugoniot(Vector<Real> uL, Vector<Real> uR, int d);
Vector<Real> MUSCL_Hancock_TVD_flux(const Array4<Real>& arr, const Array4<Real>& slopes,
				    int i, int j, int k, int iOffset, int jOffset, int kOffset,
				    int start, int len,
				    Real dx, Real dt, int d,
				    std::function<Vector<Real> (const Vector<Real>&, int)> flux_function,
				    std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver);
Vector<Real> MUSCL_Hancock_WENO_flux(const Array4<Real>& arr, const Array4<Real>& slopes,
				     int i, int j, int k, int iOffset, int jOffset, int kOffset,
				     int start, int len,
				     Real dx, Real dt, int d,
				     std::function<Vector<Real> (const Vector<Real>&, int)> flux_function,
				     std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver);
Vector<Real> monotone_flux(const Array4<Real>& arr,
 			   int i, int j, int k, int iOffset, int jOffset, int kOffset,
 			   int start, int len, int d,
 			   std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver);
Vector<Real> TVD_flux(const Array4<Real>& arr, const Array4<Real>& slopes,
		      int i, int j, int k, int iOffset, int jOffset, int kOffset,
		      int start, int len, int startSlope, int d,
		      std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver);
Vector<Real> WENO_flux(const Array4<Real>& arr, const Array4<Real>& slopes,
		       int i, int j, int k, int iOffset, int jOffset, int kOffset,
		       int start, int len, int d,
		       std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver);
Vector<Real> WENO_flux_flat(const Array4<Real>& arr, const Array4<Real>& slopes, const Array4<Real>& tau,
			    int i, int j, int k, int iOffset, int jOffset, int kOffset,
			    int start, int len, int d,
			    std::function<Vector<Real> (Vector<Real>,Vector<Real>,int)> solver);
Vector<Real> exact_flux(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset,
 			int start, int length, Real dx, Real dt, int d);
Vector<Real> EM_linearFunc(const Array4<Real>& Bc, const Array4<Real>& Ec,  int i, int j, int k, Real x, Real y, Real z,
			   const Real* dx);
Vector<Real> EM_quadraticFunc(const Array4<Real>& Bc, const Array4<Real>& Ec,  int i, int j, int k, Real x, Real y, Real z,
			   const Real* dx);
Vector<Real> get_data_zone(const Array4<Real>& arr, int i, int j, int k, int start, int length);
Vector<Real> get_data_stencil(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset, int n);
Real TVD_slope(Vector<Real> u, Vector<Real> limiter);
std::array<Real, 2> WENO3_slope(Vector<Real> u);
Real WENO3_slopeCross(Vector<Real> u, Vector<Real> slopes);
Vector<Real> get_charge_scaled(Vector<Real> u_i, Vector<Real> u_e);
void WENOcharacteristic(const Array4<Real>& arr, const Array4<Real>& slopes,
			int i, int j, int k, int iOffset, int jOffset, int kOffset,
			int start, int len, int d);
////////////////////////////////////////////////////////////
// 2 fluid plasma with generalized Ohm's law functions
////////////////////////////////////////////////////////////
Vector<Real> fluidGOLFlux(const Vector<Real>& u_i, int d);
Vector<Real> TVD_fluxGOL(const Array4<Real>& arr, const Array4<Real>& slopes,
			 int i, int j, int k, int iOffset, int jOffset, int kOffset,
			 int start, int len, int startSlope, Real dx, Real dt, int d);

Vector<Real> flux_LF(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, double dt, int d,
		     std::function<Vector<Real> (const Vector<Real>&, int)> flux_function);
Vector<Real> flux_RI(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
		     double dx, double dt, int d,
		     std::function<Vector<Real> (const Vector<Real>&, int)> flux_function);
Vector<Real> flux_FORCE(const Vector<Real>& u_i, const Vector<Real>& u_iPlus1, 
			double dx, double dt, int d,
			std::function<Vector<Real> (const Vector<Real>&, int)> flux_function);
Vector<Real> flux_MUSTA1(Vector<Real>& u_i, Vector<Real>& u_iPlus1, 
			 double dx, double dt, int d, Real cfl,
			 std::function<Vector<Real> (const Vector<Real>&, int)> flux_function);
Vector<Real> get_electron_var(const Vector<Real>& u_i);
#endif /*_utils_H_*/
