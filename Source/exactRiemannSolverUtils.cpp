#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#include "utils.H"

Real pressure(Real rho, Real momentum, Real E){
    // define velocity and pressure
    Real v = momentum/rho;
    // absolute value
    return (Gamma-1.0)*(E-(rho*v)*(rho*v)/(2.0*rho));
}
Real energy(Real rho, Real v, Real p){
    return p/(Gamma-1.0) + (rho*v)*(rho*v)/(2.0*rho);
}
Real speed(Real rho, Real p){
  return std::sqrt((Gamma*p)/rho);
}
Real f_K(Real rho, Real p, Real p_star){
    Real A = 2.0/((Gamma+1.0)*rho);
    Real B = (Gamma-1.0)/(Gamma+1.0)*p;
    // shock
    if (p_star>p){
        return (p_star-p)*std::sqrt(A/(p_star+B));
    // rarefaction
    } else{
        return 2*speed(rho, p)/(Gamma-1.0)*(std::pow(p_star/p, (Gamma-1.0)/(2.0*Gamma))-1.0);
    }
}
Real f(Real rho_L, Real v_L, Real p_L,
        Real rho_R,  Real v_R, Real p_R, Real p_star){
    return f_K(rho_L, p_L, p_star) + f_K(rho_R, p_R, p_star) + (v_R-v_L);
}
Real f_prime_K(Real rho, Real p, Real p_star){
    Real A = 2.0/((Gamma+1.0)*rho);
    Real B = (Gamma-1.0)/(Gamma+1.0)*p;
    // shock
    if (p_star>p){
        return std::pow(A/(B+p_star), 0.5)*(1.0-(p_star-p)/(2*(B+p_star)));
    // rarefaction
    } else{
        return 1.0/(rho*speed(rho, p))*std::pow(p_star/p, -(Gamma+1.0)/(2.0*Gamma));
    }
}
Real f_prime(Real rho_L, Real p_L, 
        Real rho_R, Real p_R, Real p_star){
    return f_prime_K(rho_L, p_L, p_star) + f_prime_K(rho_R, p_R, p_star);
}
Real get_p_star(const Vector<Real>& w_L,
		const Vector<Real>& w_R, 
		int d, Real tolerance = 1e-6){
  Real rho_L = w_L[0], v_L = w_L[1+d]/w_L[0], p_L = get_pressure(w_L);
  Real rho_R = w_R[0], v_R = w_R[1+d]/w_R[0], p_R = get_pressure(w_R);
  Real p_old, p_new = 0.5*(p_L+p_R);
  do{
    p_old = p_new;
    p_new = p_old - f(rho_L, v_L, p_L, rho_R, v_R, p_R, p_old)/
      f_prime(rho_L, p_L, rho_R, p_R, p_old);
  } while (std::abs(p_old-p_new)/p_old > tolerance);
  return p_new;
}
Real get_rho_star(const Vector<Real>& w, Real p_star){
  Real rho = w[0], p = get_pressure(w);
    // shock
    if (p_star>p){
        return rho * (((p_star/p)+(Gamma-1.0)/(Gamma+1.0))/((Gamma-1.0)/(Gamma+1.0)*(p_star/p)+1.0));
    // rarefaction
    } else{
        return rho * std::pow(p_star/p, 1.0/Gamma);
    }
}
Real get_S_L(const Vector<Real>& w, Real p_star, int d){
    Real rho = w[0], v = w[1+d]/w[0], p = get_pressure(w);
    return v-speed(rho, p)*std::sqrt(((Gamma+1.0)/(2.0*Gamma)*p_star/p+(Gamma-1.0)/(2.0*Gamma)));
}
Real get_S_R(const Vector<Real>& w, Real p_star, int d){
  Real rho = w[0], v = w[1+d]/w[0], p = get_pressure(w);
  return v+speed(rho, p)*std::sqrt(((Gamma+1.0)/(2.0*Gamma)*p_star/p+(Gamma-1.0)/(2.0*Gamma)));
}
Real get_speed_star(const Vector<Real>& w, Real p_star){
    Real rho = w[0], p = get_pressure(w);
    Real c_s = speed(rho, p);
    return c_s*std::pow(p_star/p, (Gamma-1.0)/(2.0*Gamma));
}
Real get_v_star(const Vector<Real>& w_L, 
                const Vector<Real>& w_R, Real p_star, int d){  
  Real rho_L = w_L[0], v_L = w_L[1+d]/w_L[0], p_L = get_pressure(w_L);
  Real rho_R = w_R[0], v_R = w_R[1+d]/w_R[0], p_R = get_pressure(w_R);
  return 0.5*(v_L+v_R)+0.5*(f_K(rho_R, p_R, p_star)-
			    f_K(rho_L, p_L, p_star));
}
Vector<Real> get_w_L_fan(const Vector<Real>& w, Real S, int d){
    Real rho = w[0], v = w[1+d]/w[0], p = get_pressure(w);
    Vector<Real> w_fan(5,0.0);
    w_fan[0] = rho*std::pow(2.0/(Gamma+1.0)+(Gamma-1.0)/
                ((Gamma+1.0)*speed(rho, p))*(v-S),
                2.0/(Gamma-1.0));
    w_fan[1+d] = 2.0/(Gamma+1.0)*(speed(rho, p)+
                (Gamma-1.0)/2.0 * v + S);
    w_fan[ENER_I] = p*std::pow(2.0/(Gamma+1.0)+(Gamma-1.0)/
                ((Gamma+1.0)*speed(rho, p))*(v-S),
			       2.0*Gamma/(Gamma-1.0));
    return w_fan;
}
Vector<Real> get_w_R_fan(const Vector<Real>& w, Real S, int d){
    Real rho = w[0], v = w[1+d]/w[0], p = get_pressure(w);
    Vector<Real> w_fan(5,0.0);
    w_fan[0] = rho*std::pow(2.0/(Gamma+1.0)-(Gamma-1.0)/
                ((Gamma+1.0)*speed(rho, p))*(v-S),
                2.0/(Gamma-1.0));
    w_fan[1+d] = 2.0/(Gamma+1.0)*(-speed(rho, p)+
                (Gamma-1.0)/2.0 * v + S);
    w_fan[ENER_I] = p*std::pow(2.0/(Gamma+1.0)-(Gamma-1.0)/
                ((Gamma+1.0)*speed(rho, p))*(v-S),
                2.0*Gamma/(Gamma-1.0));
    return w_fan;
}

Vector<Real> exact_flux(const Array4<Real>& arr, int i, int j, int k, int iOffset, int jOffset, int kOffset,
 			int start, int length, Real dx, Real dt, int d){
  
  Vector<Real> w_L_i;
  Vector<Real> w_R_i;

  for (int n = start; n<start+length; n++)
    {
      w_L_i.push_back(arr(i-iOffset,j-jOffset,k-kOffset,n));
      w_R_i.push_back(arr(i,j,k,n));
    }

  Vector<Real> flux(5,0.0);
  
  Real x = 0.0;
  Real S = x/dt;

  Real p_star_i = get_p_star(w_L_i, w_R_i, d);
  Real rho_L_star_i = get_rho_star(w_L_i, p_star_i);
  Real rho_R_star_i = get_rho_star(w_R_i, p_star_i);
  Real v_star_i = get_v_star(w_L_i, w_R_i, p_star_i, d);
  
  Real S_L_i = get_S_L(w_L_i, p_star_i, d);
  Real S_R_i = get_S_R(w_R_i, p_star_i, d);
    
  Real speed_star_L_i = get_speed_star(w_L_i, p_star_i);
  Real speed_star_R_i = get_speed_star(w_R_i, p_star_i);

  Real S_HL_i = w_L_i[1]/w_L_i[0]-speed(w_L_i[0], w_L_i[2]);
  Real S_TL_i = v_star_i-speed_star_L_i;
  Real S_TR_i = v_star_i+speed_star_R_i;
  Real S_HR_i = w_R_i[1]/w_R_i[0]+speed(w_R_i[0], w_R_i[2]);

  // LEFT
  if (S <= v_star_i){
    // left shock
    if (p_star_i>get_pressure(w_L_i)){
      if (S < S_L_i){
	Vector<Real> function = fluidFlux(w_L_i, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }
	//output << x << " " << w_L_i[0] << " " << w_L_i[1] << " " << w_L_i[2]
	//     << " " << w_L_i[2]/(w_L_i[0]*(gamma-1.0)) << endl;
      } else{
	Vector<Real> w_star(5,0.0);
	w_star[0] = rho_L_star_i;
	w_star[1+d] = rho_L_star_i*v_star_i;
	w_star[2-d] = rho_L_star_i*w_L_i[2-d]/w_L_i[0];
	w_star[3] = rho_L_star_i*w_L_i[3]/w_L_i[0];
	w_star[ENER_I] = get_energy({rho_L_star_i,v_star_i,w_L_i[2-d]/w_L_i[0],w_L_i[3]/w_L_i[0],p_star_i});
	Vector<Real> function = fluidFlux(w_star, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }
	//output << x << " " << rho_L_star_i << " " << v_star << " " << p_star 
	//     << " " << p_star/(rho_L_star_i*(gamma-1.0)) << endl;
      }
      // left rarefaction
    } else{
      if (S < S_HL_i){
	Vector<Real> function = fluidFlux(w_L_i, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }
	//output << x << " " << w_L_i[0] << " " << w_L_i[1] << " " << w_L_i[2] 
	//     << " " << w_L_i[2]/(w_L_i[0]*(gamma-1.0)) << endl;
      } else if (S >= S_HL_i && S <= S_TL_i){
	Vector<Real> w_L_fan_i = get_w_L_fan(w_L_i, S, d);
	
	w_L_fan_i[ENER_I] = get_energy({w_L_fan_i[0],w_L_fan_i[1+d],w_L_i[2-d]/w_L_i[0],
					w_L_i[3]/w_L_i[0],w_L_fan_i[ENER_I]});

	w_L_fan_i[1+d] *= w_L_fan_i[0];
	w_L_fan_i[2-d] = w_L_fan_i[0]*w_L_i[2-d]/w_L_i[0];
	w_L_fan_i[3] = w_L_fan_i[0]*w_L_i[3]/w_L_i[0];

	Vector<Real> function = fluidFlux(w_L_fan_i, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }
	//output << x << " " << w_L_ifan[0] << " " << w_L_ifan[1] << " " << w_L_ifan[2] 
	//     << " " << w_L_ifan[2]/(w_L_ifan[0]*(gamma-1.0)) << endl;
      } else{
	Vector<Real> w_star(5,0.0);
	w_star[0] = rho_L_star_i;
	w_star[1+d] = rho_L_star_i*v_star_i;
	w_star[2-d] = rho_L_star_i*w_L_i[2-d]/w_L_i[0];
	w_star[3] = rho_L_star_i*w_L_i[3]/w_L_i[0];
	w_star[ENER_I] = get_energy({rho_L_star_i,v_star_i,w_L_i[2-d]/w_L_i[0],w_L_i[3]/w_L_i[0],p_star_i});
	Vector<Real> function = fluidFlux(w_star, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }
	
	//output << x << " " << rho_L_star_i << " " << v_star << " " << p_star 
	//     << " " << p_star/(rho_L_star_i*(gamma-1.0)) << endl;
      }
    }
    // RIGHT
        } else {
    // right shock
    if (p_star_i>get_pressure(w_R_i)){
      if (S >= v_star_i && S <= S_R_i){
	Vector<Real> w_star(5,0.0);
	w_star[0] = rho_R_star_i;
	w_star[1+d] = rho_R_star_i*v_star_i;
	w_star[2-d] = rho_R_star_i*w_R_i[2-d]/w_R_i[0];
	w_star[3] = rho_R_star_i*w_R_i[3]/w_R_i[0];
	w_star[ENER_I] = get_energy({rho_R_star_i,v_star_i,w_R_i[2-d]/w_R_i[0],w_R_i[3]/w_R_i[0],p_star_i});
	Vector<Real> function = fluidFlux(w_star, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }
	//output << x << " " << rho_R_star_i << " " << v_star << " " << p_star 
	//     << " " << p_star/(rho_R_star_i*(gamma-1.0)) << endl;
      } else{
	Vector<Real> function = fluidFlux(w_R_i, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }
	//output << x << " " << w_R_i[0] << " " << w_R_i[1] << " " << w_R_i[2] 
	//          << " " << w_R_i[2]/(w_R_i[0]*(gamma-1.0)) << endl;
      }
      // right rarefaction
    } else{
      if (S >= v_star_i && S <= S_TR_i){
	Vector<Real> w_star(5,0.0);
	w_star[0] = rho_R_star_i;
	w_star[1+d] = rho_R_star_i*v_star_i;
	w_star[2-d] = rho_R_star_i*w_R_i[2-d]/w_R_i[0];
	w_star[3] = rho_R_star_i*w_R_i[3]/w_R_i[0];
	w_star[ENER_I] = get_energy({rho_R_star_i,v_star_i,w_R_i[2-d]/w_R_i[0],w_R_i[3]/w_R_i[0],p_star_i});
	Vector<Real> function = fluidFlux(w_star, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }		
	//output << x << " " << rho_R_star_i << " " << v_star << " " << p_star 
	//     << " " << p_star/(rho_R_star_i*(gamma-1.0)) << endl;
      } else if (S >= S_TR_i && S <= S_HR_i){
	Vector<Real> w_R_fan_i = get_w_R_fan(w_R_i, S, d);	
	
	w_R_fan_i[ENER_I] = get_energy({w_R_fan_i[0],w_R_fan_i[1+d],w_R_i[2-d]/w_R_i[0],
				    w_R_i[3]/w_R_i[0],w_R_fan_i[ENER_I]});

	w_R_fan_i[1+d] *= w_R_fan_i[0];
	w_R_fan_i[2-d] = w_R_fan_i[0]*w_R_i[2-d]/w_R_i[0];
	w_R_fan_i[3] = w_R_fan_i[0]*w_R_i[3]/w_R_i[0];

	Vector<Real> function = fluidFlux(w_R_fan_i, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }
	//output << x << " " << w_R_ifan[0] << " " << w_R_ifan[1] << " " << w_R_ifan[2]
	//     << " " << w_R_ifan[2]/(w_R_ifan[0]*(gamma-1.0)) << endl;
      } else{
	Vector<Real> function = fluidFlux(w_R_i, d);
	for (int n = 0; n<=ENER_I; n++)
	  {
	    flux[n] = function[n];
	  }	
	//output << x << " " << w_R_i[0] << " " << w_R_i[1] << " " << w_R_i[2] 
	//     << " " << w_R_i[2]/(w_R_i[0]*(gamma-1.0))<< endl;
      }
    }
  }

    
  return flux;
}
