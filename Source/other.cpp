// this are different limiter functions used from Laguna et al. and Abgrall et at.
// not very useful

Real minmod(Real r_iPlu1, Real r_i)
{
  // MinMod
  // same sign
  if (r_iPlu1*r_i >= 0.0)
    return (r_iPlu1>=0.0 ? 1.0 : -1.0)*std::min(std::abs(r_iPlu1),std::abs(r_i));
  else
    return 0.0;
}
/*Real get_r(const Real& q_i, Real& q_iPlus1, const Real& dx){
  return (q_iPlus1-q_i)/dx;
  }*/
Real u_recons_without_limit(Real u_iMinus1, Real u_i, Real u_iPlus1, Real dx, Real dir)
{
  return u_i + computeDerivative(u_iMinus1, u_iPlus1, dx) * 0.5*dx*dir;
}
Real venkatakrishnan(Real a, Real b, Real dx)
{
  b = (b>=0.0 ? 1.0 : -1.0)*(std::abs(b) + 1e-12);
  Real e = std::sqrt(3.0*dx*dx*dx);
  return 1.0/b * ((a*a + e*e)*b + 2.0*b*b*a)/(a*a + 2.0*b*b + b*a + e*e);
}
Real get_phi(Real u_iMinus1, Real u_i, Real u_iPlus1, Real dx, Real dir)
{
  Real grad_u = computeDerivative(u_iMinus1, u_iPlus1, dx);
  if (grad_u*dir > 0.0)
    return venkatakrishnan(std::max(u_i, std::max(u_iMinus1,u_iPlus1))-u_i, grad_u*0.5*dx*dir, dx);
  else if (grad_u*dir < 0.0)
    return venkatakrishnan(std::min(u_i, std::min(u_iMinus1,u_iPlus1))-u_i, grad_u*0.5*dx*dir, dx);
  else
    return 1.0;
}
Vector<Real> half_step_evolution_L(Vector<Real> u_iMinus1, Vector<Real> u_i, Vector<Real> u_iPlus1,
				   Vector<int> limiting_idx, Real dx, Real dt, int d){
  
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

    Vector<Real> func_iR = func(u_iR, d);
    Vector<Real> func_iL = func(u_iL, d);
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
				   Vector<int> limiting_idx, Real dx, Real dt, int d){
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

    Vector<Real> func_iR = func(u_iR, d);
    Vector<Real> func_iL = func(u_iL, d);
    Vector<Real> diff;
    Vector<Real> half_step_evolution;
    for (int i = 0; i<u_i.size(); i++)
      {
	diff.push_back(func_iR[i]-func_iL[i]);
	half_step_evolution.push_back(u_iR[i] - 0.5*(dt/dx)*diff[i]);
      }
    return half_step_evolution;
}
