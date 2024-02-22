/*
This is an exact solver for the two-fluid plasma equations.
It is heavily based on https://github.com/ammarhakim/gkyl
Especifically the function gkylMomentSrcExact() from
gkyl/Updater/MomentSrcCommon.cpp 
 */

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

#define sq(x) ((x) * (x))
#define cube(x) ((x) * (x) * (x))
template <typename T> inline static T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/* Updater for the parallel component of the exact source term.
 *
 * @param q_par Normalized parallel electric field and current along the
 *    background B field direction, i.e., [Ez, Jz0, Jz1, ..., Jzs, ...], where s
 *    denotes a species. The content of q_par will be modified in-place.
 * @param dt Time step.
 * @param wp Plasma frequency for each species.
 */
static void
update_par(Eigen::VectorXd &q_par, const double dt, const Eigen::VectorXd &wp)
{
  unsigned nFluids = wp.size();
  double wp_tot2 = 0.;
  for (unsigned n=0; n < nFluids; ++n)
  {
    wp_tot2 += sq(wp[n]);
  }
  double wp_tot = std::sqrt(wp_tot2);

  // TODO reuse v0, v1
  // eigenvector with w=-w_p at t=0
  Eigen::VectorXd v0 = Eigen::VectorXd::Zero(nFluids + 1);
  v0[0] = wp_tot;
  
  // eigenvector with w=w_p at t=0
  Eigen::VectorXd v1(nFluids + 1);
  v1[0] = 0.;
  v1.tail(nFluids).noalias() = wp;

  double coeff0, coeff1;
  // coeff0 = q_par.dot(v0) / v0.squaredNorm();
  // coeff1 = q_par.dot(v1) / v1.squaredNorm();
  coeff0 = q_par[0] / wp_tot;
  coeff1 = q_par.dot(v1) / wp_tot2;

  double cost = std::cos(wp_tot * dt);
  double sint = std::sin(wp_tot * dt);

  // incremental changes to the eigenvectors
  v0 *= -1.;
  v1 *= -1.;
  v0[0] += wp_tot * cost;
  v1[0] += -wp_tot * sint;
  v0.tail(nFluids).noalias() += wp * sint;
  v1.tail(nFluids).noalias() += wp * cost;

  // accumulate incremental changes
  q_par.noalias() += coeff0 * v0  + coeff1 * v1;
}
/*
 * Finding roots of a polynomial as eigenvalues of its companion matrix.
 *
 * @param coeffs The coefficients of the polynomial from the higher to lower order.
 * @return roots The roots.
 */
static Eigen::VectorXd roots(const std::vector<double> &coeffs)
{
  int N = coeffs.size() - 1;

  // companion matrix of the polynomial
  Eigen::MatrixXd M = Eigen::MatrixXd::Constant(N, N, 0);
  
  for(int n = 1; n < N; ++n){
    M(n, n-1) = 1.;
  }

  for(int n = 0; n < N; ++n){
    M(n, N-1) = -coeffs[N-n] / coeffs.front();
  }

  Eigen::VectorXd roots = M.eigenvalues().real();

  return roots;
}
/* Updater for the perpendicular component of the exact source term.
 *
 * @param q_par Normalized perpendicular electric field and current along the
 *    background B field direction, i.e., [Ex, Ey; Jx0, Jy0, Jx1, Jy1,  ...,
 *    Jxs, Jys, ...], where s denotes a species.
 * @param dt Time step.
 * @param wp Plasma frequency for each species.
 * @param Wc Signed cyclotron frequency for each species.
 * @param q_perp_. Updated vector.
 */
static void
update_perp(const Eigen::VectorXd &q_perp, const double dt,
            const Eigen::VectorXd &wp, const Eigen::VectorXd &Wc,
            Eigen::VectorXd &q_perp_)
{
  unsigned nFluids = wp.size();

  // TODO reuse v0, v1
  // compute all eigenvalues
  Eigen::VectorXd eigs(nFluids + 1);
  if (nFluids == 2)
  {
    if (false) {
      std::vector<double> poly_coeff = {
        1.,

        Wc[0] + Wc[1],

        Wc[0] * Wc[1] - sq(wp[0]) - sq(wp[1]),

        -Wc[0] * sq(wp[1]) - Wc[1] * sq(wp[0])
      };
      eigs = roots(poly_coeff);
    } else {
      // analytic solution based on
      // http://web.cs.iastate.edu/~cs577/handouts/polyroots.pdf
      double p = Wc[0] + Wc[1];
      double q = Wc[0] * Wc[1] - sq(wp[0]) - sq(wp[1]);
      double r = -Wc[0] * sq(wp[1]) - Wc[1] * sq(wp[0]);

      double a = (3 * q - sq(p)) / 3;
      double b = (2 * cube(p) - 9 * p * q + 27 * r) / 27;

      double det = sq(b) / 4 + cube(a) / 27;

      if (det < 0) {
        double tmp = 2 * std::sqrt(-a / 3);
        double phi = std::acos(-sgn(b) * std::sqrt(b * b / 4 / (-cube(a) / 27)));
        eigs[0] = tmp * std::cos((phi) / 3) - p / 3;
        eigs[1] = tmp * std::cos((phi + 2 * M_PI) / 3) - p / 3;
        eigs[2] = tmp * std::cos((phi + 4 * M_PI) / 3) - p / 3;
      } else if (det == 0) {
        double tmp = sgn(b) * std::sqrt(-a / 3);
        eigs[0] = -2 * tmp;
        eigs[1] = tmp;
        eigs[2] = tmp;
      } else {
        assert(false);
      }
    }
  } else if (nFluids == 3) {
    std::vector<double> poly_coeff = {
      1.,
      
      Wc[0] + Wc[1] + Wc[2],
      
      Wc[0] * Wc[1] + Wc[0] * Wc[2] +
      Wc[1] * Wc[2] - sq(wp[0]) - sq(wp[1]) - sq(wp[2]),

      Wc[0] * Wc[1] * Wc[2] - Wc[0] * sq(wp[1]) - Wc[0] * sq(wp[2]) -
      Wc[1] * sq(wp[0]) - Wc[1] * sq(wp[2]) - Wc[2] * sq(wp[0]) -
      Wc[2] * sq(wp[1]),
      
      -Wc[0] * Wc[1] * sq(wp[2]) -
      Wc[0] * Wc[2] * sq(wp[1]) - Wc[1] * Wc[2] * sq(wp[0])
    };
    eigs = roots(poly_coeff);
  } else {
    assert(false);
  }

  // compute the two eigenvector for each eigenvalue and accumulate their
  // contributions to the final parallel state vector
  q_perp_.setZero();
  Eigen::VectorXd v0(2 * (nFluids + 1));
  Eigen::VectorXd v1(2 * (nFluids + 1));
  for (unsigned i = 0; i < nFluids + 1; ++i)
  {
    double w = eigs[i];

    // compute the two eigenvectors for w at t=0
    v0[0] = 0.;
    v0[1] = 1.;
    v1[0] = 1.;
    v1[1] = 0.;
    for (unsigned n = 0; n < nFluids; ++n)
    {
      unsigned nn = n + 1;
      double tmp = wp[n] / (w + Wc[n]);
      v0[2 * nn] = tmp;
      v0[2 * nn + 1] = 0.;
      v1[2 * nn] = 0;
      v1[2 * nn + 1] = -tmp;
    }

    // compute eigencoefficients
    double coeff0, coeff1;
    coeff0 = q_perp.dot(v0) / v0.squaredNorm();
    coeff1 = q_perp.dot(v1) / v1.squaredNorm();

    // compute the two eigenvectors for w at t=dt
    double cost = std::cos(w * dt);
    double sint = std::sin(w * dt);
    v0[0] = -sint;
    v0[1] = cost;
    v1[0] = cost;
    v1[1] = sint;
    for (unsigned n = 0; n < nFluids; ++n)
    {
      unsigned nn = n + 1;
      double tmp = wp[n] / (w + Wc[n]);
      v0[2 * nn] = tmp * cost;
      v0[2 * nn + 1] = tmp * sint;
      v1[2 * nn] = tmp * sint;
      v1[2 * nn + 1] = -tmp * cost;
    }

    // accumulate contribution from the two eigenvectors
    q_perp_.noalias() += coeff0 * v0  + coeff1 * v1;
  }
}
/* Rotate around axis by angle.
 *
 * @param v Original vector.
 * @param cosa cos(angle).
 * @param cosah cos(angle/2).
 * @param sinah sin(angle/2).
 * @param axis Normalized axis vector.
 *
 * @return Rotated vector.
 */
static inline Eigen::Vector3d
rotate(const Eigen::Vector3d &v, const double cosa, const double cosah,
       const double sinah, const Eigen::Vector3d &axis)
{
  return v * cosa + 2. * sq(sinah) * (v.dot(axis)) * axis
       + 2. * cosah * sinah * v.cross(axis);
}
void CAMReXmp::sourceUpdateExact(Array4<Real>& arr, int i, int j, int k, Real dt)
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

  // get kinetic energies
  Real v_i = get_magnitude(v_x_i, v_y_i, v_z_i);
  Real v_e = get_magnitude(v_x_e, v_y_e, v_z_e);
  Real kin_i = 0.5*rho_i*v_i*v_i;
  Real kin_e = 0.5*rho_e*v_e*v_e;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Eigen::Vector3d B(B_x, B_y, B_z);
  Real Bmag = B.norm();

  Eigen::Vector3d E(E_x, E_y, E_z);
  Real Enorm = 1.; // nominal normalization
  Eigen::Vector3d E_ = E / Enorm;

  int nFluids = 2;
  Eigen::VectorXd qbym(nFluids); 
  Eigen::VectorXd Wc(nFluids); 
  Eigen::VectorXd wp(nFluids);
  Eigen::VectorXd Pnorm(nFluids); 
  std::vector<Eigen::Vector3d> J_(nFluids);

  qbym[0] = r_i/l_r, qbym[1] = r_e/l_r;
  Real epsilon0 = lambda_d*lambda_d;
  for (unsigned n=0; n < nFluids; ++n)
  {
    //double *f = ff[n];

    //qbym[n] = fd[n].charge / fd[n].mass;
    Wc[n] = qbym[n] * Bmag;
    double wp2 = arr(i,j,k,RHO_I+(NUM_STATE_FLUID/2)*n) * sq(qbym[n]) / epsilon0;
    wp[n] = std::sqrt(wp2);

    Pnorm[n] = std::sqrt(epsilon0 * arr(i,j,k,RHO_I+(NUM_STATE_FLUID/2)*n));
    //if (fd[n].charge < 0.)
    if (n==1)
    {
      Pnorm[n] *= -1.;
    }
    J_[n][0] = arr(i,j,k,MOMX_I+(NUM_STATE_FLUID/2)*n) / Pnorm[n];
    J_[n][1] = arr(i,j,k,MOMY_I+(NUM_STATE_FLUID/2)*n) / Pnorm[n];
    J_[n][2] = arr(i,j,k,MOMZ_I+(NUM_STATE_FLUID/2)*n) / Pnorm[n];
  }
  
  
  double chargeDens = 0.0;
  for (unsigned n = 0; n < nFluids; ++n)
  {
    //double *f = ff[n];
    chargeDens += qbym[n] * arr(i,j,k,RHO_I+(NUM_STATE_FLUID/2)*n);
  }
  
  //if (Bmag > 0.)  // TODO set threshold
  {
    double angle = std::acos(B[2] / Bmag);
    Eigen::Vector3d axis = B.cross(Eigen::Vector3d::UnitZ());
    axis.normalize();  // AngleAxisd needs a normalized axis vector
    double cosa = std::cos(angle);
    double cosah = std::cos(angle / 2.);
    double sinah = std::sin(angle / 2.);
    
    E_ = rotate(E_, cosa, cosah, -sinah, axis);
    for (unsigned n=0; n < nFluids; ++n)
      J_[n] = rotate(J_[n], cosa, cosah, -sinah, axis);
    
    // parallel component of initial condition
    Eigen::VectorXd q_par(nFluids + 1);
    q_par[0] = E_[2];
    for (unsigned n=0; n < nFluids; ++n)
    {
      q_par[n + 1] = J_[n][2];
    }

    // perpendicular component of initial condition
    Eigen::VectorXd q_perp(2 * (nFluids + 1));
    q_perp[0] = E_[0];
    q_perp[1] = E_[1];
    for (unsigned n=0; n < nFluids; ++n)
    {
      q_perp[2*(n + 1)] = J_[n][0];
      q_perp[2*(n + 1) + 1] = J_[n][1];
    }
    
    // update parallel component
    update_par(q_par, dt, wp);
    // update perpendicular component
    Eigen::VectorXd q_perp_(2 * (nFluids + 1));
    update_perp(q_perp, dt, wp, Wc, q_perp_);

    // fill the full vectors
    E_[2] = q_par[0];
    for (unsigned n = 0; n < nFluids; ++n)
    {
      J_[n][2] = q_par[n + 1];
    }
    E_[0] = q_perp_[0];
    E_[1] = q_perp_[1];
    for (unsigned n = 0; n < nFluids; ++n)
    {
      unsigned nn = n + 1;
      J_[n][0] = q_perp_[2 * nn];
      J_[n][1] = q_perp_[2 * nn + 1];
    }

    // rotate back
    E_ = rotate(E_, cosa, cosah, sinah, axis);
    for (unsigned n=0; n < nFluids; ++n)
      J_[n] = rotate(J_[n], cosa, cosah, sinah, axis);

    // fill state vector
    arr(i,j,k,EX) = E_[0] * Enorm;
    arr(i,j,k,EY) = E_[1] * Enorm;
    arr(i,j,k,EZ) = E_[2] * Enorm;
    for (unsigned n = 0; n < nFluids; ++n)
    {
      //double *f = ff[n];
      arr(i,j,k,MOMX_I+(NUM_STATE_FLUID/2)*n) = J_[n][0] * Pnorm[n];
      arr(i,j,k,MOMY_I+(NUM_STATE_FLUID/2)*n) = J_[n][1] * Pnorm[n];
      arr(i,j,k,MOMZ_I+(NUM_STATE_FLUID/2)*n) = J_[n][2] * Pnorm[n];
    }
    
  } /*else {
    Eigen::VectorXd q_par(nFluids + 1);

    for (unsigned d = 0; d < 3; ++d)
    {
      q_par[0] = E_[d];
      for (unsigned n=0; n < nFluids; ++n)
      {
        q_par[n + 1] = J_[n][d];
      }

      // FIXME avoid repeated calculation of eigenvectors
      update_par(q_par, dt, wp);

      // re-normalize back
      em[EX + d] = q_par[0] * Enorm;
      for (unsigned n = 0; n < nFluids; ++n)
      {
        double *f = ff[n];
        f[MX + d] = q_par[n + 1] * Pnorm[n];
      }
    }
  }
    */
  /*
  if (sd->hasSigma) {
    em[EX] = em[EX]*std::exp(-sigma[0]*dt/sd->epsilon0);
    em[EY] = em[EY]*std::exp(-sigma[0]*dt/sd->epsilon0);
    em[EZ] = em[EZ]*std::exp(-sigma[0]*dt/sd->epsilon0);
  }

  //------------> update correction potential
  double crhoc = sd->chi_e * chargeDens/epsilon0;
  em[PHIE] += dt * crhoc;
  */  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // get new kinetic energies
  Real mom_i_new = get_magnitude(arr(i,j,k,MOMX_I), arr(i,j,k,MOMY_I), arr(i,j,k,MOMZ_I));
  Real mom_e_new = get_magnitude(arr(i,j,k,MOMX_E), arr(i,j,k,MOMY_E), arr(i,j,k,MOMZ_E));
  Real kin_i_new = 0.5*mom_i_new*mom_i_new/rho_i;
  Real kin_e_new = 0.5*mom_e_new*mom_e_new/rho_e;

  arr(i,j,k,ENER_I) += kin_i_new-kin_i;
  arr(i,j,k,ENER_E) += kin_e_new-kin_e;
  
}
