/*
 * rosenbrock4.cpp
 *
 * Copyright 2010-2012 Mario Mulansky
 * Copyright 2011-2012 Karsten Ahnert
 * Copyright 2012 Andreas Angelopoulos
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <iostream>
#include <fstream>
#include <utility>

#include <boost/numeric/odeint.hpp>

#include <boost/phoenix/core.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;


// indices for the 5-moment source terms
enum StateVariable {
		    VX = 0, VY = 1, VZ = 2, ENER = 3
};
// charge to mass ratio and larmor radius
double r, lr;
// density and EM fields
double rho, Bx, By, Bz, Ex, Ey, Ez;

//[ stiff_system_definition
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct stiff_system
{
    void operator()( const vector_type &x , vector_type &dxdt , double /* t */ )
    {
        dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ];
        dxdt[ 1 ] = x[ 0 ];
    }
};

struct stiff_system_jacobi
{
    void operator()( const vector_type & /* x */ , matrix_type &J , const double & /* t */ , vector_type &dfdt )
    {
        J( 0 , 0 ) = -101.0;
        J( 0 , 1 ) = -100.0;
        J( 1 , 0 ) = 1.0;
        J( 1 , 1 ) = 0.0;
        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
    }
};

struct five_moment_source
{
  void operator()( const vector_type &x , vector_type &dxdt , double /* t */ )
  {
    dxdt[ VX ] = r/lr*(Ex + x[VY]*Bz - x[VZ]*By);
    dxdt[ VY ] = r/lr*(Ey + x[VZ]*Bx - x[VX]*Bz);
    dxdt[ VZ ] = r/lr*(Ez + x[VX]*By - x[VY]*Bx);
    dxdt[ ENER ] = r/lr*rho*(x[VX]*Ex + x[VY]*Ey + x[VZ]*Ez);
  }
};

struct five_moment_source_jacobi
{
  void operator()( const vector_type & x , matrix_type &J , const double & /* t */ , vector_type &dfdt )
  {
    J( VX , VX ) = 0.0;
    J( VX , VY ) = r/lr*Bz;
    J( VX , VZ ) = -r/lr*By;
    J( VX , ENER ) = 0.0;
    J( VY , VX ) = -r/lr*Bz;
    J( VY , VY ) = 0.0;
    J( VY , VZ ) = r/lr*Bx;
    J( VY , ENER ) = 0.0;
    J( VZ , VX ) = r/lr*By;
    J( VZ , VY ) = -r/lr*Bx;
    J( VZ , VZ ) = 0.0;
    J( VZ , ENER ) = 0.0;
    J( ENER , VX ) = r/lr*rho*Ex;
    J( ENER , VY ) = r/lr*rho*Ey;
    J( ENER , VZ ) = r/lr*rho*Ez;
    J( ENER , ENER ) = 0.0;

    dfdt[VX] = 0.0;
    dfdt[VY] = 0.0;
    dfdt[VZ] = 0.0;
    dfdt[ENER] = 0.0;
  }
};

size_t stiffSolver(std::vector<double>& u, double chargeToMassRatio, double larmorRadius, double dt)
{
  r = chargeToMassRatio;
  lr = larmorRadius;
  rho = u[0], Bx = u[5], By = u[6], Bz = u[7], Ex = u[8], Ey = u[9], Ez = u[10];

  vector_type x( 4 , 0.0 );
  x[VX]=u[1], x[VY]=u[2], x[VZ]=u[3], x[ENER]=u[4];
  
  size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ) ,
					 make_pair( five_moment_source() , five_moment_source_jacobi() ) ,
					 x , 0.0 , dt , dt/10.0);
  typedef rosenbrock4< double > stepper_type;
  typedef rosenbrock4_controller< stepper_type > controlled_stepper_type;
  /*
  size_t num_of_steps = integrate_const( controlled_stepper_type() ,
					 make_pair( five_moment_source() , five_moment_source_jacobi() ) ,
					 x , 0.0 , dt , dt/100.0);*/
  u[1] = x[VX], u[2] = x[VY], u[3] = x[VZ], u[4] = x[ENER];
  return num_of_steps;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// indices for the 5-moment source terms
enum StateVariableFull {
			MOMX_I = 0, MOMY_I = 1, MOMZ_I = 2, MOMX_E = 3, MOMY_E = 4, MOMZ_E = 5,
			EX = 6, EY = 7, EZ = 8
};
// charge to mass ratios and Debye length
double ri, re, ld;
// densities and energies
double rho_i, rho_e, E_i, E_e;

struct five_moment_source_full
{
  void operator()( const vector_type &x , vector_type &dxdt , double /* t */ )
  {
    dxdt[ MOMX_I ] = ri/lr*(x[EX]*rho_i + x[MOMY_I]*Bz - x[MOMZ_I]*By);
    dxdt[ MOMY_I ] = ri/lr*(x[EY]*rho_i + x[MOMZ_I]*Bx - x[MOMX_I]*Bz);
    dxdt[ MOMZ_I ] = ri/lr*(x[EZ]*rho_i + x[MOMX_I]*By - x[MOMY_I]*Bx);
    dxdt[ MOMX_E ] = re/lr*(x[EX]*rho_e + x[MOMY_E]*Bz - x[MOMZ_E]*By);
    dxdt[ MOMY_E ] = re/lr*(x[EY]*rho_e + x[MOMZ_E]*Bx - x[MOMX_E]*Bz);
    dxdt[ MOMZ_E ] = re/lr*(x[EZ]*rho_e + x[MOMX_E]*By - x[MOMY_E]*Bx);
    dxdt[  EX    ] = -1.0/(ld*ld*lr)*(ri*x[MOMX_I] + re*x[MOMX_E]);
    dxdt[  EY    ] = -1.0/(ld*ld*lr)*(ri*x[MOMY_I] + re*x[MOMY_E]);
    dxdt[  EZ    ] = -1.0/(ld*ld*lr)*(ri*x[MOMZ_I] + re*x[MOMZ_E]);
  }
};

struct five_moment_source_full_jacobi
{
  void operator()( const vector_type & x , matrix_type &J , const double & /* t */ , vector_type &dfdt )
  {
    J( MOMX_I , MOMX_I ) = 0.0;
    J( MOMX_I , MOMY_I ) = ri/lr*Bz;
    J( MOMX_I , MOMZ_I ) = -ri/lr*By;
    J( MOMX_I , MOMX_E ) = 0.0;
    J( MOMX_I , MOMY_E ) = 0.0;
    J( MOMX_I , MOMZ_E ) = 0.0;    
    J( MOMX_I , EX     ) = ri/lr*rho_i;
    J( MOMX_I , EY     ) = 0.0;
    J( MOMX_I , EZ     ) = 0.0;
    
    J( MOMY_I , MOMX_I ) = -ri/lr*Bz;
    J( MOMY_I , MOMY_I ) = 0.0;
    J( MOMY_I , MOMZ_I ) = ri/lr*Bx;
    J( MOMY_I , MOMX_E ) = 0.0;
    J( MOMY_I , MOMY_E ) = 0.0;
    J( MOMY_I , MOMZ_E ) = 0.0;
    J( MOMY_I , EX     ) = 0.0;
    J( MOMY_I , EY     ) = ri/lr*rho_i;
    J( MOMY_I , EZ     ) = 0.0;
    
    J( MOMZ_I , MOMX_I ) = ri/lr*By;
    J( MOMZ_I , MOMY_I ) = -ri/lr*Bx;
    J( MOMZ_I , MOMZ_I ) = 0.0;
    J( MOMZ_I , MOMX_E ) = 0.0;
    J( MOMZ_I , MOMY_E ) = 0.0;
    J( MOMZ_I , MOMZ_E ) = 0.0;
    J( MOMZ_I , EX     ) = 0.0;
    J( MOMZ_I , EY     ) = 0.0;
    J( MOMZ_I , EZ     ) = ri/lr*rho_i;
    
    J( MOMX_E , MOMX_I ) = 0.0;
    J( MOMX_E , MOMY_I ) = 0.0;
    J( MOMX_E , MOMZ_I ) = 0.0;
    J( MOMX_E , MOMX_E ) = 0.0;
    J( MOMX_E , MOMY_E ) = re/lr*Bz;
    J( MOMX_E , MOMZ_E ) = -re/lr*By;
    J( MOMX_E , EX     ) = re/lr*rho_e;
    J( MOMX_E , EY     ) = 0.0;
    J( MOMX_E , EZ     ) = 0.0;
    
    J( MOMY_E , MOMX_I ) = 0.0;
    J( MOMY_E , MOMY_I ) = 0.0;
    J( MOMY_E , MOMZ_I ) = 0.0;
    J( MOMY_E , MOMX_E ) = -re/lr*Bz;
    J( MOMY_E , MOMY_E ) = 0.0;
    J( MOMY_E , MOMZ_E ) = re/lr*Bx;
    J( MOMY_E , EX     ) = 0.0;
    J( MOMY_E , EY     ) = re/lr*rho_e;
    J( MOMY_E , EZ     ) = 0.0;
    
    J( MOMZ_E , MOMX_I ) = 0.0;
    J( MOMZ_E , MOMY_I ) = 0.0;
    J( MOMZ_E , MOMZ_I ) = 0.0;
    J( MOMZ_E , MOMX_E ) = re/lr*By;
    J( MOMZ_E , MOMY_E ) = -re/lr*Bx;
    J( MOMZ_E , MOMZ_E ) = 0.0;
    J( MOMZ_E , EX     ) = 0.0;
    J( MOMZ_E , EZ     ) = re/lr*rho_e;
    J( MOMZ_E , EZ     ) = 0.0;

    J( EX     , MOMX_I ) = -1.0/(ld*ld*lr)*ri;
    J( EX     , MOMY_I ) = 0.0;
    J( EX     , MOMZ_I ) = 0.0;
    J( EX     , MOMX_E ) = -1.0/(ld*ld*lr)*re;
    J( EX     , MOMY_E ) = 0.0;
    J( EX     , MOMZ_E ) = 0.0;
    J( EX     , EX     ) = 0.0;
    J( EX     , EY     ) = 0.0;
    J( EX     , EZ     ) = 0.0;

    J( EY     , MOMX_I ) = 0.0;
    J( EY     , MOMY_I ) = -1.0/(ld*ld*lr)*ri;
    J( EY     , MOMZ_I ) = 0.0;
    J( EY     , MOMX_E ) = 0.0;
    J( EY     , MOMY_E ) = -1.0/(ld*ld*lr)*re;
    J( EY     , MOMZ_E ) = 0.0;
    J( EY     , EX     ) = 0.0;
    J( EY     , EY     ) = 0.0;
    J( EY     , EZ     ) = 0.0;
 
    J( EZ     , MOMX_I ) = 0.0;
    J( EZ     , MOMY_I ) = 0.0;   
    J( EZ     , MOMZ_I ) = -1.0/(ld*ld*lr)*ri;
    J( EZ     , MOMX_E ) = 0.0;
    J( EZ     , MOMY_E ) = 0.0;   
    J( EZ     , MOMZ_E ) = -1.0/(ld*ld*lr)*re;
    J( EZ     , EX     ) = 0.0;
    J( EZ     , EY     ) = 0.0;
    J( EZ     , EZ     ) = 0.0;
    

    dfdt[MOMX_I] = 0.0;
    dfdt[MOMY_I] = 0.0;
    dfdt[MOMZ_I] = 0.0;
    dfdt[MOMX_E] = 0.0;
    dfdt[MOMY_E] = 0.0;
    dfdt[MOMZ_E] = 0.0;
    dfdt[  EX  ] = 0.0;
    dfdt[  EY  ] = 0.0;
    dfdt[  EZ  ] = 0.0;
  }
};

size_t stiffSolverFull(std::vector<double>& u, double chargeToMassRatioI, double chargeToMassRatioE,
		   double larmorRadius, double debyeLength, double dt)
{
  ri = chargeToMassRatioI;
  re = chargeToMassRatioE;
  lr = larmorRadius;
  ld = debyeLength;
  rho_i = u[0], E_i = u[4], rho_e = u[5], E_e = u[9], Bx = u[10], By = u[11], Bz = u[12];

  vector_type x( 9 , 0.0 );
  x[MOMX_I]=u[1], x[MOMY_I]=u[2], x[MOMZ_I]=u[3], x[MOMX_E]=u[6], x[MOMY_E]=u[7], x[MOMZ_E]=u[8],
    x[EX] = u[13], x[EY] = u[14], x[EZ] = u[15];
  
  size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ) ,
					 make_pair( five_moment_source_full() , five_moment_source_full_jacobi() ) ,
					 x , 0.0 , dt , dt/10.0);
  typedef rosenbrock4< double > stepper_type;
  typedef rosenbrock4_controller< stepper_type > controlled_stepper_type;
  /*
  size_t num_of_steps = integrate_const( controlled_stepper_type() ,
					 make_pair( five_moment_source() , five_moment_source_jacobi() ) ,
					 x , 0.0 , dt , dt/100.0);*/
  u[1]=x[MOMX_I], u[2]=x[MOMY_I], u[3]=x[MOMZ_I], u[6]=x[MOMX_E], u[7]=x[MOMY_E], u[8]=x[MOMZ_E],
    u[13]=x[EX], u[14]=x[EY], u[15]=x[EZ];

  return num_of_steps;
}
