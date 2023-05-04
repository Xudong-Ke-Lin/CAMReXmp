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
