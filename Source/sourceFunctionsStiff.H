#include <iostream>
#include <fstream>
#include <utility>

/*
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint/stepper/generation/make_dense_output.hpp>
#include <boost/numeric/odeint/stepper/rosenbrock4.hpp>
*/
#include <boost/numeric/odeint.hpp>

#include <boost/phoenix/core.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

//using namespace std;
//using namespace boost::numeric::odeint;
//namespace phoenix = boost::phoenix;


//[ stiff_system_definition
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;


size_t stiffSolver(std::vector<double>& u, double chargeToMassRatio, double larmorRadius, double dt);
size_t stiffSolverFull(std::vector<double>& u, double chargeToMassRatioI, double chargeToMassRatioE,
		       double larmorRadius, double debyeLength, double dt);
