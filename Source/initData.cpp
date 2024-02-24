#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

// needed for EMwaveTM
//#include <complex>

using namespace amrex;

//
//Initialize grid data at problem start-up.
//
void
CAMReXmp::initData ()
{
  //
  // Loop over grids, call FORTRAN function to init with data.
  //
  const Real* dx  = geom.CellSize();
  // Position of the bottom left corner of the domain
  const Real* prob_lo = geom.ProbLo();
  // Create a multifab which can store the initial data
  MultiFab& S_new = get_new_data(Phi_Type);
  Real cur_time   = state[Phi_Type].curTime();

#if (AMREX_SPACEDIM >= 2)  
  // Set up a multifab that will contain the electromagnetic fields
  MultiFab& S_EM_X = get_new_data(EM_X_Type);
  MultiFab& S_EM_Y = get_new_data(EM_Y_Type);
  MultiFab& S_EM_XY = get_new_data(EM_XY_Type);
  
  BoxArray ba = S_new.boxArray();
  const DistributionMapping& dm = S_new.DistributionMap();
  
  S_EM_X.define(convert(ba,IntVect{AMREX_D_DECL(1,0,0)}), dm, 6, 2);
  S_EM_Y.define(convert(ba,IntVect{AMREX_D_DECL(0,1,0)}), dm, 6, 2);
  S_EM_XY.define(convert(ba,IntVect{AMREX_D_DECL(1,1,0)}), dm, 6, 2);
#endif
  
  // amrex::Print works like std::cout, but in parallel only prints from the root processor
  if (verbose) {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

  // Slightly messy way to ensure uninitialised data is not used.
  // AMReX has an XDim3 object, but a function needs to be written to
  // convert Real* to XDim3
  const Real dX = dx[0];
  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);
  const Real dZ = (amrex::SpaceDim > 2 ? dx[2] : 0.0);

  const Real probLoX = prob_lo[0];
  const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0);
  const Real probLoZ = (amrex::SpaceDim > 2 ? prob_lo[2] : 0.0);

  ParmParse pp;
  pp.query("test",test);
  pp.query("Gamma", Gamma);
  pp.query("r_i",r_i);
  //pp.query("r_e",r_e);
  pp.query("m",m);
  pp.query("l_r",l_r);
  pp.query("c",c);
  pp.query("lambda_d",lambda_d);
  r_e = -r_i*m;

  implicitMaxwellSolverSetUp();
  /*
  if (MaxwellMethod=="IM"){  
    // set up boundary conditions and coefficients for the implicit solver
    implicitMaxwellSolverSetUp();
  }
  */
  Real rho, v_x, v_y, v_z, p, B_x, B_y, B_z;

#if (AMREX_SPACEDIM >= 2)     
  // Set values for the x-components of the EM fields at the x-faces
  for (MFIter mfi(S_EM_X); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_EM_X.array(mfi);
    
    for(int k = lo.z; k <= hi.z; k++)
    {
      const Real z = probLoZ + (double(k)+0.5) * dZ;
      for(int j = lo.y; j <= hi.y; j++)
      {
	const Real y = probLoY + (double(j)+0.5) * dY;
	for(int i = lo.x; i <= hi.x; i++)
	{
	  // only x-face has no 0.5 shift
	  const Real x = probLoX + (double(i)) * dX;
	
	  if (test=="BrioWu"){
	    Real B_x_L = 0.75, B_y_L = 1.0, B_z_L = 0.0;
	    Real B_x_R = 0.75, B_y_R = -1.0, B_z_R = 0.0;
	    
	    if (x<=(geom.ProbLo()[0]+geom.ProbHi()[0])/2.0){
	      B_x = B_x_L, B_y = B_y_L, B_z = B_z_L;
	    } else{
	      B_x = B_x_R, B_y = B_y_R, B_z = B_z_R;
	    }
	    /*if (y<=(geom.ProbLo()[1]+geom.ProbHi()[1])/2.0){
	      B_x = B_y_L, B_y = B_x_L, B_z = B_z_L;
	    } else{
	      B_x = B_y_R, B_y = B_x_R, B_z = B_z_R;
	      }*/

	    arr(i,j,k,BX_LOCAL) = B_x;	    
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;	    
	    
	  } else if (test=="BrioWu2DCart"){
	    Real theta = std::atan(std::abs(y/x));
	    // second quadrant
	    if (x<=0.0 && y>=0.0)
	      theta = M_PI-theta;
	    // third quadrant
	    else if (x<=0.0 && y<=0.0)
	      theta = M_PI+theta;
	    // fourth quadrant
	    else if (x>=0.0 && y<=0.0)
	      theta = 2*M_PI-theta;

	    Real rCyl = std::sqrt(x*x+y*y);	    
	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real omega = xi0/(0.1*2.0*M_PI);

	    // function
	    auto B_theta = std::cyl_bessel_j(1,omega*rCyl);	    
	    
	    B_x = -B_theta*std::sin(theta);
	    B_y = B_theta*std::cos(theta);
	    B_z = 0.0;

	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;	    
	    arr(i,j,k,EZ_LOCAL) = 0.0;	   

	  } else if (test=="BrioWu1DCyl"){

	    // cylindrical coordinates
	    const Real rCyl = x;

	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real omega = xi0/(0.1*2.0*M_PI);

	    // function
	    auto B_theta = std::cyl_bessel_j(1,omega*rCyl);	    

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = B_theta;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;	    
	    arr(i,j,k,EZ_LOCAL) = 0.0;	   

	  } else if (test=="OT" || test=="OTideal"){
	    v_x = -std::sin(y);
	    v_y = std::sin(x);
	    v_z = 0.0;	    
	    B_x = -std::sin(y);
	    B_y = std::sin(2.0*x);
	    B_z = 0.0;
	    Vector<Real> vcrossB = cross_product({v_x,v_y,v_z},{B_x,B_y,B_z});
	    
	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = -vcrossB[0];
	    arr(i,j,k,EY_LOCAL) = -vcrossB[1];
	    arr(i,j,k,EZ_LOCAL) = -vcrossB[2];

	  } else if (test=="Harris_sheet"){
	    Real Lx = geom.ProbHi()[0]-geom.ProbLo()[0], Ly = geom.ProbHi()[1]-geom.ProbLo()[1];
	    Real lambda = 0.5, B0 = 1.0, B1 = 0.1;
	    B_x = B0*std::tanh(y/lambda) - B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*y/Ly);
	    B_y = B1*(2*M_PI/Lx)*std::cos(M_PI*y/Ly)*std::sin(2*M_PI*x/Lx);

	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="Harris_sheet_full"){
	    Real Lx = geom.ProbHi()[0]-geom.ProbLo()[0], Ly = geom.ProbHi()[1]-geom.ProbLo()[1];
	    Real lambda = 0.5, B0 = 1.0, B1 = 0.1;
	    if (y>=0.0)
	      B_x = B0*std::tanh((y-Ly/4)/lambda) - 2.0*B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*(y-0.25*Ly)/(Ly/2.0));
	    else
	      B_x = - B0*std::tanh((y+Ly/4)/lambda) - 2.0*B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*(y-0.25*Ly)/(Ly/2.0));
	    B_y = B1*(2*M_PI/Lx)*std::cos(M_PI*(y-0.25*Ly)/(Ly/2.0))*std::sin(2*M_PI*x/Lx);

	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="blast"){
	    Real B0 = 100.0/std::sqrt(4.0*M_PI);
	    
	    arr(i,j,k,BX_LOCAL) = B0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="zpinch1d"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    
	    Real Rp = 1.0/8.0, J0 = 1.0/10.0;
	      
	    Real B_theta_in = J0*(0.5*rCyl-16.0*rCyl*rCyl*rCyl), B_theta_out=J0*(0.5*Rp-16.0*Rp*Rp*Rp)*Rp/rCyl;
	    if (rCyl<Rp){
	      B_z = -B_theta_in;
	    } else{
	      B_z = -B_theta_out;
	    }

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    // minus sign since in cylindrical coordinates
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="zpinch2d"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    Real zCyl = y;
	    
	    Real Rp = 1.0/4.0, J0 = 1.0, epsilon = 1.0/100.0, K = 1.0;
	    /*Real epsilon0 = lambda_d*lambda_d;
	      Real mu0 = 1.0/(c*c*epsilon0);*/
	    Real mu0 = 1.0;

	    Real B_phi_in = 0.5*rCyl*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl)), B_phi_out = 0.5*(Rp*Rp/rCyl)*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl)); 
	    
	    if (rCyl<Rp){
	      // minus sign since in cylindrical coordinates
	      B_z = -B_phi_in;
	    } else{
	      B_z = -B_phi_out;
	    }

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="zpinch2dTrue"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    Real zCyl = y;

	    // parameters
	    Real epsilon0 = 1.0, mu0 = 1.0, n0 = 1.0, P0 = 1.5625*0.01, beta = 1.0/10.0, alpha = 100.0;
	    Real m_i = 1.0;
	    Real q_i = r_i*m_i;
	    // perturbation
	    Real epsilon = 1.0/100.0, K = 1.0;

	    Real oneR2A = 1.0 + rCyl*rCyl*alpha;
	    Real oneBR2AB = 1.0 + beta + rCyl*rCyl*alpha*beta;
	    
	    Real C1 = -2.0 + rCyl*rCyl*alpha*(-2.0+(P0*alpha*epsilon0)/(n0*n0*q_i*q_i*std::pow(oneBR2AB,2)));
	    Real C2 = std::sqrt(P0*((rCyl*rCyl*alpha*C1)/std::pow(oneR2A,2)+2.0*std::log(oneR2A))/(rCyl*rCyl*alpha*mu0));

	    Real Btheta = -C2*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl));
	    Real Er = -(P0*rCyl*alpha)/(q_i*std::pow(oneR2A,2)*(n0/oneR2A+n0*beta));
	    //Er *= (1.0 + epsilon*std::sin(2*M_PI*K*zCyl));
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    // minus sign since in cylindrical coordinates
	    arr(i,j,k,BZ_LOCAL) = -Btheta;
	    arr(i,j,k,EX_LOCAL) = Er;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="convergence"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = -std::sin(2*M_PI*x);

	  } else if (test=="convergence2D"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = -c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="EMwave"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;	   
	    arr(i,j,k,EX_LOCAL) = -c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EY_LOCAL) = c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="EMwave1d"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;	   
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = c*std::cos(2.0*M_PI*(x));
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="EMwaveTM"){

	    using namespace std::complex_literals;
	    
	    // cylindrical coordinates
	    const Real rCyl = x;
	    const Real zCyl = y;

	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real kr = xi0;
	    Real kz = M_PI;
	    Real omega = std::sqrt(kr*kr + kz*kz);

	    // function
	    auto Bessel0ExpFunc = std::cyl_bessel_j(0,kr*rCyl)*std::exp(1i*kz*zCyl);
	    auto Bessel1ExpFunc = std::cyl_bessel_j(1,kr*rCyl)*std::exp(1i*kz*zCyl);
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    // negative sign for theta components in cyl. coord.
	    arr(i,j,k,BZ_LOCAL) = -std::real(-1i*omega*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
	    arr(i,j,k,EX_LOCAL) = std::real(-1i*kz*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
	    arr(i,j,k,EY_LOCAL) = std::real(Bessel0ExpFunc);
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="EMwaveTE1d"){

	    using namespace std::complex_literals;
	    
	    // cylindrical coordinates
	    const Real rCyl = x;
	    Real L = 2.0*M_PI; //geom.ProbHi()[0]-geom.ProbLo()[0];
	    // parameters
	    Real z1 = 2.40482555769577;
	    Real omega = z1/L;
	    
	    // function
	    auto Bessel0SinFunc = std::cyl_bessel_j(0,omega*rCyl)*std::sin(omega*0.0);
	    auto Bessel1CosFunc = std::cyl_bessel_j(1,omega*rCyl)*std::cos(omega*0.0);
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = -Bessel0SinFunc;	    
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    // negative sign for theta components in cyl. coord.
	    arr(i,j,k,EZ_LOCAL) = -Bessel1CosFunc;
	    
	  } else if (test=="gaussianEM"){

	    Real lambda = 1.5, chi = 1.5, a = -2.5, b = -2.5;
	    Real COS = std::cos(2.0*M_PI*(x+y)/lambda);
	    Real SIN = std::sin(2.0*M_PI*(x+y)/lambda);
	    Real FACTOR = ((x-a)*(x-a)+(y-b)*(y-b))/(chi*chi);
	    Real EXP = std::exp(-FACTOR);
	    Real epsilon = 5.0 - 4.0*std::tanh((std::sqrt(x*x+y*y)-0.75)/0.08);

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = -COS*EXP + lambda*SIN*EXP*(y-b)/(chi*chi*M_PI);
	    arr(i,j,k,EX_LOCAL) *= c/(epsilon*std::sqrt(2.0));
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	  } 
	}
      }
    }
  }
  // Set values for the y-components of the EM fields at the y-faces
  for (MFIter mfi(S_EM_Y); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_EM_Y.array(mfi);
    
    for(int k = lo.z; k <= hi.z; k++)
    {
      const Real z = probLoZ + (double(k)+0.5) * dZ;
      for(int j = lo.y; j <= hi.y; j++)
      {
  	// only y-face has no 0.5 shift
  	const Real y = probLoY + (double(j)) * dY;
  	for(int i = lo.x; i <= hi.x; i++)
  	{
  	  const Real x = probLoX + (double(i)+0.5) * dX;

  	  if (test=="BrioWu"){
  	    Real B_x_L = 0.75, B_y_L = 1.0, B_z_L = 0.0;
  	    Real B_x_R = 0.75, B_y_R = -1.0, B_z_R = 0.0;
	    
  	    if (x<=(geom.ProbLo()[0]+geom.ProbHi()[0])/2.0){
  	      B_x = B_x_L, B_y = B_y_L, B_z = B_z_L;
  	    } else{
  	      B_x = B_x_R, B_y = B_y_R, B_z = B_z_R;
	    }
  	    /*if (y<=(geom.ProbLo()[1]+geom.ProbHi()[1])/2.0){
  	      B_x = B_y_L, B_y = B_x_L, B_z = B_z_L;
  	    } else{
  	      B_x = B_y_R, B_y = B_x_R, B_z = B_z_R;
	      }*/
	    
  	    arr(i,j,k,BX_LOCAL) = B_x;
  	    arr(i,j,k,BY_LOCAL) = B_y;
  	    arr(i,j,k,BZ_LOCAL) = 0.0;
  	    arr(i,j,k,EX_LOCAL) = 0.0;
  	    arr(i,j,k,EY_LOCAL) = 0.0;
  	    arr(i,j,k,EZ_LOCAL) = 0.0;	    
	    
	  } else if (test=="BrioWu2DCart"){	
	    Real theta = std::atan(std::abs(y/x));
	    // second quadrant
	    if (x<=0.0 && y>=0.0)
	      theta = M_PI-theta;
	    // third quadrant
	    else if (x<=0.0 && y<=0.0)
	      theta = M_PI+theta;
	    // fourth quadrant
	    else if (x>=0.0 && y<=0.0)
	      theta = 2*M_PI-theta;

	    Real rCyl = std::sqrt(x*x+y*y);
	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real omega = xi0/(0.1*2.0*M_PI);

	    // function
	    auto B_theta = std::cyl_bessel_j(1,omega*rCyl);	    
	    
	    B_x = -B_theta*std::sin(theta);
	    B_y = B_theta*std::cos(theta);
	    B_z = 0.0;

	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;	    
	    arr(i,j,k,EZ_LOCAL) = 0.0;	   

	  } else if (test=="BrioWu1DCyl"){

	    // cylindrical coordinates
	    const Real rCyl = x;

	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real omega = xi0/(0.1*2.0*M_PI);

	    // function
	    auto B_theta = std::cyl_bessel_j(1,omega*rCyl);	    

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = B_theta;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;	    
	    arr(i,j,k,EZ_LOCAL) = 0.0;	   

  	  } else if (test=="OT" || test=="OTideal"){
  	    v_x = -std::sin(y);
  	    v_y = std::sin(x);
  	    v_z = 0.0;	    
  	    B_x = -std::sin(y);
  	    B_y = std::sin(2.0*x);
  	    B_z = 0.0;
  	    Vector<Real> vcrossB = cross_product({v_x,v_y,v_z},{B_x,B_y,B_z});	    

	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = -vcrossB[0];
	    arr(i,j,k,EY_LOCAL) = -vcrossB[1];
	    arr(i,j,k,EZ_LOCAL) = -vcrossB[2];
	    	    
  	  } else if (test=="Harris_sheet"){
  	    Real Lx = geom.ProbHi()[0]-geom.ProbLo()[0], Ly = geom.ProbHi()[1]-geom.ProbLo()[1];
  	    Real lambda = 0.5, B0 = 1.0, B1 = 0.1;
  	    B_x = B0*std::tanh(y/lambda) - B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*y/Ly);
	    B_y = B1*(2*M_PI/Lx)*std::cos(M_PI*y/Ly)*std::sin(2*M_PI*x/Lx);	    

  	    arr(i,j,k,BX_LOCAL) = B_x;
  	    arr(i,j,k,BY_LOCAL) = B_y;
  	    arr(i,j,k,BZ_LOCAL) = 0.0;
  	    arr(i,j,k,EX_LOCAL) = 0.0;
  	    arr(i,j,k,EY_LOCAL) = 0.0;
  	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
  	  } else if (test=="Harris_sheet_full"){
  	    Real Lx = geom.ProbHi()[0]-geom.ProbLo()[0], Ly = geom.ProbHi()[1]-geom.ProbLo()[1];
  	    Real lambda = 0.5, B0 = 1.0, B1 = 0.1;
	    if (y>=0.0)
	      B_x = B0*std::tanh((y-Ly/4)/lambda) - 2.0*B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*(y-0.25*Ly)/(Ly/2.0));
	    else
	      B_x = - B0*std::tanh((y+Ly/4)/lambda) - 2.0*B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*(y-0.25*Ly)/(Ly/2.0));
	    B_y = B1*(2*M_PI/Lx)*std::cos(M_PI*(y-0.25*Ly)/(Ly/2.0))*std::sin(2*M_PI*x/Lx);

  	    arr(i,j,k,BX_LOCAL) = B_x;
  	    arr(i,j,k,BY_LOCAL) = B_y;
  	    arr(i,j,k,BZ_LOCAL) = 0.0;
  	    arr(i,j,k,EX_LOCAL) = 0.0;
  	    arr(i,j,k,EY_LOCAL) = 0.0;
  	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
  	  } else if (test=="blast"){
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="zpinch1d"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    
	    Real Rp = 1.0/8.0, J0 = 1.0/10.0;
	      
	    Real B_theta_in = J0*(0.5*rCyl-16.0*rCyl*rCyl*rCyl), B_theta_out=J0*(0.5*Rp-16.0*Rp*Rp*Rp)*Rp/rCyl;
	    if (rCyl<Rp){
	      B_z = -B_theta_in;
	    } else{
	      B_z = -B_theta_out;
	    }

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    // minus sign since in cylindrical coordinates
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="zpinch2d"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    Real zCyl = y;
	    
	    Real Rp = 1.0/4.0, J0 = 1.0, epsilon = 1.0/100.0, K = 1.0;
	    /*Real epsilon0 = lambda_d*lambda_d;
	      Real mu0 = 1.0/(c*c*epsilon0);*/
	    Real mu0 = 1.0;

	    Real B_phi_in = 0.5*rCyl*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl)), B_phi_out = 0.5*(Rp*Rp/rCyl)*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl)); 
	    
	    if (rCyl<Rp){
	      // minus sign since in cylindrical coordinates
	      B_z = -B_phi_in;
	    } else{
	      B_z = -B_phi_out;
	    }

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="zpinch2dTrue"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    Real zCyl = y;

	    // parameters
	    Real epsilon0 = 1.0, mu0 = 1.0, n0 = 1.0, P0 = 1.5625*0.01, beta = 1.0/10.0, alpha = 100.0;
	    Real m_i = 1.0;
	    Real q_i = r_i*m_i;
	    // perturbation
	    Real epsilon = 1.0/100.0, K = 1.0;

	    Real oneR2A = 1.0 + rCyl*rCyl*alpha;
	    Real oneBR2AB = 1.0 + beta + rCyl*rCyl*alpha*beta;
	    
	    Real C1 = -2.0 + rCyl*rCyl*alpha*(-2.0+(P0*alpha*epsilon0)/(n0*n0*q_i*q_i*std::pow(oneBR2AB,2)));
	    Real C2 = std::sqrt(P0*((rCyl*rCyl*alpha*C1)/std::pow(oneR2A,2)+2.0*std::log(oneR2A))/(rCyl*rCyl*alpha*mu0));

	    Real Btheta = -C2*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl));
	    Real Er = -(P0*rCyl*alpha)/(q_i*std::pow(oneR2A,2)*(n0/oneR2A+n0*beta));
	    //Er *= (1.0 + epsilon*std::sin(2*M_PI*K*zCyl));
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    // minus sign since in cylindrical coordinates
	    arr(i,j,k,BZ_LOCAL) = -Btheta;
	    arr(i,j,k,EX_LOCAL) = Er;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="convergence"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = std::sin(2*M_PI*x);
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = -std::sin(2*M_PI*x);

	  } else if (test=="convergence2D"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="EMwave"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = -c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EY_LOCAL) = c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="EMwave1d"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = c*std::cos(2.0*M_PI*(x));
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="EMwaveTM"){

	    using namespace std::complex_literals;
	    
	    // cylindrical coordinates
	    const Real rCyl = x;
	    const Real zCyl = y;

	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real kr = xi0;
	    Real kz = M_PI;
	    Real omega = std::sqrt(kr*kr + kz*kz);

	    // function
	    auto Bessel0ExpFunc = std::cyl_bessel_j(0,kr*rCyl)*std::exp(1i*kz*zCyl);
	    auto Bessel1ExpFunc = std::cyl_bessel_j(1,kr*rCyl)*std::exp(1i*kz*zCyl);
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    // negative sign for theta components in cyl. coord.
	    arr(i,j,k,BZ_LOCAL) = -std::real(-1i*omega*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
	    arr(i,j,k,EX_LOCAL) = std::real(-1i*kz*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
	    arr(i,j,k,EY_LOCAL) = std::real(Bessel0ExpFunc);
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="EMwaveTE1d"){

	    using namespace std::complex_literals;
	    
	    // cylindrical coordinates
	    const Real rCyl = x;
	    Real L = 2.0*M_PI; //geom.ProbHi()[0]-geom.ProbLo()[0];
	    // parameters
	    Real z1 = 2.40482555769577;
	    Real omega = z1/L;
	    
	    // function
	    auto Bessel0SinFunc = std::cyl_bessel_j(0,omega*rCyl)*std::sin(omega*0.0);
	    auto Bessel1CosFunc = std::cyl_bessel_j(1,omega*rCyl)*std::cos(omega*0.0);
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = -Bessel0SinFunc;	    
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    // negative sign for theta components in cyl. coord.
	    arr(i,j,k,EZ_LOCAL) = -Bessel1CosFunc;
	    
	  } else if (test=="gaussianEM"){

	    Real lambda = 1.5, chi = 1.5, a = -2.5, b = -2.5;
	    Real COS = std::cos(2.0*M_PI*(x+y)/lambda);
	    Real SIN = std::sin(2.0*M_PI*(x+y)/lambda);
	    Real FACTOR = ((x-a)*(x-a)+(y-b)*(y-b))/(chi*chi);
	    Real EXP = std::exp(-FACTOR);
	    Real epsilon = 5.0 - 4.0*std::tanh((std::sqrt(x*x+y*y)-0.75)/0.08);

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = COS*EXP - lambda*SIN*EXP*(x-a)/(chi*chi*M_PI);
	    arr(i,j,k,EY_LOCAL) *= c/(epsilon*std::sqrt(2.0));
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } 
  	}
      }
    }
  }  
  for (MFIter mfi(S_EM_XY); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_EM_XY.array(mfi);
    
    for(int k = lo.z; k <= hi.z; k++)
    {
      const Real z = probLoZ + (double(k)+0.5) * dZ;
      for(int j = lo.y; j <= hi.y; j++)
      {
	const Real y = probLoY + (double(j)) * dY;
	for(int i = lo.x; i <= hi.x; i++)
	{
	  // only x-face has no 0.5 shift
	  const Real x = probLoX + (double(i)) * dX;
	
	  if (test=="BrioWu"){
	    Real B_x_L = 0.75, B_y_L = 1.0, B_z_L = 0.0;
	    Real B_x_R = 0.75, B_y_R = -1.0, B_z_R = 0.0;
	    
	    if (x<=(geom.ProbLo()[0]+geom.ProbHi()[0])/2.0){
	      B_x = B_x_L, B_y = B_y_L, B_z = B_z_L;
	    } else{
	      B_x = B_x_R, B_y = B_y_R, B_z = B_z_R;
	    }
	    /*if (y<=(geom.ProbLo()[1]+geom.ProbHi()[1])/2.0){
	      B_x = B_y_L, B_y = B_x_L, B_z = B_z_L;
	    } else{
	      B_x = B_y_R, B_y = B_x_R, B_z = B_z_R;
	      }*/

	    arr(i,j,k,BX_LOCAL) = B_x;	    
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;	    
	    
	  } else if (test=="BrioWu2DCart"){
	    Real theta = std::atan(std::abs(y/x));
	    // second quadrant
	    if (x<=0.0 && y>=0.0)
	      theta = M_PI-theta;
	    // third quadrant
	    else if (x<=0.0 && y<=0.0)
	      theta = M_PI+theta;
	    // fourth quadrant
	    else if (x>=0.0 && y<=0.0)
	      theta = 2*M_PI-theta;

	    Real rCyl = std::sqrt(x*x+y*y);
	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real omega = xi0/(0.1*2.0*M_PI);

	    // function
	    auto B_theta = std::cyl_bessel_j(1,omega*rCyl);	    
	    
	    B_x = -B_theta*std::sin(theta);
	    B_y = B_theta*std::cos(theta);
	    B_z = 0.0;

	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;	    
	    arr(i,j,k,EZ_LOCAL) = 0.0;	   

	  } else if (test=="BrioWu1DCyl"){

	    // cylindrical coordinates
	    const Real rCyl = x;

	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real omega = xi0/(0.1*2.0*M_PI);

	    // function
	    auto B_theta = std::cyl_bessel_j(1,omega*rCyl);	    

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = B_theta;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;	    
	    arr(i,j,k,EZ_LOCAL) = 0.0;	   

	  } else if (test=="OT" || test=="OTideal"){
	    v_x = -std::sin(y);
	    v_y = std::sin(x);
	    v_z = 0.0;	    
	    B_x = -std::sin(y);
	    B_y = std::sin(2.0*x);
	    B_z = 0.0;
	    Vector<Real> vcrossB = cross_product({v_x,v_y,v_z},{B_x,B_y,B_z});	    
	    
	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = -vcrossB[0];
	    arr(i,j,k,EY_LOCAL) = -vcrossB[1];
	    arr(i,j,k,EZ_LOCAL) = -vcrossB[2];

	  } else if (test=="Harris_sheet"){
	    Real Lx = geom.ProbHi()[0]-geom.ProbLo()[0], Ly = geom.ProbHi()[1]-geom.ProbLo()[1];
	    Real lambda = 0.5, B0 = 1.0, B1 = 0.1;
	    B_x = B0*std::tanh(y/lambda) - B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*y/Ly);
	    B_y = B1*(2*M_PI/Lx)*std::cos(M_PI*y/Ly)*std::sin(2*M_PI*x/Lx);
	    B_z = 0.0;
	    
	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="Harris_sheet_full"){
	    Real Lx = geom.ProbHi()[0]-geom.ProbLo()[0], Ly = geom.ProbHi()[1]-geom.ProbLo()[1];
	    Real lambda = 0.5, B0 = 1.0, B1 = 0.1;
	    if (y>=0.0)
	      B_x = B0*std::tanh((y-Ly/4)/lambda) - 2.0*B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*(y-0.25*Ly)/(Ly/2.0));
	    else
	      B_x = - B0*std::tanh((y+Ly/4)/lambda) - 2.0*B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*(y-0.25*Ly)/(Ly/2.0));
	    B_y = B1*(2*M_PI/Lx)*std::cos(M_PI*(y-0.25*Ly)/(Ly/2.0))*std::sin(2*M_PI*x/Lx);
	    B_z = 0.0;
	    
	    arr(i,j,k,BX_LOCAL) = B_x;
	    arr(i,j,k,BY_LOCAL) = B_y;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="blast"){
	    Real B0 = 100.0/std::sqrt(4.0*M_PI);
	    
	    arr(i,j,k,BX_LOCAL) = B0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="zpinch1d"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    
	    Real Rp = 1.0/8.0, J0 = 1.0/10.0;
	      
	    Real B_theta_in = J0*(0.5*rCyl-16.0*rCyl*rCyl*rCyl), B_theta_out=J0*(0.5*Rp-16.0*Rp*Rp*Rp)*Rp/rCyl;
	    if (rCyl<Rp){
	      B_z = -B_theta_in;
	    } else{
	      B_z = -B_theta_out;
	    }

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    // minus sign since in cylindrical coordinates
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="zpinch2d"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    Real zCyl = y;
	    
	    Real Rp = 1.0/4.0, J0 = 1.0, epsilon = 1.0/100.0, K = 1.0;
	    /*Real epsilon0 = lambda_d*lambda_d;
	      Real mu0 = 1.0/(c*c*epsilon0);*/
	    Real mu0 = 1.0;

	    Real B_phi_in = 0.5*rCyl*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl)), B_phi_out = 0.5*(Rp*Rp/rCyl)*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl)); 
	    
	    if (rCyl<Rp){
	      // minus sign since in cylindrical coordinates
	      B_z = -B_phi_in;
	    } else{
	      B_z = -B_phi_out;
	    }

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = B_z;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="zpinch2dTrue"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    Real zCyl = y;

	    // parameters
	    Real epsilon0 = 1.0, mu0 = 1.0, n0 = 1.0, P0 = 1.5625*0.01, beta = 1.0/10.0, alpha = 100.0;
	    Real m_i = 1.0;
	    Real q_i = r_i*m_i;
	    // perturbation
	    Real epsilon = 1.0/100.0, K = 1.0;

	    Real oneR2A = 1.0 + rCyl*rCyl*alpha;
	    Real oneBR2AB = 1.0 + beta + rCyl*rCyl*alpha*beta;
	    
	    Real C1 = -2.0 + rCyl*rCyl*alpha*(-2.0+(P0*alpha*epsilon0)/(n0*n0*q_i*q_i*std::pow(oneBR2AB,2)));
	    Real C2 = std::sqrt(P0*((rCyl*rCyl*alpha*C1)/std::pow(oneR2A,2)+2.0*std::log(oneR2A))/(rCyl*rCyl*alpha*mu0));

	    Real Btheta = -C2*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl));
	    Real Er = -(P0*rCyl*alpha)/(q_i*std::pow(oneR2A,2)*(n0/oneR2A+n0*beta));
	    //Er *= (1.0 + epsilon*std::sin(2*M_PI*K*zCyl));
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    // minus sign since in cylindrical coordinates
	    arr(i,j,k,BZ_LOCAL) = -Btheta;
	    arr(i,j,k,EX_LOCAL) = Er;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="convergence"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = -std::sin(2*M_PI*x);

	  } else if (test=="convergence2D"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = -c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="EMwave"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = std::cos(2.0*M_PI*(x+y));
	    arr(i,j,k,EX_LOCAL) = -c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EY_LOCAL) = c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="EMwave1d"){

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = std::cos(2.0*M_PI*(x));
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = c*std::cos(2.0*M_PI*(x));
	    arr(i,j,k,EZ_LOCAL) = 0.0;

	  } else if (test=="EMwaveTM"){

	    using namespace std::complex_literals;
	    
	    // cylindrical coordinates
	    const Real rCyl = x;
	    const Real zCyl = y;

	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real kr = xi0;
	    Real kz = M_PI;
	    Real omega = std::sqrt(kr*kr + kz*kz);

	    // function
	    auto Bessel0ExpFunc = std::cyl_bessel_j(0,kr*rCyl)*std::exp(1i*kz*zCyl);
	    auto Bessel1ExpFunc = std::cyl_bessel_j(1,kr*rCyl)*std::exp(1i*kz*zCyl);
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    // negative sign for theta components in cyl. coord.
	    arr(i,j,k,BZ_LOCAL) = -std::real(-1i*omega*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
	    arr(i,j,k,EX_LOCAL) = std::real(-1i*kz*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
	    arr(i,j,k,EY_LOCAL) = std::real(Bessel0ExpFunc);
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	    
	  } else if (test=="EMwaveTE1d"){

	    using namespace std::complex_literals;
	    
	    // cylindrical coordinates
	    const Real rCyl = x;
	    Real L = 2.0*M_PI; //geom.ProbHi()[0]-geom.ProbLo()[0];
	    // parameters
	    Real z1 = 2.40482555769577;
	    Real omega = z1/L;
	    
	    // function
	    auto Bessel0SinFunc = std::cyl_bessel_j(0,omega*rCyl)*std::sin(omega*0.0);
	    auto Bessel1CosFunc = std::cyl_bessel_j(1,omega*rCyl)*std::cos(omega*0.0);
	    
	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = -Bessel0SinFunc;	    
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = 0.0;
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    // negative sign for theta components in cyl. coord.
	    arr(i,j,k,EZ_LOCAL) = -Bessel1CosFunc;
	    
	  } else if (test=="gaussianEM"){

	    Real lambda = 1.5, chi = 1.5, a = -2.5, b = -2.5;
	    Real COS = std::cos(2.0*M_PI*(x+y)/lambda);
	    Real SIN = std::sin(2.0*M_PI*(x+y)/lambda);
	    Real FACTOR = ((x-a)*(x-a)+(y-b)*(y-b))/(chi*chi);
	    Real EXP = std::exp(-FACTOR);
	    Real epsilon = 5.0 - 4.0*std::tanh((std::sqrt(x*x+y*y)-0.75)/0.08);

	    arr(i,j,k,BX_LOCAL) = 0.0;
	    arr(i,j,k,BY_LOCAL) = 0.0;
	    arr(i,j,k,BZ_LOCAL) = 0.0;
	    arr(i,j,k,EX_LOCAL) = -COS*EXP + lambda*SIN*EXP*(y-b)/(chi*chi*M_PI);
	    arr(i,j,k,EX_LOCAL) *= c/(epsilon*std::sqrt(2.0));
	    arr(i,j,k,EY_LOCAL) = 0.0;
	    arr(i,j,k,EZ_LOCAL) = 0.0;
	  } 
	}
      }
    }
  }
#endif  
  // Set values for cell-centred fluid and EM fields
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_new.array(mfi);
    
    for(int k = lo.z; k <= hi.z; k++)
    {
      const Real z = probLoZ + (double(k)+0.5) * dZ;
      for(int j = lo.y; j <= hi.y; j++)
      {
	const Real y = probLoY + (double(j)+0.5) * dY;
	for(int i = lo.x; i <= hi.x; i++)
	{
	  const Real x = probLoX + (double(i)+0.5) * dX;

	  if (test=="BrioWu"){
	    Real rho_L = 1.0, v_x_L = 0.0, v_y_L = 0.0, v_z_L = 0.0, p_L = 0.5;
	    Real B_x_L = 0.75, B_y_L = 1.0, B_z_L = 0.0;
	    Real rho_R = 0.125, v_x_R = 0.0, v_y_R = 0.0, v_z_R = 0.0, p_R = 0.05;
	    Real B_x_R = 0.75, B_y_R = -1.0, B_z_R = 0.0;
	    
	    if (x<=(geom.ProbLo()[0]+geom.ProbHi()[0])/2.0){
	      rho = rho_L, v_x = v_x_L, v_y = v_y_L, v_z = v_z_L, p = p_L, B_x = B_x_L, B_y = B_y_L, B_z = B_z_L;
	    } else{
	      rho = rho_R, v_x = v_x_R, v_y = v_y_R, v_z = v_z_R, p = p_R, B_x = B_x_R, B_y = B_y_R, B_z = B_z_R;
	    }
	    /*
	    if (y<=(geom.ProbLo()[1]+geom.ProbHi()[1])/2.0){
              rho = rho_L, v_x = v_y_L, v_y = v_x_L, v_z = v_z_L, p = p_L, B_x = B_y_L, B_y = B_x_L, B_z = B_z_L;                           
            } else{                                                                                                       
              rho = rho_R, v_x = v_y_R, v_y = v_x_R, v_z = v_z_R, p = p_R, B_x = B_y_R, B_y = B_x_R, B_z = B_z_R;   
	      }*/	    
	    
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = B_x;
	    arr(i,j,k,BY) = B_y;
	    arr(i,j,k,BZ) = B_z;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;	    
	    arr(i,j,k,EZ) = 0.0;	   
	    arr(i,j,k,DIVB) = 0.0;
	    arr(i,j,k,DIVE) = 0.0;	    
	    
	    // set speed of light
	    //c = 100.0;
	    
	  } else if (test=="BrioWu2DCart"){
	    Real rho_in = 1.0, v_x_in = 0.0, v_y_in = 0.0, v_z_in = 0.0, p_in = 0.5;
	    Real rho_out = 0.125, v_x_out = 0.0, v_y_out = 0.0, v_z_out = 0.0, p_out = 0.05;
	    
	    if (x*x+y*y<=0.4*0.4){
	      rho = rho_in, v_x = v_x_in, v_y = v_y_in, v_z = v_z_in, p = p_in;
	    } else{
	      rho = rho_out, v_x = v_x_out, v_y = v_y_out, v_z = v_z_out, p = p_out;
	    }
	    
	    Real theta = std::atan(std::abs(y/x));
	    // second quadrant
	    if (x<=0.0 && y>=0.0)
	      theta = M_PI-theta;
	    // third quadrant
	    else if (x<=0.0 && y<=0.0)
	      theta = M_PI+theta;
	    // fourth quadrant
	    else if (x>=0.0 && y<=0.0)
	      theta = 2*M_PI-theta;

	    Real rCyl = std::sqrt(x*x+y*y);
	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real omega = xi0/(0.1*2.0*M_PI);

	    // function
	    auto B_theta = std::cyl_bessel_j(1,omega*rCyl);	    
	    
	    B_x = -B_theta*std::sin(theta);
	    B_y = B_theta*std::cos(theta);
	    B_z = 0.0;

	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = B_x;
	    arr(i,j,k,BY) = B_y;
	    arr(i,j,k,BZ) = B_z;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;	    
	    arr(i,j,k,EZ) = 0.0;
	    arr(i,j,k,DIVB) = 0.0;
	    arr(i,j,k,DIVE) = 0.0;	    
	    
	  } else if (test=="BrioWu1DCyl"){
	    Real rho_in = 1.0, v_x_in = 0.0, v_y_in = 0.0, v_z_in = 0.0, p_in = 0.5;
	    Real rho_out = 0.125, v_x_out = 0.0, v_y_out = 0.0, v_z_out = 0.0, p_out = 0.05;
	    
	    if (x<=0.4){
	      rho = rho_in, v_x = v_x_in, v_y = v_y_in, v_z = v_z_in, p = p_in;
	    } else{
	      rho = rho_out, v_x = v_x_out, v_y = v_y_out, v_z = v_z_out, p = p_out;
	    }	   

	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    // cylindrical coordinates
	    const Real rCyl = x;

	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real omega = xi0/(0.1*2.0*M_PI);

	    // function
	    auto B_theta = std::cyl_bessel_j(1,omega*rCyl);
	    
	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = 0.0;
	    arr(i,j,k,BZ) = B_theta;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;	    
	    arr(i,j,k,EZ) = 0.0;	   
	    
	  } else if (test=="OT"){

	    rho = Gamma*Gamma;
	    v_x = -std::sin(y);
	    v_y = std::sin(x);
	    v_z = 0.0;
	    p = Gamma;
	    B_x = -std::sin(y);
	    B_y = std::sin(2.0*x);
	    B_z = 0.0;
	    Vector<Real> vcrossB = cross_product({v_x,v_y,v_z},{B_x,B_y,B_z});

	    arr(i,j,k,0) = rho;	    
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = B_x;
	    arr(i,j,k,BY) = B_y;
	    arr(i,j,k,BZ) = B_z;
	    arr(i,j,k,EX) = -vcrossB[0];
	    arr(i,j,k,EY) = -vcrossB[1];
	    arr(i,j,k,EZ) = -vcrossB[2];
	    arr(i,j,k,DIVB) = 0.0;	    
	    arr(i,j,k,DIVE) = 0.0;

	    // set speed of light
	    //c = std::max(c, 3.0*std::max(std::abs(v_x)+get_speed(w_i),std::abs(v_y)+get_speed(w_i)));
	    
	  } else if (test=="OTideal"){

	    rho = Gamma*Gamma/2.0;
	    v_x = -std::sin(y);
	    v_y = std::sin(x);
	    v_z = 0.0;
	    p = Gamma/2.0;
	    B_x = -std::sin(y);
	    B_y = std::sin(2.0*x);
	    B_z = 0.0;
	    Vector<Real> vcrossB = cross_product({v_x,v_y,v_z},{B_x,B_y,B_z});

	    arr(i,j,k,0) = rho;	    
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = B_x;
	    arr(i,j,k,BY) = B_y;
	    arr(i,j,k,BZ) = B_z;
	    arr(i,j,k,EX) = -vcrossB[0];
	    arr(i,j,k,EY) = -vcrossB[1];
	    arr(i,j,k,EZ) = -vcrossB[2];
	    arr(i,j,k,DIVB) = 0.0;	    
	    arr(i,j,k,DIVE) = 0.0;
	    
	  } else if (test=="Harris_sheet"){
	    Real Lx = geom.ProbHi()[0]-geom.ProbLo()[0], Ly = geom.ProbHi()[1]-geom.ProbLo()[1];
	    Real lambda = 0.5, n0 = 1.0, nInf = 0.2, B0 = 1.0, B1 = 0.1;
	    Real n = n0*(nInf + 1.0/(std::cosh(y/lambda)*std::cosh(y/lambda)));
	    p = B0*B0/2.0 * n/n0;
	    v_x = 0.0, v_y = 0.0, v_z = 0.0;
	    Real v_z_e = 1.0/(std::cosh(y/lambda)*std::cosh(y/lambda)) * 1.0/(n*lambda);
	    B_x = B0*std::tanh(y/lambda) - B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*y/Ly);
	    B_y = B1*(2*M_PI/Lx)*std::cos(M_PI*y/Ly)*std::sin(2*M_PI*x/Lx);
	    B_z = 0.0;	    

	    arr(i,j,k,0) = n;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p*5.0/6.0};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = n/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z_e;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z_e, p/6.0};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = B_x;
	    arr(i,j,k,BY) = B_y;
	    arr(i,j,k,BZ) = B_z;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;
	    arr(i,j,k,EZ) = 0.0;
	    arr(i,j,k,DIVB) = 0.0;	    
	    arr(i,j,k,DIVE) = 0.0;

	  } else if (test=="Harris_sheet_full"){
	    Real Lx = geom.ProbHi()[0]-geom.ProbLo()[0], Ly = geom.ProbHi()[1]-geom.ProbLo()[1];
	    Real lambda = 0.5, n0 = 1.0, nInf = 0.2, B0 = 1.0, B1 = 0.1;
	    Real n = n0*(nInf + 1.0/(std::cosh((y-Ly/4)/lambda)*std::cosh((y-Ly/4)/lambda)) +  1.0/(std::cosh((y+Ly/4)/lambda)*std::cosh((y+Ly/4)/lambda)));
	    p = B0*B0/2.0 * n/n0;
	    v_x = 0.0, v_y = 0.0, v_z = 0.0;
	    Real v_z_e = (1.0/(std::cosh((y-Ly/4)/lambda)*std::cosh((y-Ly/4)/lambda)) +  1.0/(std::cosh((y+Ly/4)/lambda)*std::cosh((y+Ly/4)/lambda))) * 1.0/(n*lambda);
	    if (y<0.0)
	      v_z_e *= -1.0;
	    if (y>=0.0)
	      B_x = B0*std::tanh((y-Ly/4)/lambda) - 2.0*B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*(y-0.25*Ly)/(Ly/2.0));
	    else
	      B_x = - B0*std::tanh((y+Ly/4)/lambda) - 2.0*B1*(M_PI/Ly)*std::cos(2*M_PI*x/Lx)*std::sin(M_PI*(y-0.25*Ly)/(Ly/2.0));
	    B_y = B1*(2*M_PI/Lx)*std::cos(M_PI*(y-0.25*Ly)/(Ly/2.0))*std::sin(2*M_PI*x/Lx);
	    B_z = 0.0;	    

	    arr(i,j,k,0) = n;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p*5.0/6.0};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = n/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z_e;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z_e, p/6.0};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = B_x;
	    arr(i,j,k,BY) = B_y;
	    arr(i,j,k,BZ) = B_z;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;
	    arr(i,j,k,EZ) = 0.0;
	    arr(i,j,k,DIVB) = 0.0;	    
	    arr(i,j,k,DIVE) = 0.0;

	  } else if (test=="blast"){

	    rho = 0.5, v_x = 0.0, v_y = 0.0, v_z = 0.0;
	    Real p_in = 500.0, p_out = 0.05;
	    Real B0 = 100.0/std::sqrt(4.0*M_PI);
	    
	    if (x*x+y*y <= 0.1*0.1)
	      p = p_in;
	    else
	      p = p_out;
	    	
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = B0;
	    arr(i,j,k,BY) = 0.0;
	    arr(i,j,k,BZ) = 0.0;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;
	    arr(i,j,k,EZ) = 0.0;

	  } else if (test=="zpinch1d"){

	    // cylindrical coordinates
	    Real rCyl = x;

	    Real m_e = 1.0/m;
	    Real q_e = r_e*m_e;
	    
	    Real Rp = 1.0/8.0, alpha = 1.0/10.0, J0 = 1.0/10.0;
	    Real p0 = J0*J0/(1-alpha)*(0.25*Rp*Rp-12.0*std::pow(Rp,4)+(4.0/3.0)*128.0*std::pow(Rp,6));
	      
	    Real j_z_e_in = J0*(1.0-64.0*rCyl*rCyl), j_z_e_out = 0.0;
	    Real B_theta_in = J0*(0.5*rCyl-16.0*rCyl*rCyl*rCyl), B_theta_out=J0*(0.5*Rp-16.0*Rp*Rp*Rp)*Rp/rCyl;
	    Real p_in = p0 - J0*J0*(0.25*rCyl*rCyl-12.0*std::pow(rCyl,4)+(4.0/3.0)*128.0*std::pow(rCyl,6));
	    Real p_out = p0 - J0*J0*(0.25*Rp*Rp-12.0*std::pow(Rp,4)+(4.0/3.0)*128.0*std::pow(Rp,6));

	    Real v_r_in = 0.0, v_theta_in = 0.0, v_z_in = 0.0;
	    Real v_r_out = 0.0, v_theta_out = 0.0, v_z_out = 0.0;
	    
	    // number density
	    Real n_in = p_in/p0, n_out = p_out/p0;

	    // electron velocity
	    Real v_z_e;

	    if (rCyl<Rp){
	      rho = n_in, v_z_e = j_z_e_in/(q_e*n_in), p = p_in/2.0, B_z = -B_theta_in;
	      v_x = v_r_in, v_y = v_z_in, v_z = -v_theta_in;
	    } else{
	      rho = n_out, v_z_e = j_z_e_out/(q_e*n_out), p = p_out/2.0, B_z = -B_theta_out;
	      v_x = v_r_out, v_y= v_z_out, v_z = -v_theta_out;
	    }

	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_z_e;	    
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_z_e, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = 0.0;
	    // minus sign since in cylindrical coordinates
	    arr(i,j,k,BZ) = B_z;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;
	    arr(i,j,k,EZ) = 0.0;

	  } else if (test=="zpinch2d"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    Real zCyl = y;

	    Real m_e = 1.0/m;
	    Real q_e = r_e*m_e/l_r;

	    Real Rp = 1.0/4.0, alpha = 1.0/10.0, J0 = 1.0, epsilon = 1.0/100.0, K = 1.0;
	    /*Real epsilon0 = lambda_d*lambda_d;
	      Real mu0 = 1.0/(c*c*epsilon0);*/
	    Real mu0 = 1.0;//*l_r;
	    Real p0 = (1+alpha)*mu0*J0*J0*Rp*Rp/4.0;
	      
	    Real j_z_e_in = J0, j_z_e_out = 0.0;
	    //Real B_phi_in = 0.5*x*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*y)), B_phi_out = 0.5*(Rp*Rp/x)*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*y));
	    Real B_phi_in = 0.5*rCyl*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl)), B_phi_out = 0.5*(Rp*Rp/rCyl)*mu0*J0*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl)); 
	    //Real p_in = p0 - 0.25*mu0*J0*J0*x*x, p_out  = alpha*mu0*J0*J0*Rp*Rp;
	    Real p_in = p0 - 0.25*mu0*J0*J0*rCyl*rCyl, p_out  = 0.25*alpha*mu0*J0*J0*Rp*Rp;

	    Real v_x_L = 0.0, v_y_L = 0.0, v_z_L = 0.0;
	    Real v_x_R = 0.0, v_y_R = 0.0, v_z_R = 0.0;
	    
	    // number density
	    Real n_in = p_in/p0, n_out = p_out/p0;

	    // electron velocity
	    Real v_z_e;
	    
	    if (rCyl<Rp){	    
	      rho = n_in, v_z_e = j_z_e_in/(q_e*n_in), p = p_in/2.0, B_z = -B_phi_in;
	      v_x = v_x_L, v_y = v_y_L, v_z = v_z_L;
	      //v_z_e = -j_z_e_in/m;	      
	    } else{
	      rho = n_out, v_z_e = j_z_e_out/(q_e*n_out), p = p_out/2.0, B_z = -B_phi_out;
	      v_x = v_x_R, v_y = v_y_R, v_z = v_z_R;
	      //v_z_e = -j_z_e_out/m;
	    }
	    //v_z_e /= r_i;
	    /*
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z_e;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z_e, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = B_y;
	    arr(i,j,k,BZ) = 0.0;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;
	    arr(i,j,k,EZ) = 0.0;
	    */
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_z_e;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_z_e, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = 0.0;
	    // minus sign since in cylindrical coordinates
	    arr(i,j,k,BZ) = B_z;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;
	    arr(i,j,k,EZ) = 0.0;
	    
	    //c_a = std::max(c_a, get_speed_a(rho, B_y));
	    //c = std::max(c, 8.0*c_a);

	    arr(i,j,k,DIVB) = std::abs(v_z_e/get_speed({rho,0.0,0.0,0.0,arr(i,j,k,ENER_I)}));

	  } else if (test=="zpinch2dTrue"){

	    // cylindrical coordinates
	    Real rCyl = x;
	    Real zCyl = y;

	    // parameters
	    Real epsilon0 = 1.0, mu0 = 1.0, n0 = 1.0, P0 = 1.5625*0.01, beta = 1.0/10.0, alpha = 100.0;
	    Real m_i = 1.0, m_e = 1.0;
	    Real q_i = r_i*m_i;

	    // perturbation
	    Real epsilon = 1.0/100.0, K = 1.0;

	    Real oneR2A = 1.0 + rCyl*rCyl*alpha;
	    Real oneBR2AB = 1.0 + beta + rCyl*rCyl*alpha*beta;
	    Real oneBR4A2B = 1.0 + beta - std::pow(rCyl,4)*alpha*alpha*beta;
	    
	    Real C1 = -2.0 + rCyl*rCyl*alpha*(-2.0+(P0*alpha*epsilon0)/(n0*n0*q_i*q_i*std::pow(oneBR2AB,2)));
	    Real C2 = std::sqrt(P0*((rCyl*rCyl*alpha*C1)/std::pow(oneR2A,2)+2.0*std::log(oneR2A))/(rCyl*rCyl*alpha*mu0));
	    Real C3 = n0*n0*q_i*q_i*oneR2A*std::pow(oneBR2AB,3);

	    Real n_i = n0/oneR2A + n0*beta;
	    Real n_e = (n0*n0*(1.0/oneR2A+beta) + (2.0*P0*alpha*oneBR4A2B*epsilon0)/(q_i*q_i*oneR2A*oneR2A*oneBR2AB*oneBR2AB))/n0;
	    Real P = P0/oneR2A + P0*beta;
	    Real P_i = P/2.0, P_e = P/2.0;
	    Real Btheta = -C2*(1.0 + epsilon*std::sin(2*M_PI*K*zCyl));

	    Real mey = (2.0*P0*rCyl*alpha*(C3+P0*alpha*oneBR4A2B*epsilon0))/(n0*n0*q_i*q_i*q_i*std::pow(oneR2A,3)*std::pow(oneBR2AB,3)*mu0*mu0*C2);
	    Real Er = -(P0*rCyl*alpha)/(q_i*std::pow(oneR2A,2)*(n0/oneR2A+n0*beta));
	    //Er *= (1.0 + epsilon*std::sin(2*M_PI*K*zCyl));
	    
	    arr(i,j,k,0) = n_i/m_i;
	    arr(i,j,k,1) = arr(i,j,k,0)*0.0;
	    arr(i,j,k,2) = arr(i,j,k,0)*0.0;
	    arr(i,j,k,3) = arr(i,j,k,0)*0.0;
	    Vector<Real> w_i{arr(i,j,k,0), 0.0, 0.0, 0.0, P_i};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = n_e/m_e;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*0.0;
	    arr(i,j,k,MOMY_E) = mey;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*0.0;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), 0.0, mey/arr(i,j,k,RHO_E), 0.0, P_e};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = 0.0;
	    // minus sign since in cylindrical coordinates
	    arr(i,j,k,BZ) = -Btheta;
	    arr(i,j,k,EX) = Er;
	    arr(i,j,k,EY) = 0.0;
	    arr(i,j,k,EZ) = 0.0;

	    arr(i,j,k,DIVB) = (mey/arr(i,j,k,RHO_E))/get_speed({n_i,0.0,0.0,0.0,arr(i,j,k,ENER_I)});

	  } else if (test=="convergence"){
	    arr(i,j,k,0) = 2.0+std::sin(2*M_PI*x);
	    arr(i,j,k,1) = arr(i,j,k,0)*1.0;
	    arr(i,j,k,2) = arr(i,j,k,0)*0.0;
	    arr(i,j,k,3) = arr(i,j,k,0)*0.0;
	    Vector<Real> w_i{arr(i,j,k,0), 1.0, 0.0, 0.0, 1.0};
	    arr(i,j,k,ENER_I) = get_energy(w_i);

	    arr(i,j,k,RHO_E) = 2.0+std::sin(2*M_PI*x);
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*1.0;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*0.0;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*0.0;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), 1.0, 0.0, 0.0, 1.0};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = std::sin(2*M_PI*x);
	    arr(i,j,k,BZ) = 0.0;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;
	    arr(i,j,k,EZ) = -std::sin(2*M_PI*x);

	  } else if (test=="convergence2D"){
	    Real vx = c/std::sqrt(2.0), vy = c/std::sqrt(2.0);
	    
	    arr(i,j,k,0) = 2.0+std::sin(2*M_PI*(x+y));
	    arr(i,j,k,1) = arr(i,j,k,0)*vx;
	    arr(i,j,k,2) = arr(i,j,k,0)*vy;
	    arr(i,j,k,3) = arr(i,j,k,0)*0.0;
	    Vector<Real> w_i{arr(i,j,k,0), vx, vy, 0.0, 1.0};
	    arr(i,j,k,ENER_I) = get_energy(w_i);

	    arr(i,j,k,RHO_E) = 2.0+std::sin(2*M_PI*(x+y));
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*vx;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*vy;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*0.0;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), vx, vy, 0.0, 1.0};
	    arr(i,j,k,ENER_E) = get_energy(w_e);

	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = 0.0;
	    arr(i,j,k,BZ) = std::cos(2.0*M_PI*(x+y));
	    arr(i,j,k,EX) = -c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EY) = c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EZ) = 0.0;

	  } else if (test=="Toro1" || test=="Toro2" || test=="Toro3" || test=="Toro4" || test=="Toro5"){
	    Real rho_L, v_x_L, v_y_L, v_z_L, p_L;
	    Real rho_R, v_x_R, v_y_R, v_z_R, p_R;

	    if (test=="Toro1"){
	      rho_L = 1.0, v_x_L = 0.0, v_y_L = 0.0, v_z_L = 0.0, p_L = 1.0;
	      rho_R = 0.125, v_x_R = 0.0, v_y_R = 0.0, v_z_R = 0.0, p_R = 0.1;
	    } else if (test=="Toro2"){
	      rho_L = 1.0, v_x_L = -2.0, v_y_L = 0.0, v_z_L = 0.0, p_L = 0.4;
	      rho_R = 1.0, v_x_R = 2.0, v_y_R = 0.0, v_z_R = 0.0, p_R = 0.4;
	    } else if (test=="Toro3"){
	      rho_L = 1.0, v_x_L = 0.0, v_y_L = 0.0, v_z_L = 0.0, p_L = 1000.0;
	      rho_R = 1.0, v_x_R = 0.0, v_y_R = 0.0, v_z_R = 0.0, p_R = 0.01;
	    } else if (test=="Toro4"){
	      rho_L = 1.0, v_x_L = 0.0, v_y_L = 0.0, v_z_L = 0.0, p_L = 0.01;
	      rho_R = 1.0, v_x_R = 0.0, v_y_R = 0.0, v_z_R = 0.0, p_R = 100.0;
	    } else if (test=="Toro5"){
	      rho_L = 5.99924, v_x_L = 19.5975, v_y_L = 0.0, v_z_L = 0.0, p_L = 460.894;
	      rho_R = 5.99242, v_x_R = -6.19633, v_y_R = 0.0, v_z_R = 0.0, p_R = 46.0950;
	    }
	    
	    if (x<=(geom.ProbLo()[0]+geom.ProbHi()[0])/2.0){
	      rho = rho_L, v_x = v_x_L, v_y = v_y_L, v_z = v_z_L, p = p_L;
	    } else{
	      rho = rho_R, v_x = v_x_R, v_y = v_y_R, v_z = v_z_R, p = p_R;
	    }
	    /*if (y<=(geom.ProbLo()[1]+geom.ProbHi()[1])/2.0){
	      rho = rho_L, v_x = v_x_L, v_y = v_y_L, v_z = v_z_L, p = p_L;
	    } else{
	      rho = rho_R, v_x = v_x_R, v_y = v_y_R, v_z = v_z_R, p = p_R;
	    }*/
	    /*if (x+y<=geom.ProbLo()[1]+geom.ProbHi()[1]){
	      rho = rho_L, v_x = v_x_L, v_y = v_y_L, v_z = v_z_L, p = p_L;
	    } else{
	      rho = rho_R, v_x = v_x_R, v_y = v_y_R, v_z = v_z_R, p = p_R;
	    }*/
	    
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    // set speed of light
	    c = 100.0;
	    
	  } else if (test=="Colella"){

	    Real rho_L = 1.0, v_x_L = 0.0, v_y_L = 0.0, v_z_L = 0.0, p_L = 1000.0;
	    Real rho_C = 1.0, v_x_C = 0.0, v_y_C = 0.0, v_z_C = 0.0, p_C = 0.01;
	    Real rho_R = 1.0, v_x_R = 0.0, v_y_R = 0.0, v_z_R = 0.0, p_R = 100.0;
	    
	    if (x<=0.1*(geom.ProbLo()[0]+geom.ProbHi()[0])){
	      rho = rho_L, v_x = v_x_L, v_y = v_y_L, v_z = v_z_L, p = p_L;
	    } else if (x>=0.9*(geom.ProbLo()[0]+geom.ProbHi()[0])){
	      rho = rho_R, v_x = v_x_R, v_y = v_y_R, v_z = v_z_R, p = p_R;
	    } else{
	      rho = rho_C, v_x = v_x_C, v_y = v_y_C, v_z = v_z_C, p = p_C;
	    }
	    
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    // set speed of light
	    c = 100.0;
	    
	  } else if (test=="cylExp"){

	    Real rho_L = 1.0, v_x_L = 0.0, v_y_L = 0.0, v_z_L = 0.0, p_L = 1.0;
	    Real rho_R = 0.125, v_x_R = 0.0, v_y_R = 0.0, v_z_R = 0.0, p_R = 0.1;
	    
	    if ((x-1.0)*(x-1.0)+(y-1.0)*(y-1.0)<=0.4*0.4){
	      rho = rho_L, v_x = v_x_L, v_y = v_y_L, v_z = v_z_L, p = p_L;
	    } else{
	      rho = rho_R, v_x = v_x_R, v_y = v_y_R, v_z = v_z_R, p = p_R;
	    }
	    
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);
	    
	  } else if (test=="sphExp"){

	    Real rho_L = 1.0, v_x_L = 0.0, v_y_L = 0.0, v_z_L = 0.0, p_L = 1.0;
	    Real rho_R = 0.125, v_x_R = 0.0, v_y_R = 0.0, v_z_R = 0.0, p_R = 0.1;
	    
	    if (x*x+y*y<=0.4*0.4){
	      rho = rho_L, v_x = v_x_L, v_y = v_y_L, v_z = v_z_L, p = p_L;
	    } else{
	      rho = rho_R, v_x = v_x_R, v_y = v_y_R, v_z = v_z_R, p = p_R;
	    }
	    
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    arr(i,j,k,RHO_E) = rho/m;
	    arr(i,j,k,MOMX_E) = arr(i,j,k,RHO_E)*v_x;
	    arr(i,j,k,MOMY_E) = arr(i,j,k,RHO_E)*v_y;
	    arr(i,j,k,MOMZ_E) = arr(i,j,k,RHO_E)*v_z;
	    Vector<Real> w_e{arr(i,j,k,RHO_E), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_E) = get_energy(w_e);
	    
	  } else if (test=="LiuLax"){

	    Real rho_1 = 0.5313, v_x_1 = 0.0, v_y_1 = 0.0, v_z_1 = 0.0, p_1 = 0.4;
	    Real rho_2 = 1.0, v_x_2 = 0.7276, v_y_2 = 0.0, v_z_2 = 0.0, p_2 = 1.0;
	    Real rho_3 = 0.8, v_x_3 = 0.0, v_y_3 = 0.0, v_z_3 = 0.0, p_3 = 1.0;
	    Real rho_4 = 1.0, v_x_4 = 0.0, v_y_4 = 0.7276, v_z_4 = 0.0, p_4 = 1.0;
	    
	    if (x>0.5 && y>0.5){
	      rho = rho_1, v_x = v_x_1, v_y = v_y_1, v_z = v_z_1, p = p_1;
	    } else if (x<0.5 && y>0.5){
	      rho = rho_2, v_x = v_x_2, v_y = v_y_2, v_z = v_z_2, p = p_2;
	    } else if (x<0.5 && y<0.5){
	      rho = rho_3, v_x = v_x_3, v_y = v_y_3, v_z = v_z_3, p = p_3;
	    } else if (x>0.5 && y<0.5){
	      rho = rho_4, v_x = v_x_4, v_y = v_y_4, v_z = v_z_4, p = p_4;
	    }
	    
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    // set speed of light
	    c = 100.0;
	    
	  } else if (test=="LiuLax2"){

	    Real rho_1 = 1.0, v_x_1 = 0.75, v_y_1 = -0.5, v_z_1 = 0.0, p_1 = 1.0;
	    Real rho_2 = 2.0, v_x_2 = 0.75, v_y_2 = 0.5, v_z_2 = 0.0, p_2 = 1.0;
	    Real rho_3 = 1.0, v_x_3 = -0.75, v_y_3 = 0.5, v_z_3 = 0.0, p_3 = 1.0;
	    Real rho_4 = 3.0, v_x_4 = -0.75, v_y_4 = -0.5, v_z_4 = 0.0, p_4 = 1.0;
	    
	    if (x>0.5 && y>0.5){
	      rho = rho_1, v_x = v_x_1, v_y = v_y_1, v_z = v_z_1, p = p_1;
	    } else if (x<0.5 && y>0.5){
	      rho = rho_2, v_x = v_x_2, v_y = v_y_2, v_z = v_z_2, p = p_2;
	    } else if (x<0.5 && y<0.5){
	      rho = rho_3, v_x = v_x_3, v_y = v_y_3, v_z = v_z_3, p = p_3;
	    } else if (x>0.5 && y<0.5){
	      rho = rho_4, v_x = v_x_4, v_y = v_y_4, v_z = v_z_4, p = p_4;
	    }
	    
	    arr(i,j,k,0) = rho;
	    arr(i,j,k,1) = arr(i,j,k,0)*v_x;
	    arr(i,j,k,2) = arr(i,j,k,0)*v_y;
	    arr(i,j,k,3) = arr(i,j,k,0)*v_z;
	    Vector<Real> w_i{arr(i,j,k,0), v_x, v_y, v_z, p};
	    arr(i,j,k,ENER_I) = get_energy(w_i);
	    
	    // set speed of light
	    c = 100.0;
	    
	  } else if (test=="EMwave"){

	    for (int n=0; n<NUM_STATE_FLUID; n++)
	      arr(i,j,k,n) = 0.0;
	    
	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = 0.0;
	    arr(i,j,k,BZ) = std::cos(2.0*M_PI*(x+y));
	    arr(i,j,k,EX) = -c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EY) = c*std::cos(2.0*M_PI*(x+y))/std::sqrt(2.0);
	    arr(i,j,k,EZ) = 0.0;

	  } else if (test=="EMwave1d"){

	    for (int n=0; n<NUM_STATE_FLUID; n++)
	      arr(i,j,k,n) = 0.0;
	    
	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = 0.0;
	    arr(i,j,k,BZ) = std::cos(2.0*M_PI*(x));
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = c*std::cos(2.0*M_PI*(x));
	    arr(i,j,k,EZ) = 0.0;

	  } else if (test=="EMwaveTM"){

	    using namespace std::complex_literals;

	    // cylindrical coordinates
	    const Real rCyl = x;
	    const Real zCyl = y;

	    // parameters
	    Real xi0 = 2.40482555769577;
	    Real kr = xi0;
	    Real kz = M_PI;
	    Real omega = std::sqrt(kr*kr + kz*kz);

	    // function
	    auto Bessel0ExpFunc = std::cyl_bessel_j(0,kr*rCyl)*std::exp(1i*kz*zCyl);
	    auto Bessel1ExpFunc = std::cyl_bessel_j(1,kr*rCyl)*std::exp(1i*kz*zCyl);
	    
	    for (int n=0; n<NUM_STATE_FLUID; n++)
	      arr(i,j,k,n) = 0.0;

	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = 0.0;
	    // negative sign for theta components in cyl. coord.
	    arr(i,j,k,BZ) = -std::real(-1i*omega*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
	    arr(i,j,k,EX) = std::real(-1i*kz*kr/(omega*omega-kz*kz)*Bessel1ExpFunc);
	    arr(i,j,k,EY) = std::real(Bessel0ExpFunc);
	    arr(i,j,k,EZ) = 0.0;

	  } else if (test=="EMwaveTE1d"){

	    using namespace std::complex_literals;
	    
	    // cylindrical coordinates
	    const Real rCyl = x;
	    Real L = 2.0*M_PI; //geom.ProbHi()[0]-geom.ProbLo()[0];
	    // parameters
	    Real z1 = 2.40482555769577;
	    Real omega = z1/L;
	    
	    // function
	    auto Bessel0SinFunc = std::cyl_bessel_j(0,omega*rCyl)*std::sin(omega*0.0);
	    auto Bessel1CosFunc = std::cyl_bessel_j(1,omega*rCyl)*std::cos(omega*0.0);
	    
	    for (int n=0; n<NUM_STATE_FLUID; n++)
	      arr(i,j,k,n) = 0.0;

	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = -Bessel0SinFunc;	    
	    arr(i,j,k,BZ) = 0.0;
	    arr(i,j,k,EX) = 0.0;
	    arr(i,j,k,EY) = 0.0;
	    // negative sign for theta components in cyl. coord.
	    arr(i,j,k,EZ) = -Bessel1CosFunc;
	    
	  } else if (test=="gaussianEM"){

	    Real lambda = 1.5, chi = 1.5, a = -2.5, b = -2.5;
	    Real COS = std::cos(2.0*M_PI*(x+y)/lambda);
	    Real SIN = std::sin(2.0*M_PI*(x+y)/lambda);
	    Real FACTOR = ((x-a)*(x-a)+(y-b)*(y-b))/(chi*chi);
	    Real EXP = std::exp(-FACTOR);
	    Real epsilon = 5.0 - 4.0*std::tanh((std::sqrt(x*x+y*y)-0.75)/0.08);
	    for (int n=0; n<NUM_STATE_FLUID; n++)
	      arr(i,j,k,n) = 0.0;
	    
	    arr(i,j,k,BX) = 0.0;
	    arr(i,j,k,BY) = 0.0;
	    arr(i,j,k,BZ) = COS*EXP - lambda*SIN*EXP*(x-a)/(chi*chi*M_PI);
	    arr(i,j,k,EX) = -COS*EXP + lambda*SIN*EXP*(y-b)/(chi*chi*M_PI);
	    arr(i,j,k,EX) *= c/(epsilon*std::sqrt(2.0));
	    arr(i,j,k,EY) = COS*EXP - lambda*SIN*EXP*(x-a)/(chi*chi*M_PI);
	    arr(i,j,k,EY) *= c/(epsilon*std::sqrt(2.0));
	    arr(i,j,k,EZ) = 0.0;

	  } else{
	    amrex::Abort("Please enter valid test in inputs file");
	  }	  
	}
      }
    }
  }

  if (verbose) {
    amrex::Print() << "Done initializing the level " << level 
		   << " data " << std::endl;
  }
}
