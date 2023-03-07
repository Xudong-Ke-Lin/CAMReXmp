#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::fluidSolver(MultiFab& Sborder, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  /*Vector<Real> dim;
  if (geom.Coord()==0)
    dim = {0,1};
  else if (geom.Coord()==1)
  dim = {0,2};  */
  /*ParmParse pp;
  std::string coor;
  pp.query("coor",coor);
  Vector<Real> dim;
  if (coor=="xz")
    dim = {0,2};
  else
    dim = {0,1};
*/
  
  MFIter::allowMultipleMFIters(true);

  int from, to, step;
  if (geom.Coord()==0){
    from = 0;
    to = amrex::SpaceDim;
    step = 1;
  } else {
    from = amrex::SpaceDim-1;
    to = -1;
    step = -1;
  }
  
  for (int d = from; d != to; d += step)
  //for (int d = amrex::SpaceDim-1; d >= 0 ; d--)
  //for (int d = 0; d < amrex::SpaceDim ; d++)   
  {

    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);

    // Loop over all the patches at this level
    for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = Sborder.array(mfi);
      const auto& fluxArr = fluxes[d].array(mfi);
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    //Vector<Real> flux = flux_SLIC(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d, cfl);
	    Vector<Real> flux = fluid_flux_HLLC(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);
	    //Vector<Real> flux = flux_HLLC(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, dim[d]);
	    for(int n=0; n<NUM_STATE_FLUID; n++)
	      //for(int n=0; n<NUM_STATE; n++)
	      {
		fluxArr(i,j,k,n) = flux[n];
	      }	   
	  }
	}
      }
      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    for(int n=0; n<NUM_STATE_FLUID; n++)
	    //for(int n=0; n<NUM_STATE; n++)
              {
		// Conservative update formula
		arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
              }
	  }
	}
      }
    }
    
    // We need to compute boundary conditions again after each update
    Sborder.FillBoundary(geom.periodicity());

    // added by 2020D 
    // Fill non-periodic physical boundaries                      
    FillDomainBoundary(Sborder, geom, bc);
    
    // The fluxes now need scaling for the reflux command.
    // This scaling is by the size of the boundary through which the flux passes, e.g. the x-flux needs scaling by the dy, dz and dt
    if(do_reflux)
    {
      Real scaleFactor = dt;
      for(int scaledir = 0; scaledir < amrex::SpaceDim; ++scaledir)
      {
	// Fluxes don't need scaling by dx[d]
	if(scaledir == d)
	{
	  continue;
	}
	scaleFactor *= dx[scaledir];
      }

      fluxes[d].mult(scaleFactor, 0, NUM_STATE);
    }
  }  
}
void CAMReXmp::fluidSolverWENO(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  MultiFab Slopes(grids, dmap, NUM_STATE_FLUID*nSlopes, NUM_GROW);

  const int iOffset = 1;
  const int jOffset = ( amrex::SpaceDim > 1 ? 1 : 0);
  const int kOffset = ( amrex::SpaceDim == 3 ? 1 : 0);

  // Compute slopes for WENO reconstruction
  for (MFIter mfi(Slopes, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_source.array(mfi);
      const auto& slopes = Slopes.array(mfi);

      for(int k = lo.z-kOffset; k <= hi.z+kOffset; k++)
	{
	  for(int j = lo.y-jOffset; j <= hi.y+jOffset; j++)
	    {
	      for(int i = lo.x-iOffset; i <= hi.x+iOffset; i++)
		{
		  for (int n = 0; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataX = WENO_data(arr, i, j, k, 2, 0, 0, n);
		      std::array<Real, 2> slopesX = WENO3_slope(dataX);
		      slopes(i,j,k,n) = slopesX[0];		      
		      slopes(i,j,k,n+NUM_STATE_FLUID) = slopesX[1];
		    }
#if (AMREX_SPACEDIM >= 2)
		  for (int n = 0; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataY = WENO_data(arr, i, j, k, 0, 2, 0, n);
		      std::array<Real, 2> slopesY = WENO3_slope(dataY);
		      slopes(i,j,k,n+2*NUM_STATE_FLUID) = slopesY[0];
		      slopes(i,j,k,n+3*NUM_STATE_FLUID) = slopesY[1];		      
		    }

		  for (int n = 0; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataXY = WENO_data(arr, i, j, k, 1, 1, 0, n);
		      Real slopesCross = WENO3_slopeCross(dataXY, {slopes(i,j,k,n),slopes(i,j,k,n+NUM_STATE_FLUID),
								   slopes(i,j,k,n+2*NUM_STATE_FLUID),slopes(i,j,k,n+3*NUM_STATE_FLUID)});
		      slopes(i,j,k,n+4*NUM_STATE_FLUID) = slopesCross;
		    }
#endif
		  
		}
	    }
	}
    }

  // Update cell-centred fluid variables and z-components of EM fields
  for (int d = 0; d < amrex::SpaceDim ; d++)   
  {
    
    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);

    // Loop over all the patches at this level
    for (MFIter mfi(S_source, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_source.array(mfi);
      const auto& fluxArr = fluxes[d].array(mfi);
      //const auto& arrEM = S_EM[d].array(mfi);      


      const auto& slopes = Slopes.array(mfi);
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = WENO_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					    0, NUM_STATE_FLUID/2, dx[d], dt, d,
					    fluidFlux, HLLC);
	    Vector<Real> flux_e = WENO_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
							  NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, dx[d], dt, d,
							  fluidFlux, HLLC);	    
	    //Vector<Real> flux = fluid_flux_HLLC(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);

	    for(int n=0; n<NUM_STATE_FLUID/2; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
		fluxArr(i,j,k,n+NUM_STATE_FLUID/2) = flux_e[n];
	      }	    
	  }
	}
      }
      /*
      for(int k = lo.z; k <= hi.z; k++)
      {
      	for(int j = lo.y; j <= hi.y; j++)
      	{
      	  for(int i = lo.x; i <= hi.x; i++)
      	  {
      	    // Update fluid variables
      	    //for(int n=0; n<NUM_STATE; n++)
      	    for(int n=0; n<NUM_STATE_FLUID; n++)
              {
      		// Conservative update formula
      		arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
              }
      	    // Initialise to zero
      	    arr(i,j,k,DIVB) = 0.0;
      	    arr(i,j,k,DIVE) = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
      	  }
      	}
	}*/
    }
      
    // We need to compute boundary conditions again after each update
    //S_source.FillBoundary(geom.periodicity());
     
    // added by 2020D 
    // Fill non-periodic physical boundaries
    //FillDomainBoundary(S_source, geom, bc);
    
    // The fluxes now need scaling for the reflux command.
    // This scaling is by the size of the boundary through which the flux passes, e.g. the x-flux needs scaling by the dy, dz and dt
    /*if(do_reflux)
    {
      Real scaleFactor = dt;
      for(int scaledir = 0; scaledir < amrex::SpaceDim; ++scaledir)
      {
	// Fluxes don't need scaling by dx[d]
	if(scaledir == d)
	{
	  continue;
	}
	scaleFactor *= dx[scaledir];
      }
      // The mult function automatically multiplies entries in a multifab by a scalar
      // scaleFactor: The scalar to multiply by
      // 0: The first data index in the multifab to multiply
      // NUM_STATE:  The total number of data indices that will be multiplied
      fluxes[d].mult(scaleFactor, 0, NUM_STATE);
      }*/
  }
  
  // Unsplit hyperbolic update  
  /*
  for (int d = 0; d < amrex::SpaceDim ; d++)   
  {

    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);
  */
    // Loop over all the patches at this level
    for (MFIter mfi(S_dest, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_dest.array(mfi);
      const auto& arrOld = S_source.array(mfi);
      const auto& fluxArrX = fluxes[0].array(mfi);
      const auto& fluxArrY = fluxes[1].array(mfi);      
  
      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    for(int n=0; n<NUM_STATE_FLUID; n++)
              {		
		// Conservative update formula
		//arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
		arr(i,j,k,n) = arrOld(i,j,k,n) - (dt / dx[0]) * (fluxArrX(i+1, j, k, n) - fluxArrX(i,j,k,n)) - (dt / dx[1]) * (fluxArrY(i, j+1, k, n) - fluxArrY(i,j,k,n)) ;
              }
	    // Update cell-centred z-components becuause it is 2D code
	    //arr(i,j,k,BZ) = arr(i,j,k,BZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, BZ) - fluxArr(i,j,k,BZ));
	    //arr(i,j,k,EZ) = arr(i,j,k,EZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, EZ) - fluxArr(i,j,k,EZ));
	    
	    // Initialise to zero
	    //arr(i,j,k,DIVB) = 0.0;
	    //arr(i,j,k,DIVE) = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
	  }
	}
      }    
    }    
    //}  

  // We need to compute boundary conditions again after each update
  S_dest.FillBoundary(geom.periodicity());
     
  // added by 2020D 
  // Fill non-periodic physical boundaries                      
  FillDomainBoundary(S_dest, geom, bc);  

}
