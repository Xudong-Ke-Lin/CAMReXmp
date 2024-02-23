#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::fluidSolverExact(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
    
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
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = exact_flux(arr, i, j, k, iOffset, jOffset, kOffset,
					     0, NUM_STATE_FLUID/2, dx[d], dt, d);	    
	      
	    for(int n=0; n<NUM_STATE_FLUID/2; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
		//fluxArr(i,j,k,n+NUM_STATE_FLUID/2) = flux_e[n];
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
#if (AMREX_SPACEDIM >= 2)
      const auto& fluxArrY = fluxes[1].array(mfi);      
#endif

      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    //for(int n=0; n<NUM_STATE_FLUID; n++)
	    for(int n=0; n<NUM_STATE_FLUID/2; n++)
              {		
		// Conservative update formula
		//arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
		arr(i,j,k,n) = arrOld(i,j,k,n) - (dt / dx[0]) * (fluxArrX(i+1, j, k, n) - fluxArrX(i,j,k,n));
#if (AMREX_SPACEDIM >= 2)
		arr(i,j,k,n) -= (dt / dx[1]) * (fluxArrY(i, j+1, k, n) - fluxArrY(i,j,k,n));
#endif
              }

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
void CAMReXmp::MUSCLHancokFluidSolverTVD(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  
  MultiFab Slopes(grids, dmap, NUM_STATE_FLUID*2, NUM_GROW);

  const int iOffset = 1;
  const int jOffset = ( amrex::SpaceDim > 1 ? 1 : 0);
  const int kOffset = ( amrex::SpaceDim == 3 ? 1 : 0);

  // Compute slopes for TVD reconstruction
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
		  Vector<Real> limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_I);
		  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      Real slopesX = TVD_slope(dataX,limiterX);
		      slopes(i,j,k,n) = slopesX;		      
		    }
		  limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_E);
		  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      Real slopesX = TVD_slope(dataX,limiterX);
		      slopes(i,j,k,n) = slopesX;		      
		    }		  
#if (AMREX_SPACEDIM >= 2)
		  Vector<Real> limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_I);
		  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE_FLUID) = slopesY;		      
		    }
		  limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_E);
		  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE_FLUID) = slopesY;		      
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

      const auto& slopes = Slopes.array(mfi);
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = MUSCL_Hancock_TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
							 0, NUM_STATE_FLUID/2, dx[d], dt, d, fluidFlux, HLLC);	    
	    Vector<Real> flux_e = MUSCL_Hancock_TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
							 NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, dx[d], dt, d, fluidFlux, HLLC);	    

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
#if (AMREX_SPACEDIM >= 2)
      const auto& fluxArrY = fluxes[1].array(mfi);      
#endif

      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    for(int n=0; n<NUM_STATE_FLUID; n++)
	      //for(int n=0; n<NUM_STATE_FLUID/2; n++)
              {		
		// Conservative update formula
		//arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
		arr(i,j,k,n) = arrOld(i,j,k,n) - (dt / dx[0]) * (fluxArrX(i+1, j, k, n) - fluxArrX(i,j,k,n));
#if (AMREX_SPACEDIM >= 2)
		arr(i,j,k,n) -= (dt / dx[1]) * (fluxArrY(i, j+1, k, n) - fluxArrY(i,j,k,n));
#endif
              }
	    // Initialise to zero
	    arr(i,j,k,DIVB) = 0.0;
	    arr(i,j,k,DIVE) = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
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
void CAMReXmp::fluidSolverNothing(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  MultiFab::Copy(S_dest,S_source,0,0,NUM_STATE_FLUID,0);
  fluxes[0] = 0.0;
  fluxes[1] = 0.0;
  return;
}
void CAMReXmp::fluidSolverTVD(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  
  MultiFab Slopes(grids, dmap, NUM_STATE_FLUID*2, NUM_GROW);

  const int iOffset = 1;
  const int jOffset = ( amrex::SpaceDim > 1 ? 1 : 0);
  const int kOffset = ( amrex::SpaceDim == 3 ? 1 : 0);

  // Compute slopes for TVD reconstruction
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
		  Vector<Real> limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_I);
		  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      Real slopesX = TVD_slope(dataX,limiterX);
		      slopes(i,j,k,n) = slopesX;		      
		    }
		  limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_E);
		  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      Real slopesX = TVD_slope(dataX,limiterX);
		      slopes(i,j,k,n) = slopesX;		      
		    }
#if (AMREX_SPACEDIM >= 2)
		  Vector<Real> limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_I);
		  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE_FLUID) = slopesY;		      
		    }
		  limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_E);
		  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE_FLUID) = slopesY;		      
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

      const auto& slopes = Slopes.array(mfi);
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					   0, NUM_STATE_FLUID/2, d*NUM_STATE_FLUID,
					   dx[d], dt, d, HLLC);
	    Vector<Real> flux_e = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					   NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, d*NUM_STATE_FLUID+NUM_STATE_FLUID/2,
					   dx[d], dt, d, HLLC);

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
#if (AMREX_SPACEDIM >= 2)
      const auto& fluxArrY = fluxes[1].array(mfi);
#endif
  
      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    for(int n=0; n<NUM_STATE_FLUID; n++)
	      //for(int n=0; n<NUM_STATE_FLUID/2; n++)
              {		
		// Conservative update formula
		//arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
		arr(i,j,k,n) = arrOld(i,j,k,n) - (dt / dx[0]) * (fluxArrX(i+1, j, k, n) - fluxArrX(i,j,k,n));
#if (AMREX_SPACEDIM >= 2)
		arr(i,j,k,n) -= (dt / dx[1]) * (fluxArrY(i, j+1, k, n) - fluxArrY(i,j,k,n));
#endif
              }
	    // Update cell-centred z-components becuause it is 2D code
	    //arr(i,j,k,BZ) = arr(i,j,k,BZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, BZ) - fluxArr(i,j,k,BZ));
	    //arr(i,j,k,EZ) = arr(i,j,k,EZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, EZ) - fluxArr(i,j,k,EZ));
	    
	    // Initialise to zero
	    arr(i,j,k,DIVB) = 0.0;
	    arr(i,j,k,DIVE) = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
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
void CAMReXmp::fluidSolverWENO(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  MultiFab Slopes(grids, dmap, NUM_STATE_FLUID*5, NUM_GROW);

  const int iDomainOffset = 1;
  const int jDomainOffset = ( amrex::SpaceDim > 1 ? 1 : 0);
  const int kDomainOffset = ( amrex::SpaceDim == 3 ? 1 : 0);

  // Compute slopes for WENO reconstruction
  for (int d = 0; d < amrex::SpaceDim ; d++)   
    {
      
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
      
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

	  for(int k = lo.z-kDomainOffset; k <= hi.z+kDomainOffset; k++)
	    {
	      for(int j = lo.y-jDomainOffset; j <= hi.y+jDomainOffset; j++)
		{
		  for(int i = lo.x-iDomainOffset; i <= hi.x+iDomainOffset; i++)
		    {
		      /*//for (int n = 0; n<NUM_STATE_FLUID; n++)
		      for(int n=0; n<NUM_STATE_FLUID/2; n++)
			{
			  Vector<Real> dataX = get_data_stencil(arr, i, j, k, 2*iOffset, 2*jOffset, 2*kOffset, n);
			  //Vector<Real> dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, n);
			  std::array<Real, 2> slopesX = WENO3_slope(dataX);		      
			  slopes(i,j,k,n+2*NUM_STATE_FLUID*d) = slopesX[0];
			  //slopes(i,j,k,n) = slopesX[0];
			  slopes(i,j,k,n+NUM_STATE_FLUID+2*NUM_STATE_FLUID*d) = slopesX[1];
			  //slopes(i,j,k,n+NUM_STATE_FLUID) = slopesX[1];
			}
		      */		      
		      //////////////////////////////////////////////////
		      // local characteristic reconstruction
		      WENOcharacteristic(arr,slopes,i,j,k,iOffset,jOffset,kOffset,0,NUM_STATE_FLUID/2,d);
		      WENOcharacteristic(arr,slopes,i,j,k,iOffset,jOffset,kOffset,NUM_STATE_FLUID/2,NUM_STATE_FLUID/2,d);
		      //////////////////////////////////////////////////
		    }
		}
	    }
	}
    }

#if (AMREX_SPACEDIM >= 2)
  for (MFIter mfi(Slopes, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      
      const Dim3 hiDomain = ubound(geom.Domain()); 
      
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_source.array(mfi);
      const auto& slopes = Slopes.array(mfi);

      for(int k = lo.z-kDomainOffset; k <= hi.z+kDomainOffset; k++)
	{
	  for(int j = lo.y-jDomainOffset; j <= hi.y+jDomainOffset; j++)
	    {
	      for(int i = lo.x-iDomainOffset; i <= hi.x+iDomainOffset; i++)
		{
		  for (int n = 0; n<NUM_STATE_FLUID; n++)
		    //for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataXY = get_data_stencil(arr, i, j, k, 1, 1, 0, n);
		      Real slopesCross = WENO3_slopeCross(dataXY, {slopes(i,j,k,n),slopes(i,j,k,n+NUM_STATE_FLUID),
			slopes(i,j,k,n+2*NUM_STATE_FLUID),slopes(i,j,k,n+3*NUM_STATE_FLUID)});
		      slopes(i,j,k,n+4*NUM_STATE_FLUID) = slopesCross;
		      
		      ////////////////////////////////////////////////////////////////////// 
		      if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
			{
			  Vector<Real> limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_I);
			  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
			    {
			      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
			      Real slopesX = TVD_slope(dataX,limiterX);
			      slopes(i,j,k,n) = slopesX;
			      slopes(i,j,k,n+NUM_STATE_FLUID) = 0.0;
			      slopes(i,j,k,n+4*NUM_STATE_FLUID) = 0.0;
			    }
			  limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_E);
			  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
			    {
			      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
			      Real slopesX = TVD_slope(dataX,limiterX);
			      slopes(i,j,k,n) = slopesX;
			      slopes(i,j,k,n+NUM_STATE_FLUID) = 0.0;
			      slopes(i,j,k,n+4*NUM_STATE_FLUID) = 0.0;
			    }
			}
		      if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
			{
			  Vector<Real> limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_I);
			  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
			    {
			      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
			      Real slopesY = TVD_slope(dataY,limiterY);
			      slopes(i,j,k,n+2*NUM_STATE_FLUID) = slopesY;
			      slopes(i,j,k,n+3*NUM_STATE_FLUID) = 0.0;
			      slopes(i,j,k,n+4*NUM_STATE_FLUID) = 0.0;
			    }
			  limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_E);
			  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
			    {
			      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
			      Real slopesY = TVD_slope(dataY,limiterY);
			      slopes(i,j,k,n+2*NUM_STATE_FLUID) = slopesY;
			      slopes(i,j,k,n+3*NUM_STATE_FLUID) = 0.0;
			      slopes(i,j,k,n+4*NUM_STATE_FLUID) = 0.0;
			    }
			}
		      //////////////////////////////////////////////////////////////////////
		    }
		}
	    }
	}
    } 
#endif
  
  MultiFab positivity(grids, dmap, 2, 1);
  //flattenerAlgorithm(positivity, S_source, Slopes, 0, NUM_STATE_FLUID/2);
  MultiFab positivity_e(grids, dmap, 2, 1);
  //flattenerAlgorithm(positivity_e, S_source, Slopes, NUM_STATE_FLUID/2, NUM_STATE_FLUID/2);
  
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

      const auto& slopes = Slopes.array(mfi);

      const auto& tau = positivity.array(mfi);
      const auto& tau_e = positivity_e.array(mfi);
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = WENO_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					    0, NUM_STATE_FLUID/2, dx[d], dt, d,
					    HLLC);
	    Vector<Real> flux_e = WENO_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					    NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, dx[d], dt, d,
					    HLLC);

	    /*Vector<Real> flux_i = WENO_flux_flat(arr, slopes, tau, i, j, k, iOffset, jOffset, kOffset,
						 0, NUM_STATE_FLUID/2, dx[d], dt, d, HLLC);
	    Vector<Real> flux_e = WENO_flux_flat(arr, slopes, tau_e, i, j, k, iOffset, jOffset, kOffset,
						 NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, dx[d], dt, d, HLLC);
	    */
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
#if (AMREX_SPACEDIM >= 2)
      const auto& fluxArrY = fluxes[1].array(mfi);      
#endif

      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    for(int n=0; n<NUM_STATE_FLUID; n++)
	      //for(int n=0; n<NUM_STATE_FLUID/2; n++)
              {		
		// Conservative update formula
		//arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
		arr(i,j,k,n) = arrOld(i,j,k,n) - (dt / dx[0]) * (fluxArrX(i+1, j, k, n) - fluxArrX(i,j,k,n));
#if (AMREX_SPACEDIM >= 2)
		arr(i,j,k,n) -= (dt / dx[1]) * (fluxArrY(i, j+1, k, n) - fluxArrY(i,j,k,n));
#endif
              }
	    // Update cell-centred z-components becuause it is 2D code
	    //arr(i,j,k,BZ) = arr(i,j,k,BZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, BZ) - fluxArr(i,j,k,BZ));
	    //arr(i,j,k,EZ) = arr(i,j,k,EZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, EZ) - fluxArr(i,j,k,EZ));
	    
	    // Initialise to zero
	    arr(i,j,k,DIVB) = 0.0;
	    arr(i,j,k,DIVE) = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
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
void CAMReXmp::flattenerAlgorithm(MultiFab& positivity, MultiFab& S_source, MultiFab& Slopes, int start, int len){
  
  const int ydim = ( amrex::SpaceDim >= 2 ? 1 : 0);
  const int zdim = ( amrex::SpaceDim == 3 ? 1 : 0);
  
  Real kappa1 = 0.4;
  Real kappa2 = 0.4;

  // Calculate flattener
  MultiFab flattener(grids, dmap, 1, 1);
  // Zones that have entered the shock
  for (MFIter mfi(S_source, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();
      
      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arr = S_source.array(mfi);
      const auto& eta = flattener.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Real c_min_nbr = 1e14;
		  // loop the neighbours for each cell
		  for (int iloc = -1; iloc <= 1; iloc++)   
		    {
		      for (int jloc = -1*ydim; jloc <= 1*ydim; jloc++)
			{
			  for (int kloc = -1*zdim; kloc <= 1*zdim; kloc++)
			    {
			      Vector<Real> u_i = get_data_zone(arr,i+iloc,j+jloc,k+kloc,start,start+len);
			      c_min_nbr = std::min(c_min_nbr,get_speed(u_i));
			    }		      		      			  
			}		      		      
		    }

		  // calculate divergence of velocity
		  Real divV = 0.0;
		  for (int d = 0; d < amrex::SpaceDim ; d++)
		    {
		      const int iOffset = ( d == 0 ? 1 : 0);
		      const int jOffset = ( d == 1 ? 1 : 0);
		      const int kOffset = ( d == 2 ? 1 : 0);

		      divV += (arr(i+iOffset,j+jOffset,k+kOffset,1+d+start)/arr(i+iOffset,j+jOffset,k+kOffset,start))
			- (arr(i-iOffset,j-jOffset,k-kOffset,1+d+start)/arr(i-iOffset,j-jOffset,k-kOffset,start)); 
		    }
		  // calculate flattener
		  eta(i,j,k,0) = std::min(1.0,std::max(0.0,-(divV+kappa1*c_min_nbr)/(kappa1*c_min_nbr)));
		}
	    }
	}
    }

  flattener.FillBoundary(geom.periodicity());
  FillDomainBoundary(flattener, geom, {bc[start]});  

  // Zones that are about to be run over by a shock
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

      const auto& arr = S_source.array(mfi);
      const auto& eta = flattener.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
      {
  	for(int j = lo.y; j <= hi.y; j++)
  	{
  	  for(int i = lo.x; i <= hi.x; i++)
  	  {
  	    Vector<Real> u_iMinus1 = get_data_zone(arr,i-iOffset,j-jOffset,k-kOffset,start,start+len);
  	    Vector<Real> u_i = get_data_zone(arr,i,j,k,start,start+len);
  	    Vector<Real> u_iPlus1 = get_data_zone(arr,i+iOffset,j+jOffset,k+kOffset,start,start+len);
  	    Real p_iMinus1 = get_pressure(u_iMinus1);
  	    Real p_i = get_pressure(u_i);
  	    Real p_iPlus1 = get_pressure(u_iPlus1);
  	    if ((eta(i,j,k,0)>0.0 && eta(i+iOffset,j+jOffset,k+kOffset,0)==0.0) && (p_i>p_iPlus1))
  	      {
  		//std::cout << "Next about to be run " << i << std::endl;
  		eta(i+iOffset,j+jOffset,k+kOffset,0) = eta(i,j,k,0);
  	      }
  	    if ((eta(i,j,k,0)>0.0 && eta(i-iOffset,j-jOffset,k-kOffset,0)==0.0) && (p_i>p_iMinus1))
  	      {
  		//std::cout << "Previous about to be run " << i << std::endl;
  		eta(i-iOffset,j-jOffset,k-kOffset,0) = eta(i,j,k,0);
  	      }
  	  }
  	}
      }
      /*
      for(int k = hi.z; k >= lo.z; k--)
      {
  	for(int j = hi.y; j >= lo.y; j--)
  	{
  	  for(int i = hi.x; i >= lo.x; i--)
  	  {
  	    Vector<Real> u_iMinus1 = get_data_zone(arr,i-iOffset,j-jOffset,k-kOffset,start,start+len);
  	    Vector<Real> u_i = get_data_zone(arr,i,j,k,start,start+len);
  	    Vector<Real> u_iPlus1 = get_data_zone(arr,i+iOffset,j+jOffset,k+kOffset,start,start+len);
  	    Real p_iMinus1 = get_pressure(u_iMinus1);
  	    Real p_i = get_pressure(u_i);
  	    Real p_iPlus1 = get_pressure(u_iPlus1);
  	    if ((eta(i,j,k,0)>0.0 && eta(i+iOffset,j+jOffset,k+kOffset,0)==0.0) && (p_i>p_iPlus1))
  	      {
  		//std::cout << "Next about to be run " << i << std::endl;
  		eta(i+iOffset,j+jOffset,k+kOffset,0) = eta(i,j,k,0);
  	      }
  	    if ((eta(i,j,k,0)>0.0 && eta(i-iOffset,j-jOffset,k-kOffset,0)==0.0) && (p_i>p_iMinus1))
  	      {
  		//std::cout << "Previous about to be run " << i << std::endl;
  		eta(i-iOffset,j-jOffset,k-kOffset,0) = eta(i,j,k,0);
  	      }
  	  }
  	}
      }
      */
    }
  }
  flattener.FillBoundary(geom.periodicity());
  FillDomainBoundary(flattener, geom, {bc[start]});  
  
  int tau_rho = 0, tau_p = 1;
  
  // Parameter used to restore positivity
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

      const auto& arr = S_source.array(mfi);
      const auto& slopes = Slopes.array(mfi);

      const auto& eta = flattener.array(mfi);
      const auto& tau = positivity.array(mfi);
      
      for(int k = lo.z; k <= hi.z; k++)
      {
  	for(int j = lo.y; j <= hi.y; j++)
  	{
  	  for(int i = lo.x; i <= hi.x; i++)
  	  {
	    // define conserved variables                     
	    Real rho_i = arr(i,j,k,start);
	    Real momX_i = arr(i,j,k,1+start);
	    Real momY_i = arr(i,j,k,2+start);
	    Real momZ_i = arr(i,j,k,MOMZ_I+start);
	    Real E_i = arr(i,j,k,ENER_I+start);
  
	    Real rho_min_nbr = 1e14;
	    Real rho_max_nbr = 1e-14;
	    // loop the neighbours for each cell
	    for (int iloc = -1; iloc <= 1; iloc++)   
	      {
		for (int jloc = -1*ydim; jloc <= 1*ydim; jloc++)
		  {
		    for (int kloc = -1*zdim; kloc <= 1*zdim; kloc++)
		      {
			rho_min_nbr = std::min(rho_min_nbr,arr(i+iloc,j+jloc,k+kloc,start));
			rho_max_nbr = std::max(rho_max_nbr,arr(i+iloc,j+jloc,k+kloc,start));
		      }		      		      			  
		  }		      		      
	      }

  	    // range of density
  	    Real rho_min_extended = rho_min_nbr*(1.0-kappa2+kappa2*eta(i,j,k,0));
  	    Real rho_max_extended = rho_max_nbr*(1.0+kappa2-kappa2*eta(i,j,k,0));
  	    rho_min_extended = std::max(rho_min_extended,1e-14);
	    
	    Real p_min_nbr = 1e14;
	    Real p_max_nbr = 1e-14;
	    // loop the neighbours for each cell
	    for (int iloc = -1; iloc <= 1; iloc++)   
	      {
		for (int jloc = -1*ydim; jloc <= 1*ydim; jloc++)
		  {
		    for (int kloc = -1*zdim; kloc <= 1*zdim; kloc++)
		      {
			Vector<Real> u_i = get_data_zone(arr,i+iloc,j+jloc,k+kloc,start,start+len);
			p_min_nbr = std::min(p_min_nbr,get_pressure(u_i));
			p_max_nbr = std::max(p_max_nbr,get_pressure(u_i));
		      }		      		      			  
		  }		      		      
	      }
  	    // range of pressure
  	    Real p_min_extended = p_min_nbr*(1.0-kappa2+kappa2*eta(i,j,k,0));
  	    Real p_max_extended = p_max_nbr*(1.0+kappa2-kappa2*eta(i,j,k,0));
  	    p_min_extended = std::max(p_min_extended,1e-14);
	    
  	    // slopes
  	    Vector<Real> ux_i,uxx_i;
#if (AMREX_SPACEDIM >= 2)
	    Vector<Real> uy_i,uyy_i,uxy_i;
#endif
	    
  	    for (int n = start; n<start+len; n++)
  	      {
  		ux_i.push_back(slopes(i,j,k,n));
  		uxx_i.push_back(slopes(i,j,k,n+NUM_STATE_FLUID));
#if (AMREX_SPACEDIM >= 2)
		uy_i.push_back(slopes(i,j,k,n+2*NUM_STATE_FLUID));
		uyy_i.push_back(slopes(i,j,k,n+3*NUM_STATE_FLUID));
		uxy_i.push_back(slopes(i,j,k,n+4*NUM_STATE_FLUID));
#endif      		
  	      }
  	    // reconstruction at nodal points
  	    Real rho_min_zone = 1e14, rho_max_zone=1e-14;
	    Vector<Vector<Real>> u_allnod;
	    for (int iloc = -1; iloc <= 1; iloc++)   
	      {
		for (int jloc = -1*ydim; jloc <= 1*ydim; jloc++)
		  {
		    for (int kloc = -1*zdim; kloc <= 1*zdim; kloc++)
		      {
			Vector<Real> u_nod;
			for (int n = 0; n<ux_i.size(); n++)
			  {			    
			    Real x = 0.5*iloc;
			    Real Lx = x, Lxx = x*x - 1.0/12.0;
			    Real var_nod = arr(i,j,k,n+start) + ux_i[n]*Lx + uxx_i[n]*Lxx;

#if (AMREX_SPACEDIM >= 2)
			    Real y = 0.5*jloc;
			    Real Ly = y, Lyy = y*y - 1.0/12.0;
			
			    var_nod +=  uy_i[n]*Ly + uyy_i[n]*Lyy + uxy_i[n]*Lx*Ly;
#endif
			    u_nod.push_back(var_nod);
			  }
			// min and max values of density at nodal points
			rho_min_zone = std::min(rho_min_zone,u_nod[0]);
			rho_max_zone = std::max(rho_max_zone,u_nod[0]);
			u_allnod.push_back(u_nod);
		      }
		  }
	      }
  	    
  	    //Real epsilon = 1e-12;
	    Real ratio_max = (rho_max_extended-rho_i)/(rho_max_zone-rho_i);
	    Real ratio_min = (rho_i-rho_min_extended)/(rho_i-rho_min_zone);
	    
  	    if ((rho_min_zone>=rho_min_extended) //(rho_min_zone>=rho_min_extended-epsilon
  		&& (rho_max_zone<=rho_max_extended)) // rho_max_zone<=rho_max_extended+epsilon
  	      tau(i,j,k,tau_rho) = 1.0;	      
  	    else
  	      {
  		tau(i,j,k,tau_rho) = std::min((rho_max_extended-rho_i)/(rho_max_zone-rho_i),(rho_i-rho_min_extended)/(rho_i-rho_min_zone));
  		tau(i,j,k,tau_rho) = std::max(tau(i,j,k,tau_rho), 0.0);    
  	      }
  	    // changes for the positivity of density	    
	    for (int n_nod = 0; n_nod<u_allnod.size(); n_nod++)
	      {
		for (int n = 0; n<u_allnod[n_nod].size(); n++)
		  {
		    u_allnod[n_nod][n] = (1.0-tau(i,j,k,tau_rho))*arr(i,j,k,n+start) + tau(i,j,k,tau_rho)*u_allnod[n_nod][n];		    
		  }
	      }
	    
	    Vector <Real> p, rho, momx, momy, momz, energy;
	    for (int n_nod = 0; n_nod<u_allnod.size(); n_nod++)
	      {
		Real p_nod = get_pressure(u_allnod[n_nod]);
		p.push_back(p_nod);
		rho.push_back(u_allnod[n_nod][0]);
		momx.push_back(u_allnod[n_nod][1]);
		momy.push_back(u_allnod[n_nod][2]);
		momz.push_back(u_allnod[n_nod][3]);
		energy.push_back(u_allnod[n_nod][4]);
	      }
  	    Vector <Real> tau_p_node;
  	    for (int n=0; n<u_allnod.size(); n++){
  	      if ((p[n]>=p_min_extended)
  		  && (p[n]<=p_max_extended))	      
  		tau_p_node.push_back(1.0);
  	      else
  		{
  		  Real e_min_extended = p_min_extended/(Gamma-1.0);
  		  // quadratic equation ax^2 + bx + c = 0
  		  Real A = 2.0*(rho[n]-rho_i)*(energy[n]-E_i) -
		    get_magnitude_squared(momx[n]-momX_i,momy[n]-momY_i,momz[n]-momZ_i);		  
  		  Real B = 2.0*rho_i*(energy[n]-E_i)+2.0*E_i*(rho[n]-rho_i) -
  		    2.0*((momx[n]-momX_i)*momX_i+(momy[n]-momY_i)*momY_i+
  			 (momz[n]-momZ_i)*momZ_i) -2.0*e_min_extended*(rho[n]-rho_i);
  		  Real C = 2.0*rho_i*E_i-get_magnitude_squared(momX_i,momY_i,momZ_i)
  		    -2.0*e_min_extended*rho_i;
  		  Real discriminant = B*B - 4.0*A*C;

  		  //if (discriminant<-epsilon)
		  if (discriminant<0.0)
		    {
		      std::cout << "Negative square root at " << i << " " << j << " " << k << " with A=" << A << " B=" << B << " C=" << C << std::endl;
		      std::cout << "Cell averages " << rho_i << " " << momX_i << " " << momY_i << " " << momZ_i << " " << E_i << std::endl;
		      std::cout << "Nodal values " << rho[n] << " " << momx[n] << " " << momy[n] << " " << momz[n] << " " << energy[n] << std::endl;
		      std::cout << "e_min_extended " << e_min_extended << std::endl;
		      
		      for (int iloc = -1; iloc <= 1; iloc++)   
			{
			  for (int jloc = -1*ydim; jloc <= 1*ydim; jloc++)
			    {
			      for (int kloc = -1*zdim; kloc <= 1*zdim; kloc++)
				{
				  Vector<Real> u_i = get_data_zone(arr,i+iloc,j+jloc,k+kloc,start,start+len);
				  std::cout << get_pressure(u_i) << " ";
				  if (iloc==0 && jloc==0 && get_pressure(u_i)<0.0)
				    amrex::Abort();
				}		      		      			  
			    }		      		      
			}
		      std::cout << std::endl;

		      tau_p_node.push_back(0.0);
		      
		      amrex::Abort("Negative square root");
		    }
		  else
		    tau_p_node.push_back(std::min((-B + std::sqrt(discriminant))/(2.0*A), (-B - std::sqrt(discriminant))/(2.0*A)));
		  
  		}
  	    }
  	    tau(i,j,k,tau_p) = *std::min_element(tau_p_node.begin(), tau_p_node.end());
  	    tau(i,j,k,tau_p) = std::max(tau(i,j,k,tau_p),0.0);
  	  }
  	}
      }
    }
  }   
  positivity.FillBoundary(geom.periodicity());
  FillDomainBoundary(positivity, geom, {bc[0],bc[0]});  

}
// void CAMReXmp::flattenerAlgorithm(MultiFab& positivity, MultiFab& S_source, MultiFab& Slopes){
//   // Calculate flattener
//   Real kappa = 0.4;
//   MultiFab flattener(grids, dmap, 1, 1);
//   // Zones that have entered the shock
//   for (int d = 0; d < amrex::SpaceDim ; d++)   
//   {
    
//     const int iOffset = ( d == 0 ? 1 : 0);
//     const int jOffset = ( d == 1 ? 1 : 0);
//     const int kOffset = ( d == 2 ? 1 : 0);

//     // Loop over all the patches at this level
//     for (MFIter mfi(S_source, true); mfi.isValid(); ++mfi)
//     {
//       const Box& bx = mfi.tilebox();

//       const Dim3 lo = lbound(bx);
//       const Dim3 hi = ubound(bx);

//       const auto& arr = S_source.array(mfi);
//       const auto& eta = flattener.array(mfi);
      
//       for(int k = lo.z; k <= hi.z; k++)
//       {
//   	for(int j = lo.y; j <= hi.y; j++)
//   	{
//   	  for(int i = lo.x; i <= hi.x; i++)
//   	  {
//   	    Vector<Real> u_iMinus1 = get_data_zone(arr,i-1,j,k,0,5);
//   	    Vector<Real> u_i = get_data_zone(arr,i,j,k,0,5);
//   	    Vector<Real> u_iPlus1 = get_data_zone(arr,i+1,j,k,0,5);
//   	    Real c_min_nbr = std::min({get_speed(u_iMinus1),get_speed(u_i),get_speed(u_iPlus1)});
//   	    Real divV = (arr(i+1,j,k,1)/arr(i+1,j,k,0)) - (arr(i-1,j,k,1)/arr(i-1,j,k,0)); 
//   	    eta(i,j,k,0) = std::min(1.0,std::max(0.0,-(divV+kappa*c_min_nbr)/(kappa*c_min_nbr)));
//   	  }
//   	}
//       }
//     }
//   }
//   flattener.FillBoundary(geom.periodicity());
//   FillDomainBoundary(flattener, geom, {bc[0]});  

//   // Zones that are abount to be run over by a shock
//   for (int d = 0; d < amrex::SpaceDim ; d++)   
//   {
    
//     const int iOffset = ( d == 0 ? 1 : 0);
//     const int jOffset = ( d == 1 ? 1 : 0);
//     const int kOffset = ( d == 2 ? 1 : 0);

//     // Loop over all the patches at this level
//     for (MFIter mfi(S_source, true); mfi.isValid(); ++mfi)
//     {
//       const Box& bx = mfi.tilebox();

//       const Dim3 lo = lbound(bx);
//       const Dim3 hi = ubound(bx);

//       const auto& arr = S_source.array(mfi);
//       const auto& eta = flattener.array(mfi);
      
//       for(int k = lo.z; k <= hi.z; k++)
//       {
//   	for(int j = lo.y; j <= hi.y; j++)
//   	{
//   	  for(int i = lo.x; i <= hi.x; i++)
//   	  {
//   	    Vector<Real> u_iMinus1 = get_data_zone(arr,i-1,j,k,0,5);
//   	    Vector<Real> u_i = get_data_zone(arr,i,j,k,0,5);
//   	    Vector<Real> u_iPlus1 = get_data_zone(arr,i+1,j,k,0,5);
//   	    Real p_iMinus1 = get_pressure(u_iMinus1);
//   	    Real p_i = get_pressure(u_i);
//   	    Real p_iPlus1 = get_pressure(u_iPlus1);
//   	    /*if (eta(i,j,k,0)>0.0)
//   	      std::cout << "Eta " << i << " " << eta(i-1,j,k,0) << " " << eta(i,j,k,0) << " " << eta(i+1,j,k,0) << " " << p_iMinus1 << " " << p_i << " " << p_iPlus1 << std::endl;*/
//   	    if ((eta(i,j,k,0)>0.0 && eta(i+1,j,k,0)==0.0) && (p_i>p_iPlus1))
//   	      {
//   		//std::cout << "Next about to be run " << i << std::endl;
//   		eta(i+1,j,k,0) = eta(i,j,k,0);
//   	      }
//   	    if ((eta(i,j,k,0)>0.0 && eta(i-1,j,k,0)==0.0) && (p_i>p_iMinus1))
//   	      {
//   		//std::cout << "Previous about to be run " << i << std::endl;
//   		eta(i-1,j,k,0) = eta(i,j,k,0);
//   	      }
//   	  }
//   	}
//       }
//       /*
//       for(int k = hi.z; k >= lo.z; k--)
//       {
//   	for(int j = hi.y; j >= lo.y; j--)
//   	{
//   	  for(int i = hi.x; i >= lo.x; i--)
//   	  {
//   	    Vector<Real> u_iMinus1 = get_data_zone(arr,i-1,j,k,0,5);
//   	    Vector<Real> u_i = get_data_zone(arr,i,j,k,0,5);
//   	    Vector<Real> u_iPlus1 = get_data_zone(arr,i+1,j,k,0,5);
//   	    Real p_iMinus1 = get_pressure(u_iMinus1);
//   	    Real p_i = get_pressure(u_i);
//   	    Real p_iPlus1 = get_pressure(u_iPlus1);
//   	    if (eta(i,j,k,0)>0.0)
//   	      std::cout << "Eta " << i << " " << eta(i-1,j,k,0) << " " << eta(i,j,k,0) << " " << eta(i+1,j,k,0) << " " << p_iMinus1 << " " << p_i << " " << p_iPlus1 << std::endl;
//   	    if ((eta(i,j,k,0)>0.0 && eta(i+1,j,k,0)==0.0) && (p_i>p_iPlus1))
//   	      {
//   		std::cout << "Next2 about to be run " << i << std::endl;
//   		eta(i+1,j,k,0) = eta(i,j,k,0);
//   	      }
//   	    if ((eta(i,j,k,0)>0.0 && eta(i-1,j,k,0)==0.0) && (p_i>p_iMinus1))
//   	      {
//   		std::cout << "Previous2 about to be run " << i << std::endl;
//   		eta(i-1,j,k,0) = eta(i,j,k,0);
//   	      }
//   	  }
//   	}
//   	}*/
//     }
//   }
//   flattener.FillBoundary(geom.periodicity());
//   FillDomainBoundary(flattener, geom, {bc[0]});  
  
//   int tau_rho = 0, tau_p = 1;
  
//   // Parameter used to restore positivity
//   for (int d = 0; d < amrex::SpaceDim ; d++)   
//   {
    
//     const int iOffset = ( d == 0 ? 1 : 0);
//     const int jOffset = ( d == 1 ? 1 : 0);
//     const int kOffset = ( d == 2 ? 1 : 0);

//     // Loop over all the patches at this level
//     for (MFIter mfi(S_source, true); mfi.isValid(); ++mfi)
//     {
//       const Box& bx = mfi.tilebox();

//       const Dim3 lo = lbound(bx);
//       const Dim3 hi = ubound(bx);

//       const auto& arr = S_source.array(mfi);
//       const auto& slopes = Slopes.array(mfi);

//       const auto& eta = flattener.array(mfi);
//       const auto& tau = positivity.array(mfi);
      
//       for(int k = lo.z; k <= hi.z; k++)
//       {
//   	for(int j = lo.y; j <= hi.y; j++)
//   	{
//   	  for(int i = lo.x; i <= hi.x; i++)
//   	  {
//   	    // range of density
//   	    Real rho_min_nbr = std::min({arr(i-1,j,k,0),arr(i,j,k,0),arr(i+1,j,k,0)});
//   	    Real rho_max_nbr = std::max({arr(i-1,j,k,0),arr(i,j,k,0),arr(i+1,j,k,0)});
//   	    Real rho_min_extended = rho_min_nbr*(1.0-kappa+kappa*eta(i,j,k,0));
//   	    Real rho_max_extended = rho_max_nbr*(1.0+kappa-kappa*eta(i,j,k,0));
//   	    rho_min_extended = std::max(rho_min_extended,1e-14);
	    
//   	    // range of pressure
//   	    Vector<Real> u_iMinus1 = get_data_zone(arr,i-1,j,k,0,5);
//   	    Vector<Real> u_i = get_data_zone(arr,i,j,k,0,5);
//   	    Vector<Real> u_iPlus1 = get_data_zone(arr,i+1,j,k,0,5);
//   	    Real p_iMinus1 = get_pressure(u_iMinus1);
//   	    Real p_i = get_pressure(u_i);
//   	    Real p_iPlus1 = get_pressure(u_iPlus1);
//   	    Real p_min_nbr = std::min({p_iMinus1,p_i,p_iPlus1});
//   	    Real p_max_nbr = std::max({p_iMinus1,p_i,p_iPlus1});
//   	    Real p_min_extended = p_min_nbr*(1.0-kappa+kappa*eta(i,j,k,0));
//   	    Real p_max_extended = p_max_nbr*(1.0+kappa-kappa*eta(i,j,k,0));
//   	    p_min_extended = std::max(p_min_extended,1e-14);
	    
//   	    // slopes
//   	    Vector<Real> ux_i,uxx_i;
//   	    for (int n = 0; n<5; n++)
//   	      {
//   		ux_i.push_back(slopes(i,j,k,n));
//   		uxx_i.push_back(slopes(i,j,k,n+NUM_STATE_FLUID));
//   	      }
//   	    // reconstruction at nodal points
//   	    Vector<Real> u_iL,u_iC,u_iR;
//   	    for (int n = 0; n<5; n++)
//   	      {
//   		Real x = -0.5;
//   		Real Lx = x, Lxx = x*x - 1.0/12.0;
//   		u_iL.push_back(arr(i,j,k,n) + ux_i[n]*Lx + uxx_i[n]*Lxx);

//   		x = 0.0;
//   		Lx = x, Lxx = x*x - 1.0/12.0;
//   		u_iC.push_back(arr(i,j,k,n) + ux_i[n]*Lx + uxx_i[n]*Lxx);
//   		x = 0.5;
//   		Lx = x, Lxx = x*x - 1.0/12.0;
//   		u_iR.push_back(arr(i,j,k,n) + ux_i[n]*Lx + uxx_i[n]*Lxx);
//   	      }
//   	    // min and max values of density at nodal points
//   	    Real rho_min_zone = std::min({u_iL[0],u_iC[0],u_iR[0]});
//   	    Real rho_max_zone = std::max({u_iL[0],u_iC[0],u_iR[0]});
//   	    Real epsilon = 1e-12;
//   	    if ((rho_min_zone>=rho_min_extended-epsilon)
//   		&& (rho_max_zone<=rho_max_extended+epsilon))
//   	      tau(i,j,k,tau_rho) = 1.0;	      
//   	    else
//   	      {
//   		//std::cout << i << " " << (rho_min_zone >= rho_min_extended) << " " << (rho_max_zone <= rho_max_extended) << std::endl;
//   		//std::cout << i << " " << rho_min_zone << " " << rho_min_extended << " " << rho_max_zone << " " <<  rho_max_extended << std::endl;
//   		//std::cout << "density " << i << " " << rho_max_extended << " " << rho_max_zone << " " << rho_min_extended << " " << rho_min_zone << " ";
//   		tau(i,j,k,tau_rho) = std::min((rho_max_extended-arr(i,j,k,0))/(rho_max_zone-arr(i,j,k,0)),(arr(i,j,k,0)-rho_min_extended)/(arr(i,j,k,0)-rho_min_zone));
//   		//std::cout << tau(i,j,k,tau_rho) << std::endl;
//   		tau(i,j,k,tau_rho) = std::max(tau(i,j,k,tau_rho), 0.0);	       
//   		//std::cout << i << " " << tau(i,j,k,tau_rho) << std::endl;
//   	      }
//   	    // changes for the positivity of density
//   	    for (int n = 0; n<5; n++)
//   	      {
//   		u_iL[n] = (1-tau(i,j,k,tau_rho))*arr(i,j,k,n) + tau(i,j,k,tau_rho)*u_iL[n];
//   		u_iC[n] = (1-tau(i,j,k,tau_rho))*arr(i,j,k,n) + tau(i,j,k,tau_rho)*u_iC[n];
//   		u_iR[n] = (1-tau(i,j,k,tau_rho))*arr(i,j,k,n) + tau(i,j,k,tau_rho)*u_iR[n];
//   	      }
//   	    Real p_iL = get_pressure(u_iL);
//   	    Real p_iC = get_pressure(u_iC);
//   	    Real p_iR = get_pressure(u_iR);
//   	    //Real p_min_zone = std::min({p_iL,p_iC,p_iR});
//   	    //Real p_max_zone = std::max({p_iL,p_iC,p_iR});
//   	    Vector <Real> p = {p_iL,p_iC,p_iR};
//   	    Vector <Real> rho = {u_iL[0],u_iC[0],u_iR[0]};
//   	    Vector <Real> momx = {u_iL[1],u_iC[1],u_iR[1]};
//   	    Vector <Real> momy = {u_iL[2],u_iC[2],u_iR[2]};
//   	    Vector <Real> momz = {u_iL[3],u_iC[3],u_iR[3]};
//   	    Vector <Real> energy = {u_iL[4],u_iC[4],u_iR[4]};
//   	    Vector <Real> tau_p_node;
//   	    for (int n=0; n<3; n++){
//   	      if ((p[n]>=p_min_extended-epsilon)
//   		  && (p[n]<=p_max_extended+epsilon))	      
//   		tau_p_node.push_back(1.0);
//   	      else
//   		{
//   		  //std::cout << "presure " << i << " " << n << " " << p[n] << " " << p_min_extended << " " << p_max_extended  <<std::endl;
//   		  Real e_min_extended = p_min_extended/(Gamma-1.0);
//   		  // quadratic equation ax^2 + bx + c = 0
//   		  Real A = 2.0*(rho[n]-arr(i,j,k,0))*(energy[n]-arr(i,j,k,4))-((momx[n]-arr(i,j,k,1))*(momx[n]-arr(i,j,k,1)) +
//   									       (momy[n]-arr(i,j,k,2))*(momy[n]-arr(i,j,k,2)) +
//   									       (momz[n]-arr(i,j,k,3))*(momz[n]-arr(i,j,k,3)));
//   		  Real B = 2.0*arr(i,j,k,0)*(energy[n]-arr(i,j,k,4))+2.0*arr(i,j,k,4)*(rho[n]-arr(i,j,k,0)) -
//   		    2.0*((momx[n]-arr(i,j,k,1))*arr(i,j,k,1)+(momy[n]-arr(i,j,k,2))*arr(i,j,k,2)+
//   			 (momz[n]-arr(i,j,k,3))*arr(i,j,k,3)) -2.0*e_min_extended*(rho[n]-arr(i,j,k,0));
//   		  Real C = 2.0*arr(i,j,k,0)*arr(i,j,k,4)-get_magnitude_squared(arr(i,j,k,1),arr(i,j,k,2),arr(i,j,k,3))
//   		    -2.0*e_min_extended*arr(i,j,k,0);
//   		  Real discriminant = B*B - 4.0*A*C;
//   		  //std::cout << discriminant << std::endl;
//   		  if (discriminant<-epsilon)
//   		    amrex::Abort("Negative square root");
//   		  tau_p_node.push_back(std::min((-B + std::sqrt(discriminant))/(2.0*A), (-B - std::sqrt(discriminant))/(2.0*A)));
		  
//   		  //std::cout << i << " " << tau_p_node[n] << " " << std::max(*std::min_element(tau_p_node.begin(),tau_p_node.end()),0.0) << std::endl;		  
//   		}
//   	    }
//   	    tau(i,j,k,tau_p) = *std::min_element(tau_p_node.begin(), tau_p_node.end());
//   	    tau(i,j,k,tau_p) = std::max(tau(i,j,k,tau_p),0.0);
//   	  }
//   	}
//       }
//     }
//   }   
//   positivity.FillBoundary(geom.periodicity());
//   FillDomainBoundary(positivity, geom, {bc[0],bc[0]});  

// }
void CAMReXmp::fluidSolverVoid(MultiFab& S_source, const Real* dx, Real dt)
{
  amrex::Abort("No fluid solver!!");
}
void CAMReXmp::fluidSolverTVD(MultiFab& S_source, const Real* dx, Real dt)
{

  MultiFab& S_dest = get_new_data(Phi_Type);
  
  // Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes[amrex::SpaceDim];
  for (int j = 0; j < amrex::SpaceDim; j++)
    {
      BoxArray ba = S_dest.boxArray();
      ba.surroundingNodes(j);
      fluxes[j].define(ba, dmap, NUM_STATE, 0);
    }
  
  MultiFab Slopes(grids, dmap, NUM_STATE_FLUID*2, NUM_GROW);

  const int iOffset = 1;
  const int jOffset = ( amrex::SpaceDim > 1 ? 1 : 0);
  const int kOffset = ( amrex::SpaceDim == 3 ? 1 : 0);

  // Compute slopes for TVD reconstruction
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
		  Vector<Real> limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_I);
		  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      Real slopesX = TVD_slope(dataX,limiterX);
		      slopes(i,j,k,n) = slopesX;		      
		    }
		  limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_E);
		  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      Real slopesX = TVD_slope(dataX,limiterX);
		      slopes(i,j,k,n) = slopesX;		      
		    }
#if (AMREX_SPACEDIM >= 2)
		  Vector<Real> limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_I);
		  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE_FLUID) = slopesY;		      
		    }
		  limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_E);
		  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE_FLUID) = slopesY;		      
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

      const auto& slopes = Slopes.array(mfi);
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					   0, NUM_STATE_FLUID/2, d*NUM_STATE_FLUID,
					   dx[d], dt, d, HLLC);
	    Vector<Real> flux_e = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					   NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, d*NUM_STATE_FLUID+NUM_STATE_FLUID/2,
					   dx[d], dt, d, HLLC);

	    for(int n=0; n<NUM_STATE_FLUID/2; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
		fluxArr(i,j,k,n+NUM_STATE_FLUID/2) = flux_e[n];
	      }	    
	  }
	}
      }
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
#if (AMREX_SPACEDIM >= 2)
      const auto& fluxArrY = fluxes[1].array(mfi);
#endif
  
      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    for(int n=0; n<NUM_STATE_FLUID; n++)
	      //for(int n=0; n<NUM_STATE_FLUID/2; n++)
              {		
		// Conservative update formula
		//arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
		arr(i,j,k,n) = arrOld(i,j,k,n) - (dt / dx[0]) * (fluxArrX(i+1, j, k, n) - fluxArrX(i,j,k,n));
#if (AMREX_SPACEDIM >= 2)
		arr(i,j,k,n) -= (dt / dx[1]) * (fluxArrY(i, j+1, k, n) - fluxArrY(i,j,k,n));
#endif
              }
	    
	    // Initialise to zero
	    arr(i,j,k,DIVB) = 0.0;
	    arr(i,j,k,DIVE) = 1.0/(lambda_d*lambda_d*l_r)*(r_i*arr(i,j,k,RHO_I) + r_e*arr(i,j,k,RHO_E));
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
void CAMReXmp::fluidSolverWENO(MultiFab& S_source, const Real* dx, Real dt)
{

  MultiFab& S_dest = get_new_data(Phi_Type);
  
  // Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes[amrex::SpaceDim];
  for (int j = 0; j < amrex::SpaceDim; j++)
    {
      BoxArray ba = S_dest.boxArray();
      ba.surroundingNodes(j);
      fluxes[j].define(ba, dmap, NUM_STATE, 0);
    }

  MultiFab Slopes(grids, dmap, NUM_STATE_FLUID*5, NUM_GROW);

  const int iDomainOffset = 1;
  const int jDomainOffset = ( amrex::SpaceDim > 1 ? 1 : 0);
  const int kDomainOffset = ( amrex::SpaceDim == 3 ? 1 : 0);

  // Compute slopes for WENO reconstruction
  for (int d = 0; d < amrex::SpaceDim ; d++)   
    {
      
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);
      
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

	  for(int k = lo.z-kDomainOffset; k <= hi.z+kDomainOffset; k++)
	    {
	      for(int j = lo.y-jDomainOffset; j <= hi.y+jDomainOffset; j++)
		{
		  for(int i = lo.x-iDomainOffset; i <= hi.x+iDomainOffset; i++)
		    {
		      /*//for (int n = 0; n<NUM_STATE_FLUID; n++)
		      for(int n=0; n<NUM_STATE_FLUID/2; n++)
			{
			  Vector<Real> dataX = get_data_stencil(arr, i, j, k, 2*iOffset, 2*jOffset, 2*kOffset, n);
			  //Vector<Real> dataX = get_data_stencil(arr, i, j, k, 2, 0, 0, n);
			  std::array<Real, 2> slopesX = WENO3_slope(dataX);		      
			  slopes(i,j,k,n+2*NUM_STATE_FLUID*d) = slopesX[0];
			  //slopes(i,j,k,n) = slopesX[0];
			  slopes(i,j,k,n+NUM_STATE_FLUID+2*NUM_STATE_FLUID*d) = slopesX[1];
			  //slopes(i,j,k,n+NUM_STATE_FLUID) = slopesX[1];
			}
		      */		      
		      //////////////////////////////////////////////////
		      // local characteristic reconstruction
		      WENOcharacteristic(arr,slopes,i,j,k,iOffset,jOffset,kOffset,0,NUM_STATE_FLUID/2,d);
		      WENOcharacteristic(arr,slopes,i,j,k,iOffset,jOffset,kOffset,NUM_STATE_FLUID/2,NUM_STATE_FLUID/2,d);
		      //////////////////////////////////////////////////
		    }
		}
	    }
	}
    }

#if (AMREX_SPACEDIM >= 2)
  for (MFIter mfi(Slopes, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);
      
      const Dim3 hiDomain = ubound(geom.Domain()); 
      
      // Indexable arrays for the data, and the directional flux
      // Based on the vertex-centred definition of the flux array, the
      // data array runs from e.g. [0,N] and the flux array from [0,N+1]
      const auto& arr = S_source.array(mfi);
      const auto& slopes = Slopes.array(mfi);

      for(int k = lo.z-kDomainOffset; k <= hi.z+kDomainOffset; k++)
	{
	  for(int j = lo.y-jDomainOffset; j <= hi.y+jDomainOffset; j++)
	    {
	      for(int i = lo.x-iDomainOffset; i <= hi.x+iDomainOffset; i++)
		{
		  for (int n = 0; n<NUM_STATE_FLUID; n++)
		    //for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataXY = get_data_stencil(arr, i, j, k, 1, 1, 0, n);
		      Real slopesCross = WENO3_slopeCross(dataXY, {slopes(i,j,k,n),slopes(i,j,k,n+NUM_STATE_FLUID),
			slopes(i,j,k,n+2*NUM_STATE_FLUID),slopes(i,j,k,n+3*NUM_STATE_FLUID)});
		      slopes(i,j,k,n+4*NUM_STATE_FLUID) = slopesCross;
		      
		      ////////////////////////////////////////////////////////////////////// 
		      if (!geom.isPeriodic(0) && (i<=1 || i>=hiDomain.x-1))
			{
			  Vector<Real> limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_I);
			  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
			    {
			      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
			      Real slopesX = TVD_slope(dataX,limiterX);
			      slopes(i,j,k,n) = slopesX;
			      slopes(i,j,k,n+NUM_STATE_FLUID) = 0.0;
			      slopes(i,j,k,n+4*NUM_STATE_FLUID) = 0.0;
			    }
			  limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, ENER_E);
			  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
			    {
			      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
			      Real slopesX = TVD_slope(dataX,limiterX);
			      slopes(i,j,k,n) = slopesX;
			      slopes(i,j,k,n+NUM_STATE_FLUID) = 0.0;
			      slopes(i,j,k,n+4*NUM_STATE_FLUID) = 0.0;
			    }
			}
		      if (!geom.isPeriodic(1) && (j<=1 || j>=hiDomain.y-1))
			{
			  Vector<Real> limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_I);
			  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
			    {
			      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
			      Real slopesY = TVD_slope(dataY,limiterY);
			      slopes(i,j,k,n+2*NUM_STATE_FLUID) = slopesY;
			      slopes(i,j,k,n+3*NUM_STATE_FLUID) = 0.0;
			      slopes(i,j,k,n+4*NUM_STATE_FLUID) = 0.0;
			    }
			  limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_E);
			  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
			    {
			      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
			      Real slopesY = TVD_slope(dataY,limiterY);
			      slopes(i,j,k,n+2*NUM_STATE_FLUID) = slopesY;
			      slopes(i,j,k,n+3*NUM_STATE_FLUID) = 0.0;
			      slopes(i,j,k,n+4*NUM_STATE_FLUID) = 0.0;
			    }
			}
		      //////////////////////////////////////////////////////////////////////
		    }
		}
	    }
	}
    } 
#endif
  
  MultiFab positivity(grids, dmap, 2, 1);
  //flattenerAlgorithm(positivity, S_source, Slopes, 0, NUM_STATE_FLUID/2);
  MultiFab positivity_e(grids, dmap, 2, 1);
  //flattenerAlgorithm(positivity_e, S_source, Slopes, NUM_STATE_FLUID/2, NUM_STATE_FLUID/2);
  
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

      const auto& slopes = Slopes.array(mfi);

      const auto& tau = positivity.array(mfi);
      const auto& tau_e = positivity_e.array(mfi);
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = WENO_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					    0, NUM_STATE_FLUID/2, dx[d], dt, d,
					    HLLC);
	    Vector<Real> flux_e = WENO_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					    NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, dx[d], dt, d,
					    HLLC);

	    /*Vector<Real> flux_i = WENO_flux_flat(arr, slopes, tau, i, j, k, iOffset, jOffset, kOffset,
						 0, NUM_STATE_FLUID/2, dx[d], dt, d, HLLC);
	    Vector<Real> flux_e = WENO_flux_flat(arr, slopes, tau_e, i, j, k, iOffset, jOffset, kOffset,
						 NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, dx[d], dt, d, HLLC);
	    */
	    for(int n=0; n<NUM_STATE_FLUID/2; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
		fluxArr(i,j,k,n+NUM_STATE_FLUID/2) = flux_e[n];
	      }	    
	  }
	}
      }
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
#if (AMREX_SPACEDIM >= 2)
      const auto& fluxArrY = fluxes[1].array(mfi);      
#endif

      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {	    
	    for(int n=0; n<NUM_STATE_FLUID; n++)
	      //for(int n=0; n<NUM_STATE_FLUID/2; n++)
              {		
		// Conservative update formula
		//arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
		arr(i,j,k,n) = arrOld(i,j,k,n) - (dt / dx[0]) * (fluxArrX(i+1, j, k, n) - fluxArrX(i,j,k,n));
#if (AMREX_SPACEDIM >= 2)
		arr(i,j,k,n) -= (dt / dx[1]) * (fluxArrY(i, j+1, k, n) - fluxArrY(i,j,k,n));
#endif
              }
	    // Update cell-centred z-components becuause it is 2D code
	    //arr(i,j,k,BZ) = arr(i,j,k,BZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, BZ) - fluxArr(i,j,k,BZ));
	    //arr(i,j,k,EZ) = arr(i,j,k,EZ) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, EZ) - fluxArr(i,j,k,EZ));
	    
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
