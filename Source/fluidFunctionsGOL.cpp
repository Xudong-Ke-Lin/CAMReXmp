#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::fluidSolverMonotoneGOL(MultiFab& S_source, const Real* dx, Real dt)
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
	    // combined fluid-GOL system
	    Vector<Real> flux_i = monotone_fluxGOL(arr, i, j, k, iOffset, jOffset, kOffset,
						   0, NUM_STATE_FLUID,d);
	    // FORCE solver
	    /*Vector<Real> flux_i = monotone_fluxGOL(arr, i, j, k, iOffset, jOffset, kOffset,
						   0, NUM_STATE_FLUID,d,dx[d],dt);
	    */
	    for(int n=0; n<NUM_STATE_FLUID; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
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
	    arr(i,j,k,DIVE) = 0.0;
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
void CAMReXmp::fluidSolverTVDGOL(MultiFab& S_source, const Real* dx, Real dt)
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

  const int iDomainOffset = 1;
  const int jDomainOffset = ( amrex::SpaceDim > 1 ? 1 : 0);
  const int kDomainOffset = ( amrex::SpaceDim == 3 ? 1 : 0);
  
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

      for(int k = lo.z-kDomainOffset; k <= hi.z+kDomainOffset; k++)
	{
	  for(int j = lo.y-jDomainOffset; j <= hi.y+jDomainOffset; j++)
	    {
	      for(int i = lo.x-iDomainOffset; i <= hi.x+iDomainOffset; i++)
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
	    // combined fluid-GOL system
	    Vector<Real> flux_i = TVD_fluxGOL(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					      0, NUM_STATE_FLUID, d*NUM_STATE_FLUID,dx[d],dt,d);

	    for(int n=0; n<NUM_STATE_FLUID; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
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
	    arr(i,j,k,DIVE) = 0.0;
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
void CAMReXmp::SLICfluidSolverTVDGOL(MultiFab& S_source, const Real* dx, Real dt)
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

  const int iDomainOffset = 1;
  const int jDomainOffset = ( amrex::SpaceDim > 1 ? 1 : 0);
  const int kDomainOffset = ( amrex::SpaceDim == 3 ? 1 : 0);

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

      for(int k = lo.z-kDomainOffset; k <= hi.z+kDomainOffset; k++)
	{
	  for(int j = lo.y-jDomainOffset; j <= hi.y+jDomainOffset; j++)
	    {
	      for(int i = lo.x-iDomainOffset; i <= hi.x+iDomainOffset; i++)
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
	    // combined fluid-GOL system
	    Vector<Real> flux_i = SLIC_TVD_fluxGOL(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
						   0, NUM_STATE_FLUID, d*NUM_STATE_FLUID,dx[d],dt,d);

	    for(int n=0; n<NUM_STATE_FLUID; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
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
	    arr(i,j,k,DIVE) = 0.0;
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// semi-implicit method for advection sub-system
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Implicit Maxwell Solver
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>

void CAMReXmp::fluidSolverGOLelecPres(MultiFab& S_source, const Real* dx, Real dt, Real time)
{

  MultiFab& S_dest = get_new_data(Phi_Type);

  MultiFab S_pressure(grids, dmap, 1, NUM_GROW);
  for(MFIter mfi(S_pressure, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arrP = S_pressure.array(mfi);
      const auto& arr = S_source.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{
		  Vector<Real> u_i = get_data_zone(arr,i,j,k,0,10);
		  Vector<Real> u_e = get_electron_var(u_i);
		  // pressure
		  arrP(i,j,k,0) = get_pressure(u_e);
		}
	    }
	}	  
    }

  // We need to compute boundary conditions again after each update
  S_pressure.FillBoundary(geom.periodicity());

  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_pressure, geom, {bc[ENER_E]});

  // Picard iteration start
  int iterFin = 2;
  MultiFab S_tmp(grids, dmap, NUM_STATE, NUM_GROW);  
  for (int iter = 0; iter<iterFin; iter++)
  {
  FillPatch(*this, S_tmp, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

  // b coefficients for linear solver
  // these are also the face-centered enthalpies
  std::array<MultiFab, BL_SPACEDIM> bcoeffs;
  for(int n = 0; n < BL_SPACEDIM; n++)
    {
      const BoxArray& ba = convert(S_source.boxArray(), IntVect::TheDimensionVector(n));
      bcoeffs[n].define(ba, S_source.DistributionMap(), 1, 0);
    }  
  for(int d = 0; d < BL_SPACEDIM; d++)
    {
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      for(MFIter mfi(bcoeffs[d], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();

	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const auto& arrB = bcoeffs[d].array(mfi);
	  const auto& arrP = S_pressure.array(mfi);
	  const auto& arr = S_tmp.array(mfi);

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      Vector<Real> u_i = get_data_zone(arr,i-iOffset,j-jOffset,k-kOffset,0,10);
		      Vector<Real> u_iPlus1 = get_data_zone(arr,i,j,k,0,10);

		      Vector<Real> u_ei = get_electron_var(u_i);
		      Vector<Real> u_eiPlus1 = get_electron_var(u_iPlus1);
		      
		      Real rho_ei = u_ei[RHO_I];
		      Real rho_eiPlus1 = u_eiPlus1[RHO_I];

		      // pressure
		      Real p_ei = arrP(i-iOffset,j-jOffset,k-kOffset,0);
		      Real p_eiPlus1 = arrP(i,j,k,0);

		      // internal energy
		      Real e_ei = p_ei/((Gamma-1));
		      Real e_eiPlus1 = p_eiPlus1/((Gamma-1));

		      // enthalpy
		      Real h_ei = e_ei+p_ei;
		      Real h_eiPlus1 = e_eiPlus1+p_eiPlus1;

		      // momentum
		      Real mom_ei = get_magnitude(u_ei[MOMX_I],u_ei[MOMY_I],u_ei[MOMZ_I]);
		      Real mom_eiPlus1 = get_magnitude(u_eiPlus1[MOMX_I],u_eiPlus1[MOMY_I],u_eiPlus1[MOMZ_I]);

		      if (std::abs(mom_ei+mom_eiPlus1)<1e-12)
			arrB(i,j,k,0) = 0.5*(h_ei/rho_ei + h_eiPlus1/rho_eiPlus1);
		      else
			arrB(i,j,k,0) = (h_ei*mom_ei/rho_ei + h_eiPlus1*mom_eiPlus1/rho_eiPlus1)/(mom_ei+mom_eiPlus1);
		    }
		}
	    }	  
	}
    }

  MultiFab Rhs(grids, dmap, 1, NUM_GROW);

  for (int d = 0; d < amrex::SpaceDim ; d++)   
  {

    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);

    for(MFIter mfi(Rhs, true); mfi.isValid(); ++mfi)
      {
	const Box& bx = mfi.tilebox();

	const Dim3 lo = lbound(bx);
	const Dim3 hi = ubound(bx);

	// old and new data
	const auto& arrOld = S_source.array(mfi);
	const auto& arrNew = S_tmp.array(mfi);
	// enthalpies
	const auto& arrH = bcoeffs[d].array(mfi);
	const auto& rhs = Rhs.array(mfi);

	for(int k = lo.z; k <= hi.z; k++)
	  {
	    for(int j = lo.y; j <= hi.y; j++)
	      {
		for(int i = lo.x; i <= hi.x; i++)
		  {
		    Real h_eiPlusHalf = arrH(i+iOffset,j+jOffset,k+kOffset,0);
		    Real h_eiMinusHalf = arrH(i,j,k,0);

		    if (d==0)
		      {
			Vector<Real> u_i = get_data_zone(arrNew,i-iOffset,j-jOffset,k-kOffset,0,10);

			Vector<Real> u_ei = get_electron_var(u_i);
			
			// kinetic energy
			Real kin_i = 0.5*get_magnitude_squared(u_ei[MOMX_I],u_ei[MOMY_I],u_ei[MOMZ_I])/u_ei[RHO_I];

			rhs(i,j,k,0) = arrOld(i,j,k,ENER_E) - kin_i;
		      }

		    // quasineutrality n_i=n_e, rho_e=rho_i*(m_e/m_i) and m_i=1
		    Real m_e = 1.0/m;
		    Real q_i = r_i;

		    rhs(i,j,k,0) -= dt/(2.0*dx[d])*m_e*(h_eiPlusHalf*(arrOld(i+iOffset,j+jOffset,k+kOffset,MOMX_I+d)-arrOld(i+iOffset,j+jOffset,k+kOffset,MOMX_E+d)/q_i)
							+ (h_eiPlusHalf-h_eiMinusHalf)*(arrOld(i,j,k,MOMX_I+d)-arrOld(i,j,k,MOMX_E+d)/q_i)
							- h_eiMinusHalf*(arrOld(i-iOffset,j-jOffset,k-kOffset,MOMX_I+d)-arrOld(i-iOffset,j-jOffset,k-kOffset,MOMX_E+d)/q_i));
		  }
	      }
	  }
      }
  }

  // For MLMG solver
  int verbose = 2;
  int bottom_verbose = 0;
  int max_iter = 100;
  //int max_fmg_iter = 0;
  int linop_maxorder = 2;
  bool agglomeration = true;
  bool consolidation = true;
  bool semicoarsening = false;
  int max_coarsening_level = 30;
  int max_semicoarsening_level = 0;

  LPInfo info;
  info.setAgglomeration(agglomeration);
  info.setConsolidation(consolidation);
  info.setSemicoarsening(semicoarsening);
  info.setMaxCoarseningLevel(max_coarsening_level);
  info.setMaxSemicoarseningLevel(max_semicoarsening_level);

  const auto tol_rel = Real(1.e-10);
  const auto tol_abs = Real(0.0);

  MLABecLaplacian mlabec({geom}, {grids}, {dmap}, info);
  mlabec.setMaxOrder(linop_maxorder);

  // Set boundary conditions for MLABecLaplacian
  std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_lobc;
  std::array<LinOpBCType, AMREX_SPACEDIM> mlmg_hibc;  
  setDomainBC(mlmg_lobc, mlmg_hibc, ENER_E);
  mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);

  // Set boundary conditions for the current patch 
  mlabec.setLevelBC(0,&S_pressure);

  Real ascalar = 1.0/(Gamma-1.0);
  Real bscalar = dt*dt;
  mlabec.setScalars(ascalar, bscalar);

  mlabec.setACoeffs(0, 1.0);
  mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoeffs));
  MLMG mlmg(mlabec);

  mlmg.setMaxIter(max_iter);
  mlmg.setMaxFmgIter(max_fmg_iter);
  mlmg.setVerbose(verbose);
  mlmg.setBottomVerbose(bottom_verbose);  

  mlmg.solve({&S_pressure}, {&Rhs}, tol_rel, tol_abs);

  // We need to compute boundary conditions again after each update
  S_pressure.FillBoundary(geom.periodicity());

  // added by 2020D 
  // Fill non-periodic physical boundaries                          
  FillDomainBoundary(S_pressure, geom, {bc[ENER_E]});

  for (int d = 0; d < amrex::SpaceDim ; d++)
  {
    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);

    for(MFIter mfi(S_tmp, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arrP = S_pressure.array(mfi);
      const auto& arr = S_tmp.array(mfi);
      const auto& arrOld = S_source.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{

		  // new momentum
		  arr(i,j,k,MOMX_I+d) = arrOld(i,j,k,MOMX_I+d)
		    - 0.5*dt/dx[d]*(arrP(i+iOffset,j+jOffset,k+kOffset,0)
				    -arrP(i-iOffset,j-jOffset,k-kOffset,0));

		  // quasineutrality n_i=n_e, rho_e=rho_i*(m_e/m_i) and m_i=1
		  Real m_e = 1.0/m;
		  Real q_i = r_i;

		  // new current
		  arr(i,j,k,MOMX_E+d) = arrOld(i,j,k,MOMX_E+d)
		    +(q_i/m_e)*0.5*dt/dx[d]*(arrP(i+iOffset,j+jOffset,k+kOffset,0)
					     -arrP(i-iOffset,j-jOffset,k-kOffset,0));

		}
	    }
	}	  
    }
  }
  MultiFab::Copy(S_dest, S_tmp, 0, 0, 10, 0);  
  }
  // fill data with new updated momentum
  FillPatch(*this, S_tmp, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  std::array<MultiFab, BL_SPACEDIM> S_enthalpy;
  for(int n = 0; n < BL_SPACEDIM; n++)
    {
      const BoxArray& ba = convert(S_source.boxArray(), IntVect::TheDimensionVector(n));
      S_enthalpy[n].define(ba, S_source.DistributionMap(), 1, 0);
    }  
  for(int d = 0; d < BL_SPACEDIM; d++)
    {
      const int iOffset = ( d == 0 ? 1 : 0);
      const int jOffset = ( d == 1 ? 1 : 0);
      const int kOffset = ( d == 2 ? 1 : 0);

      for(MFIter mfi(S_enthalpy[d], true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();

	  const Dim3 lo = lbound(bx);
	  const Dim3 hi = ubound(bx);

	  const auto& arrB = S_enthalpy[d].array(mfi);
	  const auto& arrP = S_pressure.array(mfi);
	  const auto& arr = S_tmp.array(mfi);

	  for(int k = lo.z; k <= hi.z; k++)
	    {
	      for(int j = lo.y; j <= hi.y; j++)
		{
		  for(int i = lo.x; i <= hi.x; i++)
		    {
		      Vector<Real> u_i = get_data_zone(arr,i-iOffset,j-jOffset,k-kOffset,0,10);
		      Vector<Real> u_iPlus1 = get_data_zone(arr,i,j,k,0,10);

		      Vector<Real> u_ei = get_electron_var(u_i);
		      Vector<Real> u_eiPlus1 = get_electron_var(u_iPlus1);
		      
		      Real rho_ei = u_ei[RHO_I];
		      Real rho_eiPlus1 = u_eiPlus1[RHO_I];

		      // pressure
		      Real p_ei = arrP(i-iOffset,j-jOffset,k-kOffset,0);
		      Real p_eiPlus1 = arrP(i,j,k,0);

		      // internal energy
		      Real e_ei = p_ei/((Gamma-1));
		      Real e_eiPlus1 = p_eiPlus1/((Gamma-1));

		      // enthalpy
		      Real h_ei = e_ei+p_ei;
		      Real h_eiPlus1 = e_eiPlus1+p_eiPlus1;

		      // momentum
		      Real mom_ei = get_magnitude(u_ei[MOMX_I],u_ei[MOMY_I],u_ei[MOMZ_I]);
		      Real mom_eiPlus1 = get_magnitude(u_eiPlus1[MOMX_I],u_eiPlus1[MOMY_I],u_eiPlus1[MOMZ_I]);

		      if (std::abs(mom_ei+mom_eiPlus1)<1e-12)
			arrB(i,j,k,0) = 0.5*(h_ei/rho_ei + h_eiPlus1/rho_eiPlus1);
		      else
			arrB(i,j,k,0) = (h_ei*mom_ei/rho_ei + h_eiPlus1*mom_eiPlus1/rho_eiPlus1)/(mom_ei+mom_eiPlus1);

		    }
		}
	    }	  
	}
    }
  for (int d = 0; d < amrex::SpaceDim ; d++)
  {
    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);

    for(MFIter mfi(S_tmp, true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.tilebox();

      const Dim3 lo = lbound(bx);
      const Dim3 hi = ubound(bx);

      const auto& arrP = S_pressure.array(mfi);
      const auto& arrH = S_enthalpy[d].array(mfi);
      const auto& arr = S_tmp.array(mfi);
      const auto& arrOld = S_source.array(mfi);

      for(int k = lo.z; k <= hi.z; k++)
	{
	  for(int j = lo.y; j <= hi.y; j++)
	    {
	      for(int i = lo.x; i <= hi.x; i++)
		{

		  Vector<Real> u_iMinus1 = get_data_zone(arr,i-iOffset,j-jOffset,k-kOffset,0,10);
		  Vector<Real> u_i = get_data_zone(arr,i,j,k,0,10);
		  Vector<Real> u_iPlus1 = get_data_zone(arr,i+iOffset,j+jOffset,k+kOffset,0,10);

		  Vector<Real> u_eiMinus1 = get_electron_var(u_iMinus1);
		  Vector<Real> u_ei = get_electron_var(u_i);
		  Vector<Real> u_eiPlus1 = get_electron_var(u_iPlus1);

		  Real momX_eiMinus1 = u_eiMinus1[MOMX_I+d];
		  Real momX_ei = u_ei[MOMX_I+d];
		  Real momX_eiPlus1 = u_eiPlus1[MOMX_I+d];
		  if (d==0)
		    {
		      arr(i,j,k,ENER_I) = arrOld(i,j,k,ENER_I);
		      arr(i,j,k,ENER_E) = arrOld(i,j,k,ENER_E);
		    }
		  arr(i,j,k,ENER_I) -= dt/(2.0*dx[d])*(arrH(i+iOffset,j+jOffset,k+kOffset,0)*(momX_eiPlus1+momX_ei)
						       -arrH(i,j,k,0)*(momX_ei+momX_eiMinus1));
		  arr(i,j,k,ENER_E) -= dt/(2.0*dx[d])*(arrH(i+iOffset,j+jOffset,k+kOffset,0)*(momX_eiPlus1+momX_ei)
						       -arrH(i,j,k,0)*(momX_ei+momX_eiMinus1));
		  
		}
	    }
	}	  
    }
  }
  MultiFab::Copy(S_dest, S_tmp, ENER_I, ENER_I, 1, 0);
  MultiFab::Copy(S_dest, S_tmp, ENER_E, ENER_E, 1, 0);
}
