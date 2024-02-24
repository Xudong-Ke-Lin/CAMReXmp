#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::MUSCLHancokPlasmaSolverTVD(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  
  MultiFab Slopes(grids, dmap, NUM_STATE*2, NUM_GROW);

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
		  for (int n = BX; n<BX+NUM_STATE_MAXWELL; n++)
		    {
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      Real slopesX = TVD_slope(dataX,limiterX);
		      slopes(i,j,k,n) = slopesX;		      
		    }		  
#if (AMREX_SPACEDIM >= 2)
		  Vector<Real> limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_I);
		  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE) = slopesY;		      
		    }
		  limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_E);
		  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE) = slopesY;		      
		    }
		  for (int n = BX; n<BX+NUM_STATE_MAXWELL; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, n); 
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE) = slopesY;		      
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
	    Vector<Real> flux_em = MUSCL_Hancock_TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
							  BX, NUM_STATE_MAXWELL, dx[d], dt, d, MaxwellFlux, RankineHugoniot);
	    
	    //Vector<Real> flux = fluid_flux_HLLC(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);

	    for(int n=0; n<NUM_STATE_FLUID/2; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
		fluxArr(i,j,k,n+NUM_STATE_FLUID/2) = flux_e[n];
	      }
	    for(int n=BX; n<BX+NUM_STATE_MAXWELL; n++)
	      {		
		fluxArr(i,j,k,n) = flux_em[n-BX];
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
      	    // Update fluid variables
      	    for(int n=0; n<NUM_STATE; n++)
	      //for(int n=0; n<NUM_STATE_FLUID; n++)
              {
      		// Conservative update formula
      		arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
              }
      	  }
      	}
      }
    }
      
    // We need to compute boundary conditions again after each update
    S_source.FillBoundary(geom.periodicity());
     
    // added by 2020D 
    // Fill non-periodic physical boundaries
    FillDomainBoundary(S_source, geom, bc);
    
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
      //const auto& fluxArrX = fluxes[0].array(mfi);
#if (AMREX_SPACEDIM >= 2)
      //const auto& fluxArrY = fluxes[1].array(mfi);      
#endif

      for(int k = lo.z; k <= hi.z; k++)
      {
	for(int j = lo.y; j <= hi.y; j++)
	{
	  for(int i = lo.x; i <= hi.x; i++)
	  {
	    for(int n=0; n<NUM_STATE; n++)
	    //for(int n=0; n<NUM_STATE_FLUID; n++)
	      //for(int n=0; n<NUM_STATE_FLUID/2; n++)
              {		
		// Conservative update formula
		//arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
		arr(i,j,k,n) = arrOld(i,j,k,n);// - (dt / dx[0]) * (fluxArrX(i+1, j, k, n) - fluxArrX(i,j,k,n));
#if (AMREX_SPACEDIM >= 2)
		//arr(i,j,k,n) -= (dt / dx[1]) * (fluxArrY(i, j+1, k, n) - fluxArrY(i,j,k,n));
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
void CAMReXmp::plasmaSolverTVD(MultiFab& S_dest, MultiFab& S_source, MultiFab (&fluxes)[AMREX_SPACEDIM], const Real* dx, Real dt)
{
  
  MultiFab Slopes(grids, dmap, NUM_STATE*2, NUM_GROW);

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
		  for (int n = BX; n<BX+NUM_STATE_MAXWELL; n++)
		    {
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      Real slopesX = TVD_slope(dataX,limiterX);
		      slopes(i,j,k,n) = slopesX;		      
		    }		  
#if (AMREX_SPACEDIM >= 2)
		  Vector<Real> limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_I);
		  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE) = slopesY;		      
		    }
		  limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_E);
		  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE) = slopesY;		      
		    }
		  for (int n = BX; n<BX+NUM_STATE_MAXWELL; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, n); 
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE) = slopesY;		      
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
      //const auto& arrDest = S_dest.array(mfi);
      const auto& fluxArr = fluxes[d].array(mfi);      

      const auto& slopes = Slopes.array(mfi);
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					   0, NUM_STATE_FLUID/2, d*NUM_STATE,
					   dx[d], dt, d, HLLC);	    
	    Vector<Real> flux_e = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					   NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, d*NUM_STATE+NUM_STATE_FLUID/2,
					   dx[d], dt, d, HLLC);	    
	    Vector<Real> flux_em = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					    BX, NUM_STATE_MAXWELL, d*NUM_STATE+BX,
					    dx[d], dt, d, RankineHugoniot);
	    
	    //Vector<Real> flux = fluid_flux_HLLC(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);

	    for(int n=0; n<NUM_STATE_FLUID/2; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
		fluxArr(i,j,k,n+NUM_STATE_FLUID/2) = flux_e[n];
	      }
	    for(int n=BX; n<BX+NUM_STATE_MAXWELL; n++)
	      {		
		fluxArr(i,j,k,n) = flux_em[n-BX];
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
      	    for(int n=0; n<NUM_STATE; n++)
	      //for(int n=0; n<NUM_STATE_FLUID; n++)
              {
      		// Conservative update formula
      		arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
              }
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
	    for(int n=0; n<NUM_STATE; n++)
	    //for(int n=0; n<NUM_STATE_FLUID; n++)
	      //for(int n=0; n<NUM_STATE_FLUID/2; n++)
              {		
		// Conservative update formula
		//arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset, n) - fluxArr(i,j,k,n));
		arr(i,j,k,n) = arrOld(i,j,k,n) - (dt / dx[0]) * (fluxArrX(i+1, j, k, n) - fluxArrX(i,j,k,n));
#if (AMREX_SPACEDIM >= 2)
		arr(i,j,k,n) -= (dt / dx[1]) * (fluxArrY(i, j+1, k, n) - fluxArrY(i,j,k,n));
#endif
		/*if (sourceMethod!="IM")
		  {
		    arr(i,j,k,EX) -= dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arrOld(i,j,k,1) + r_e*arrOld(i,j,k,MOMX_E));
		    arr(i,j,k,EY) -= dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arrOld(i,j,k,2) + r_e*arrOld(i,j,k,MOMY_E));
		    arr(i,j,k,EZ) -= dt*1.0/(lambda_d*lambda_d*l_r)*(r_i*arrOld(i,j,k,MOMZ_I) + r_e*arrOld(i,j,k,MOMZ_E));
		    }*/
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
void CAMReXmp::plasmaSolverTVD(MultiFab& S_source, const Real* dx, Real dt)
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
  
  MultiFab Slopes(grids, dmap, NUM_STATE*2, NUM_GROW);

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
		  for (int n = BX; n<BX+NUM_STATE_MAXWELL; n++)
		    {
		      Vector<Real> dataX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      limiterX = get_data_stencil(arr, i, j, k, 1, 0, 0, n);
		      Real slopesX = TVD_slope(dataX,limiterX);
		      slopes(i,j,k,n) = slopesX;		      
		    }		  
#if (AMREX_SPACEDIM >= 2)
		  Vector<Real> limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_I);
		  for (int n = 0; n<NUM_STATE_FLUID/2; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE) = slopesY;		      
		    }
		  limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, ENER_E);
		  for (int n = NUM_STATE_FLUID/2; n<NUM_STATE_FLUID; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE) = slopesY;		      
		    }
		  for (int n = BX; n<BX+NUM_STATE_MAXWELL; n++)
		    {
		      Vector<Real> dataY = get_data_stencil(arr, i, j, k, 0, 1, 0, n);
		      limiterY = get_data_stencil(arr, i, j, k, 0, 1, 0, n); 
		      Real slopesY = TVD_slope(dataY,limiterY);
		      slopes(i,j,k,n+NUM_STATE) = slopesY;		      
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
      //const auto& arrDest = S_dest.array(mfi);
      const auto& fluxArr = fluxes[d].array(mfi);      

      const auto& slopes = Slopes.array(mfi);
      
      for(int k = lo.z; k <= hi.z+kOffset; k++)
      {
	for(int j = lo.y; j <= hi.y+jOffset; j++)
	{
	  for(int i = lo.x; i <= hi.x+iOffset; i++)
	  {
	    Vector<Real> flux_i = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					   0, NUM_STATE_FLUID/2, d*NUM_STATE,
					   dx[d], dt, d, HLLC);	    
	    Vector<Real> flux_e = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					   NUM_STATE_FLUID/2, NUM_STATE_FLUID/2, d*NUM_STATE+NUM_STATE_FLUID/2,
					   dx[d], dt, d, HLLC);	    
	    Vector<Real> flux_em = TVD_flux(arr, slopes, i, j, k, iOffset, jOffset, kOffset,
					    BX, NUM_STATE_MAXWELL, d*NUM_STATE+BX,
					    dx[d], dt, d, RankineHugoniot);
	    
	    //Vector<Real> flux = fluid_flux_HLLC(arr, i, j, k, iOffset, jOffset, kOffset, dx[d], dt, d);

	    for(int n=0; n<NUM_STATE_FLUID/2; n++)
	      {		
		fluxArr(i,j,k,n) = flux_i[n];
		fluxArr(i,j,k,n+NUM_STATE_FLUID/2) = flux_e[n];
	      }
	    for(int n=BX; n<BX+NUM_STATE_MAXWELL; n++)
	      {		
		fluxArr(i,j,k,n) = flux_em[n-BX];
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
	    for(int n=0; n<NUM_STATE; n++)
	    //for(int n=0; n<NUM_STATE_FLUID; n++)
	      //for(int n=0; n<NUM_STATE_FLUID/2; n++)
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
