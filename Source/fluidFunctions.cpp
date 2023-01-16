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
