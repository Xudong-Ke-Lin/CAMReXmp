#include <CAMReXmp.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

void CAMReXmp::StrangSecond(const Real* dx, Real dt, Real time)
{

  // copy old data to new data
  {    
    MultiFab& S_new = get_new_data(Phi_Type);
    
    // input states and fill the data
    MultiFab S_input(grids, dmap, NUM_STATE, NUM_GROW);
    FillPatch(*this, S_input, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    MultiFab::Copy(S_new, S_input, 0, 0, NUM_STATE, 0);

#if (AMREX_SPACEDIM >= 2)     
    /*if (MaxwellOrder!=0)
      {		
	MultiFab& S_EM_X_new = get_new_data(EM_X_Type);
	MultiFab& S_EM_Y_new = get_new_data(EM_Y_Type);

	Array<MultiFab,AMREX_SPACEDIM> S_EM_input;
	S_EM_input[0].define(convert(grids,IntVect{AMREX_D_DECL(1,0,0)}), dmap, 6, NUM_GROW);
	FillPatch(*this, S_EM_input[0], NUM_GROW, time, EM_X_Type, 0, 6);
	S_EM_input[1].define(convert(grids,IntVect{AMREX_D_DECL(0,1,0)}), dmap, 6, NUM_GROW);
	FillPatch(*this, S_EM_input[1], NUM_GROW, time, EM_Y_Type, 0, 6);

	MultiFab::Copy(S_EM_X_new, S_EM_input[0], 0, 0, 6, 0);
	MultiFab::Copy(S_EM_Y_new, S_EM_input[1], 0, 0, 6, 0);

	if (MaxwellTimeMethod=="IM" && MaxwellDivMethod=="FDTD")
	  {
	    MultiFab& S_EM_XY_new = get_new_data(EM_XY_Type);
	    MultiFab S_EM_edge_input;
	    S_EM_edge_input.define(convert(grids,IntVect{AMREX_D_DECL(1,1,0)}), dmap, 6, NUM_GROW);
	    FillPatch(*this, S_EM_edge_input, NUM_GROW, time, EM_XY_Type, 0, 6);
	    MultiFab::Copy(S_EM_XY_new, S_EM_edge_input, 0, 0, 6, 0);
	  }
      }
    */
#endif	        
  }

  if (sourceMethod!="no")
    sourceUpdate(0.5*dt, time+dt);

  if (geom.Coord()==1)
    sourceUpdateCyl(dx, 0.5*dt, time+dt);

#if (AMREX_SPACEDIM >= 2)  
  /*if (((sourceMethod=="IM" || sourceMethod=="STIFF" || sourceMethod=="EXACT")
       || (geom.Coord()==1 && MaxwellOrder!=0)) && MaxwellDivMethod!="HDC")
    {
      elecFieldCellAve(time+dt);
    }
  */
#endif

  // two-fluid with GOL solver
  RK2GOL(dx,dt,time+dt);
  // original two-fluid
  //RK2(dx,dt,time+dt,RHO_I,5,RusanovEuler);
  //RK2(dx,dt,time+dt,RHO_E,5,RusanovEuler);
  // Maxwell solver
  if (MaxwellTimeMethod=="EX" && MaxwellDivMethod=="HDC")
    RK2(dx,dt,time+dt,BX,NUM_STATE_MAXWELL,RankineHugoniot);
  if (MaxwellTimeMethod=="IM" && MaxwellDivMethod=="NONE")
    MaxwellSolverCN(dx,dt,time+dt);

#if (AMREX_SPACEDIM >= 2)  
  //if (MaxwellTimeMethod=="IM" && MaxwellDivMethod=="FDTD")
  //MaxwellSolverFDTDCN(dx,dt,time+dt);
  //MaxwellSolverFDTDCNAMReX(dx,dt,time+dt);
#endif
  
  if (sourceMethod!="no") 
    sourceUpdate(0.5*dt, time+dt);

  if (geom.Coord()==1)
    sourceUpdateCyl(dx, 0.5*dt, time+dt);

#if (AMREX_SPACEDIM >= 2)  
  //if (sourceMethod=="IM" && MaxwellDivMethod!="HDC")
  /*if (((sourceMethod=="IM" || sourceMethod=="STIFF" || sourceMethod=="EXACT")
       || (geom.Coord()==1 && MaxwellOrder!=0)) && MaxwellDivMethod!="HDC")
    {
      elecFieldCellAve(time+dt);
      if (projectionStep!= 0 &&
	  parent->levelSteps(0)!=0 &&
	  (parent->levelSteps(0))%projectionStep==0)
	{
	  Projection(dx,time+dt);      
	}
    }
  */
  if (projectionStep!= 0 &&
      parent->levelSteps(0)!=0 &&
      (parent->levelSteps(0))%projectionStep==0)
    {
      Projection(dx,time+dt);      
    }
#endif
}
