
#include <new>
#include <iostream>
#include <iomanip>
#include <regex>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>

using namespace amrex;

amrex::LevelBld* getLevelBld ();

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    Real dRunTime1 = amrex::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

    {
        ParmParse pp;

        max_step  = -1;
        strt_time =  0.0;
        stop_time = -1.0;

        pp.query("max_step",max_step);
        pp.query("strt_time",strt_time);
        pp.query("stop_time",stop_time);
    }

    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0) {
	amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
        Amr amr(getLevelBld());

	amr.init(strt_time,stop_time);
	
	while ( amr.okToContinue() &&
  	       (amr.levelSteps(0) < max_step || max_step < 0) &&
	       (amr.cumTime() < stop_time || stop_time < 0.0) )

	{
	    //
	    // Do a coarse timestep.  Recursively calls timeStep()
	    //
	    amr.coarseTimeStep(stop_time);
	}

	// Write final checkpoint and plotfile
	if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
	    amr.checkPoint();
	}

	if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
	    amr.writePlotFile();
	}

    }

    Real dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    // new code added
    // used to append run times of a simulation in a file
    std::string filename;
    ParmParse ppa("amr");
    ppa.query("plot_file",filename);
    std::regex pattern("plt_");
    filename = std::regex_replace(filename, pattern, "");
    std::regex pattern2("_plt");
    filename = std::regex_replace(filename, pattern2, "_");
#ifdef AMREX_USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (rank==0)
      {
	std::ofstream output(filename + "timing.dat", std::ios_base::app);
	output << dRunTime2 << std::endl;
      }
#elif  _OPENMP
    #pragma omp critical                                                                                                         
      {
	std::cout << "This is OMP" << std::endl;
	std::ofstream output(filename + "_timing.dat", std::ios_base::app);
	output << dRunTime2 << std::endl;
      }
#else
      std::ofstream output(filename + "_timing.dat", std::ios_base::app);
      output << dRunTime2 << std::endl;       
#endif
    
    amrex::Finalize();

    return 0;
}
