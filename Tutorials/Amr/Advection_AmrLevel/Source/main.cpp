
#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>

#include <AmrLevelAdv.H>

#ifdef BL_USE_CATALYST_INSITU
#ifdef AMREX_PARTICLES
#include <AMReX_AmrParticlesCatalystDataAdaptor.H>
#else
#include <AMReX_AmrCatalystDataAdaptor.H>
#endif
#endif

using namespace amrex;

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
	Amr amr;

	amr.init(strt_time,stop_time);

#if defined(BL_USE_CATALYST_INSITU)
#ifdef AMREX_PARTICLES
	amrex::AmrParticlesCatalystDataAdaptor *insitu_data_adaptor = new amrex::AmrParticlesCatalystDataAdaptor;
#else
    amrex::AmrCatalystDataAdaptor *insitu_data_adaptor = new amrex::AmrCatalystDataAdaptor;
#endif
    insitu_data_adaptor->Initialize();
#endif

	while ( amr.okToContinue() &&
  	       (amr.levelSteps(0) < max_step || max_step < 0) &&
	       (amr.cumTime() < stop_time || stop_time < 0.0) )

	{
	    //
	    // Do a coarse timestep.  Recursively calls timeStep()
	    //
	    amr.coarseTimeStep(stop_time);
#ifdef BL_USE_CATALYST_INSITU
#ifdef AMREX_PARTICLES
	    amrex::ParticleContainer<AMREX_SPACEDIM> * thePC = static_cast<amrex::ParticleContainer<AMREX_SPACEDIM> *> (AmrLevelAdv::theTracerPC());
	    std::vector<std::string> real_comp_names;
	    std::vector<std::string> int_comp_names;
	    insitu_data_adaptor->CoProcess(&amr, *thePC, real_comp_names, int_comp_names);
#else
        insitu_data_adaptor->CoProcess(&amr);
#endif
#endif
	}

	// Write final checkpoint and plotfile
	if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
	    amr.checkPoint();
	}

	if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
	    amr.writePlotFile();
	}

#ifdef BL_USE_CATALYST_INSITU
    insitu_data_adaptor->Finalize();
    delete insitu_data_adaptor;
#endif

    }

    Real dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    amrex::Finalize();

    return 0;
}
