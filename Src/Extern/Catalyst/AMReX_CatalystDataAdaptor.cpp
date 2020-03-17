#include "AMReX_CatalystDataAdaptor.H"

#include <AMReX_ParmParse.H>

#include <chrono>
#include <timer/Timer.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkNew.h>

namespace amrex {

CatalystDataAdaptor::CatalystDataAdaptor() :
    Processor(nullptr),
    counter(0), enabled(0), frequency(1)
{
    timer::Initialize();
}

CatalystDataAdaptor::~CatalystDataAdaptor()
{
    timer::Finalize();
}

int
CatalystDataAdaptor::Initialize()
{
    auto t0 = std::chrono::high_resolution_clock::now();
    timer::MarkEvent event("CatalystDataAdaptor::Initialize");

    // read config from ParmParse
    ParmParse pp("catalyst");

    pp.query("enabled", enabled);

    if (!enabled)
        return 0;

    pp.query("frequency", frequency);
    pp.query("script", script);

    amrex::Print() << "Catalyst Begin initialize..." << std::endl;

    // Check for invalid values
    if (script.empty())
    {
        amrex::ErrorStream() << "Error: Missing Catalyst Python script." << std::endl;
        return -1;
    }

    if (frequency < 1)
    {
        amrex::ErrorStream() << "Error: Frequency must be greater or equal to 1." << std::endl;
        return -1;
    }

    this->Processor = vtkCPProcessor::New();
    this->Processor->Initialize();
    vtkNew<vtkCPPythonScriptPipeline> pipeline;
    pipeline->Initialize(script.c_str());
    this->Processor->AddPipeline(pipeline);

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    amrex::Print() << "Catalyst initialize complete (" << dt.count() << " sec)" << std::endl;

    return 0;
}

bool
CatalystDataAdaptor::doCoProcess()
{
    bool ret = analysis_adaptor && (frequency > 0) && ((counter % frequency) == 0);
    counter += 1;
    return ret;
}

int
CatalystDataAdaptor::Finalize()
{
    if (!this->Processor)
        return 0;

    amrex::Print() << "Catalyst Begin finalize..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    timer::MarkEvent event("CatalystDataAdaptor::Finalize");
    this->Processor->Delete();
    this->Processor = nullptr;

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    amrex::Print() << "Catalyst finalize complete (" << dt.count() << " sec)" << std::endl;

    return 0;
}

}
