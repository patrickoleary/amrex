#include "AMReX_CatalystDataAdaptor.H"

#include <chrono>

#include <vtkCPProcessor.h>            // ParaView::Catalyst
#include <vtkCPPythonScriptPipeline.h> // ParaView::PythonCatalyst
#include <vtkNew.h>                    // VTK::CommonCore

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

namespace amrex {

CatalystDataAdaptor::CatalystDataAdaptor() :
    Processor(nullptr),
    counter(0), enabled(0), frequency(1)
{
}

CatalystDataAdaptor::~CatalystDataAdaptor()
{
}

int
CatalystDataAdaptor::Initialize()
{
    auto t0 = std::chrono::high_resolution_clock::now();

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
    if (!enabled) return 0;
    bool ret = (frequency > 0) && ((counter % frequency) == 0);
    counter += 1;
    return ret;
}

int
CatalystDataAdaptor::Finalize()
{
    if (!enabled) return 0;
    if (!this->Processor)
        return 0;

    amrex::Print() << "Catalyst Begin finalize..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    this->Processor->Delete();
    this->Processor = nullptr;

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    amrex::Print() << "Catalyst finalize complete (" << dt.count() << " sec)" << std::endl;

    return 0;
}

}
