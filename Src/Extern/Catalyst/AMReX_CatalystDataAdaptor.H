#ifndef AMReX_CatalystDataAdaptor_H
#define AMReX_CatalystDataAdaptor_H

#include <string>
#include <AMReX_Vector.H>

class vtkCPProcessor;

namespace amrex
{

/// Base class for coupling to Catalyst in situ. implements Initialze and Finalize
/// derived classes implement a CoProcess method. see AmrCatalystDataAdaptor and
/// AmrMeshCatalystDataAdaptor.
class CatalystDataAdaptor
{
public:
    CatalystDataAdaptor();
    virtual ~CatalystDataAdaptor();

    CatalystDataAdaptor(const CatalystDataAdaptor&) = delete;
    void operator=(const CatalystDataAdaptor&) = delete;

    // Runtime control parameters.
    //
    // These may set explicitly or via ParmParse. ParmParse values
    // take precedence over set calls.
    //
    //   catalyst.enabled  -- integer 0 or 1, when 0 Catalyst is disabled
    //   catalyst.frequency -- integer >= 1, how often to CoProcess
    //   catalyst.script -- file path, path to the Catalyst Python script
    //
    void setEnabled(int aenabled) { enabled = aenabled; }
    void setFrequency(int afrequency) { frequency = afrequency; }
    void setScript(const std::string &ascript) { script = ascript; }

    // Initialize Catalyst
    // Use of ParmParse for any unset control parameters
    virtual int Initialize();

    // Finalize Catalyst
    virtual int Finalize();

protected:
    // Do we CoProcess?
    bool doCoProcess();

protected:
    vtkCPProcessor* Processor;
    int counter;
    int enabled;
    int frequency;
    std::string script;
};

}

#endif
