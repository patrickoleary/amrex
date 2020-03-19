#ifndef AMReX_AmrCatalystDataAdaptor_H
#define AMReX_AmrCatalystDataAdaptor_H

#include "AMReX_CatalystDataAdaptor.H"

namespace amrex
{
class Amr;
class DescriptorList;

/// Instrumentation for simulations that use amrex::Amr
/// inherits Intialize and Finalize from amrex::CatalystDataAdaptor
class AmrCatalystDataAdaptor : public CatalystDataAdaptor
{
public:
    AmrCatalystDataAdaptor();
    ~AmrCatalystDataAdaptor();

    AmrCatalystDataAdaptor(const AmrCatalystDataAdaptor&) = delete;
    void operator=(const AmrCatalystDataAdaptor&) = delete;

    // invoke CoProcess with data from an Amr instance
    int CoProcess(amrex::Amr *dataSource);

    int GetGrid();

private:
    struct InternalsType;
    InternalsType *Internals;
    int BuildGrid(int rank);
    int AddGhostCellsArray(int rank);
    int AddArray(int rank, int association, const std::string &arrayName);
    int AddArrays(int rank, const DescriptorList &descriptors);
};

}

#endif
