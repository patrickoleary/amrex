#ifndef AMReX_AmrMeshCatalystDataAdaptor_H
#define AMReX_AmrMeshCatalystDataAdaptor_H

#include "AMReX_CatalystDataAdaptor.H"

#include <string>
#include <vector>

#include <AMReX_Vector.H>

class vtkOverlappingAMR;

namespace amrex
{
    class AmrMesh;
    class AmrMeshDescriptorMap;
    class MultiFab;

/// Instrumentation for simulations that use amrex::AmrMesh/Core
/// inherits intialize and finalize from amrex::CatalystDataAdaptor
class AmrMeshCatalystDataAdaptor : public CatalystDataAdaptor
{
public:
    AmrMeshCatalystDataAdaptor();
    ~AmrMeshCatalystDataAdaptor();

    AmrMeshCatalystDataAdaptor(const AmrMeshCatalystDataAdaptor&) = delete;
    void operator=(const AmrMeshCatalystDataAdaptor&) = delete;

    // invoke CoProcess with data from an AmrMesh instance.
    int CoProcess(unsigned int step, double time, amrex::AmrMesh *mesh,
                  const std::vector<amrex::Vector<amrex::MultiFab>*> &states,
                  const std::vector<std::vector<std::string>> &names);

private:
    amrex::AmrMesh *Mesh;
    std::vector<amrex::Vector<amrex::MultiFab> *> States;
    amrex::AmrMeshDescriptorMap *descriptorMap;
    vtkOverlappingAMR *amrMesh;

    int BuildGrid(int rank);
    int AddGhostCellsArray(int rank);
    int AddArray(int rank, int association, const std::string &arrayName);
    int AddArrays(int rank,
                  const std::vector<amrex::Vector<amrex::MultiFab> *> &states,
                  const std::vector<std::vector<std::string>> &names);
};

}

#endif
