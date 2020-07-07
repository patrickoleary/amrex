#ifndef AMReX_AmrMeshParticlesCatalystDataAdaptor_H
#define AMReX_AmrMeshParticlesCatalystDataAdaptor_H

#include "AMReX_CatalystDataAdaptor.H"

#include <vector>

#include <AMReX_Vector.H>

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#endif

class vtkOverlappingAMR;
#ifdef AMREX_PARTICLES
class vtkPolyData;
#endif

namespace amrex
{
    class AmrMesh;
    class AmrMeshDescriptorMap;
    class MultiFab;

/// Instrumentation for simulations that use amrex::Amr
/// inherits Intialize and Finalize from amrex::CatalystDataAdaptor
    class AmrMeshParticlesCatalystDataAdaptor : public CatalystDataAdaptor
    {
    public:
        AmrMeshParticlesCatalystDataAdaptor();
        ~AmrMeshParticlesCatalystDataAdaptor();

        AmrMeshParticlesCatalystDataAdaptor(const AmrMeshParticlesCatalystDataAdaptor&) = delete;
        void operator=(const AmrMeshParticlesCatalystDataAdaptor&) = delete;

        int CoProcess(unsigned int step, double time, amrex::AmrMesh *meshDataSource,
                      const std::vector<amrex::Vector<amrex::MultiFab>*> &states,
                      const std::vector<std::vector<std::string>> &names);
#ifdef AMREX_PARTICLES
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int CoProcess(unsigned int step, double time, amrex::AmrMesh *meshDataSource,
                      const std::vector<amrex::Vector<amrex::MultiFab>*> &states,
                      const std::vector<std::vector<std::string>> &names,
                      ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer,
                      std::vector<std::string> real_comp_names,
                      std::vector<std::string> int_comp_names);
#endif
    private:
        int BuildGrid(int rank);
        int AddGhostCellsArray(int rank);
        int AddArray(int rank, int association, const std::string &arrayName);
        int AddArrays(int rank,
                      const std::vector<amrex::Vector<amrex::MultiFab> *> &states,
                      const std::vector<std::vector<std::string>> &names);
#ifdef AMREX_PARTICLES
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int CountParticles(int rank,
                           ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer);
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int BuildParticles(int rank,
                           ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer);
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int AddParticlesArrays(int rank,
                               ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer,
                               std::vector<std::string> real_comp_names,
                               std::vector<std::string> int_comp_names);
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int AddParticlesIDArray(int rank,
                                ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer);
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int AddParticlesCPUArray(int rank,
                                 const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer);
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int AddParticlesAOSRealArray(int rank,
                                     int index,
                                     const std::string &arrayName,
                                     const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer);
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int AddParticlesAOSIntArray(int rank,
                                    int index,
                                    const std::string &arrayName,
                                    const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer);
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int AddParticlesSOARealArray(int rank,
                                     int index,
                                     const std::string &arrayName,
                                     const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer);
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int AddParticlesSOAIntArray(int rank,
                                    int index,
                                    const std::string &arrayName,
                                    const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer);
#endif

        amrex::AmrMesh *Mesh;
        std::vector<amrex::Vector<amrex::MultiFab> *> States;
        amrex::AmrMeshDescriptorMap *descriptorMap;
        vtkOverlappingAMR *amrMesh;

#ifdef AMREX_PARTICLES
        int numberOfParticles;
        vtkPolyData *particleData;
#endif

    };

}

#ifdef AMREX_PARTICLES
#include "AMReX_AmrMeshParticlesCatalystDataAdaptorI.H"
#endif

#endif
