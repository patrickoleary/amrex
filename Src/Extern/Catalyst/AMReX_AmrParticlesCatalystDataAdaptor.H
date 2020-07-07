#ifndef AMReX_AmrParticlesCatalystDataAdaptor_H
#define AMReX_AmrParticlesCatalystDataAdaptor_H

#include "AMReX_CatalystDataAdaptor.H"

#include <vector>

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#endif

class vtkOverlappingAMR;
#ifdef AMREX_PARTICLES
class vtkPolyData;
#endif

namespace amrex
{
    class Amr;
    class AmrDescriptorMap;
    class DescriptorList;

/// Instrumentation for simulations that use amrex::Amr
/// inherits Intialize and Finalize from amrex::CatalystDataAdaptor
    class AmrParticlesCatalystDataAdaptor : public CatalystDataAdaptor
    {
    public:
        AmrParticlesCatalystDataAdaptor();
        ~AmrParticlesCatalystDataAdaptor();

        AmrParticlesCatalystDataAdaptor(const AmrParticlesCatalystDataAdaptor&) = delete;
        void operator=(const AmrParticlesCatalystDataAdaptor&) = delete;

        int CoProcess(Amr *meshDataSource);
#ifdef AMREX_PARTICLES
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int CoProcess(Amr *meshDataSource,
                      ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer,
                      std::vector<std::string> real_comp_names,
                      std::vector<std::string> int_comp_names);
#endif
    private:
        int BuildGrid(int rank);
        int AddGhostCellsArray(int rank);
        int AddArrays(int rank, const DescriptorList &descriptors);
        int AddArray(int rank, int association, const std::string &arrayName);
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

        amrex::Amr *amr;
        amrex::AmrDescriptorMap *descriptorMap;
        vtkOverlappingAMR *amrMesh;

#ifdef AMREX_PARTICLES
        int numberOfParticles;
        vtkPolyData *particleData;
#endif

    };

}

#ifdef AMREX_PARTICLES
#include "AMReX_AmrParticlesCatalystDataAdaptorI.H"
#endif

#endif
