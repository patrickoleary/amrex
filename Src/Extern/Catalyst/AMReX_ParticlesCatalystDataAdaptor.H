#ifndef AMReX_ParticlesCatalystDataAdaptor_H
#define AMReX_ParticlesCatalystDataAdaptor_H

#include "AMReX_CatalystDataAdaptor.H"

#include <vector>

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>

class vtkPolyData;
#endif

namespace amrex
{

/// Instrumentation for simulations that use amrex::Amr
/// inherits Intialize and Finalize from amrex::CatalystDataAdaptor
    class ParticlesCatalystDataAdaptor : public CatalystDataAdaptor
    {
    public:
        ParticlesCatalystDataAdaptor() {}
        ~ParticlesCatalystDataAdaptor() {}

        ParticlesCatalystDataAdaptor(const ParticlesCatalystDataAdaptor&) = delete;
        void operator=(const ParticlesCatalystDataAdaptor&) = delete;

#ifdef AMREX_PARTICLES
        template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
        int CoProcess(int step, amrex::Real time,
                      ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer,
                      std::vector<std::string> real_comp_names,
                      std::vector<std::string> int_comp_names);
#endif
    private:
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

        int numberOfParticles;
        vtkPolyData *particleData;
#endif

    };

}

#ifdef AMREX_PARTICLES
#include "AMReX_ParticlesCatalystDataAdaptorI.H"
#endif

#endif
