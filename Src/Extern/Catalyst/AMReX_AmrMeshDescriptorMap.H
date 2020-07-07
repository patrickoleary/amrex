#ifndef AMREX_AMRMESHDESCRIPTORMAP_H
#define AMREX_AMRMESHDESCRIPTORMAP_H

#include <AMReX_StateDescriptor.H>
#include <AMReX_StateMap.H>
#include <AMReX_Vector.H>

namespace amrex {
    class MultiFab;

    class AmrMeshDescriptorMap : public amrex::StateMap {
    public:
        int Initialize(const std::vector<amrex::Vector<amrex::MultiFab> *> &states,
                       const std::vector<std::vector<std::string>> &names);
    };
}

#endif
