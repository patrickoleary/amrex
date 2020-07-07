#ifndef AMREX_AMRDESCRIPTORMAP_H
#define AMREX_AMRDESCRIPTORMAP_H

#include <AMReX_StateDescriptor.H>
#include <AMReX_StateMap.H>

namespace amrex {
    class AmrDescriptorMap : public amrex::StateMap {
    public:
        int Initialize(const DescriptorList &descriptors);
    };
}

#endif
