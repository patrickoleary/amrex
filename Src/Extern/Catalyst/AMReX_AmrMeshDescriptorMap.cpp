#include "AMReX_AmrMeshDescriptorMap.H"

#include <AMReX_MultiFab.H>

#include <vtkDataObject.h>          // VTK::CommonDataModel

namespace amrex {

// --------------------------------------------------------------------------
    int AmrMeshDescriptorMap::Initialize(const std::vector<amrex::Vector<amrex::MultiFab> *> &states,
                                         const std::vector<std::vector<std::string>> &names)
    {
        int nStates = states.size();
        for (int i = 0; i < nStates; ++i)
        {
            amrex::MultiFab& state = states[i]->at(0);
            int nComp = state.nComp();

            for (int j = 0; j < nComp; ++j)
            {
                const std::string &arrayName = names[i][j];

                if (state.is_cell_centered())
                {
                    this->Map[vtkDataObject::CELL][arrayName] = std::make_pair(i,j);
                }
                else if (state.is_nodal())
                {
                    this->Map[vtkDataObject::POINT][arrayName] = std::make_pair(i,j);
                }
            }
        }
        return 0;
    }

}
