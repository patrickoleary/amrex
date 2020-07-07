#include "AMReX_AmrDescriptorMap.H"

#include <AMReX_IndexType.H>

#include <vtkDataObject.h>          // VTK::CommonDataModel

namespace amrex {

// --------------------------------------------------------------------------
    int AmrDescriptorMap::Initialize(const DescriptorList &descriptors)
    {
        int ndesc = descriptors.size();
        for (int i = 0; i < ndesc; ++i)
        {
            const StateDescriptor &desc = descriptors[i];
            int ncomp = desc.nComp();
            IndexType itype = desc.getType();
            for (int j = 0; j < ncomp; ++j)
            {
                std::string arrayName = desc.name(j);
                if (itype.cellCentered())
                {
                    this->Map[vtkDataObject::CELL][arrayName] = std::make_pair(i,j);
                }
                else if (itype.nodeCentered())
                {
                    this->Map[vtkDataObject::POINT][arrayName] = std::make_pair(i,j);
                }
            }
        }
        return 0;
    }
}
