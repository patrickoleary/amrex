#include "AMReX_StateMap.H"

#include <vtkDataObject.h>

#include <AMReX_Print.H>

namespace amrex {

    // --------------------------------------------------------------------------
    int StateMap::GetIndex(const std::string &name, int centering,
        int &desc, int &comp)
    {
        auto cit = this->Map.find(centering);

        if (cit == this->Map.end())
        {
            amrex::Print() << "Error @ StateMap::GetIndex: "
                    << "No " << this->GetAttributesName(centering) << " arrays"
                    << std::endl;
            return -1;
        }

        auto nit = cit->second.find(name);
        if (nit == cit->second.end())
        {
            amrex::Print() << "Error @ StateMap::GetIndex: "
                    << "No array named " << name  << " in "
                    << this->GetAttributesName(centering) << " centered data"
                    << std::endl;
            return -1;
        }

        desc = nit->second.first;
        comp = nit->second.second;

        return 0;
    }

    // --------------------------------------------------------------------------
    int StateMap::GetName(int centering, int id, std::string &name)
    {
        auto cit = this->Map.find(centering);

        if (cit == this->Map.end())
        {
            amrex::Print() << "Error @ StateMap::GetName: "
                    << "No " << this->GetAttributesName(centering) << " arrays"
                    << std::endl;
            return -1;
        }

        if (id >= cit->second.size())
        {
            amrex::Print() << "Error @ StateMap::GetName: "
                    << "Array index " << id << " out of bounds " << cit->second.size()
                    << std::endl;
            return -1;
        }

        auto nit = cit->second.begin();
        while (id--)
            ++nit;

        name = nit->first;

        return 0;
    }

    //-----------------------------------------------------------------------------
    const char *StateMap::GetAttributesName(int association)
    {
        switch (association)
        {
            case vtkDataObject::POINT:
                return "point";
                break;
            case vtkDataObject::CELL:
                return "cell";
                break;
        }
        amrex::Print() << "Warning @ StateMap::GetAttributesName: "
                << "Invalid association"
                << std::endl;
        return "";
    }

}
