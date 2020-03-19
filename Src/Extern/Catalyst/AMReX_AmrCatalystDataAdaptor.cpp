#include "AMReX_AmrCatalystDataAdaptor.H"

#include <chrono>
#include <map>
#include <utility>

#include <vtkAMRBox.h>                 // VTK::CommonDataModel
#include <vtkCellData.h>               // VTK::CommonDataModel
#include <vtkCPDataDescription.h>      // ParaView::Catalyst
#include <vtkCPInputDataDescription.h> // ParaView::Catalyst
#include <vtkDataSetAttributes.h>      // VTK::CommonDataModel
#include <vtkDoubleArray.h>            // VTK::CommonCore
#include <vtkFloatArray.h>             // VTK::CommonCore
#include <vtkMultiProcessController.h> // VTK::ParallelCore
#include <vtkObject.h>                 // VTK::CommonCore
#include <vtkObjectFactory.h>          // VTK::CommonCore
#include <vtkOverlappingAMR.h>         // VTK::CommonDataModel
#include <vtkPointData.h>              // VTK::CommonDataModel
#include <vtkCPProcessor.h>            // ParaView::Catalyst
#include <vtkUniformGrid.h>            // VTK::CommonDataModel
#include <vtkUnsignedCharArray.h>      // VTK::CommonCore

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_Box.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_StateData.H>
#include <AMReX_MFIter.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IndexType.H>
#include <AMReX_VTKGhostUtils.H>
#include <AMReX_StateMap.H>
#include <AMReX_Print.H>

// return the number of levels currently in use
static
unsigned int numActiveLevels(
        amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels)
{
    unsigned int nLevels = levels.size();
    for (int i = 0; i < nLevels; ++i)
    {
        if (!levels[i])
        {
            nLevels = i;
            break;
        }
    }
    return nLevels;
}

namespace amrex {
    
    // helper to track names and centerings of the avaliable arrays
    class DescriptorMap : public amrex::StateMap
    {
    public:
        int Initialize(const DescriptorList &descriptors);
    };

    // --------------------------------------------------------------------------
    int DescriptorMap::Initialize(const DescriptorList &descriptors)
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

    // data adaptor's internal data
    struct AmrCatalystDataAdaptor::InternalsType
    {
        InternalsType() : amr(nullptr), amrMesh(nullptr) {}

        amrex::Amr *amr;
        vtkOverlappingAMR *amrMesh;
        amrex::DescriptorMap descriptorMap;
    };

    //-----------------------------------------------------------------------------
    AmrCatalystDataAdaptor::AmrCatalystDataAdaptor() :
            Internals(new AmrCatalystDataAdaptor::InternalsType())
    {
    }

//-----------------------------------------------------------------------------
    AmrCatalystDataAdaptor::~AmrCatalystDataAdaptor()
    {
        delete this->Internals;
    }

//-----------------------------------------------------------------------------
    int AmrCatalystDataAdaptor::CoProcess(Amr *aamr)
    {
        if (doCoProcess())
        {
            amrex::Print() << "Catalyst Begin CoProcess..." << std::endl;
            auto t0 = std::chrono::high_resolution_clock::now();

            amrex::Print() << "AMRCatalystDataAdaptor::CoProcess" << std::endl;

            this->Internals->amr = aamr;

            vtkNew<vtkCPDataDescription> dataDescription;
            dataDescription->AddInput("amr-input");
            dataDescription->SetTimeData(this->Internals->amr->cumTime(), this->Internals->amr->levelSteps(0));
            if (this->Processor->RequestDataDescription(dataDescription) != 0)
            {
                vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
                int numberOfProcesses = controller->GetNumberOfProcesses();
                int processId = controller->GetLocalProcessId();

                amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels = this->Internals->amr->getAmrLevels();
                this->Internals->descriptorMap.Initialize(levels[0]->get_desc_lst());

                this->BuildGrid(processId);
                this->AddGhostCellsArray(processId);
                this->AddArrays(processId, levels[0]->get_desc_lst());

                vtkCPInputDataDescription* inputDataDescription = dataDescription->GetInputDescriptionByName("amrex-input");
                inputDataDescription->SetGrid(this->Internals->amrMesh);

                // Set whole extent
                int wholeExtent[6] = { 0, numberOfProcesses, 0, 1, 0, 1 };
                inputDataDescription->SetWholeExtent(wholeExtent);

                // CoProcess
                this->Processor->CoProcess(dataDescription);

                // Release data
                this->Internals->amr = nullptr;
                this->Internals->amrMesh = nullptr;
                this->Internals->descriptorMap.Clear();
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            auto dt = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
            amrex::Print() << "Catalyst CoProcess complete (" << dt.count() << " sec)" << std::endl;
        }

        return 0;
    }

//-----------------------------------------------------------------------------
    int AmrCatalystDataAdaptor::BuildGrid(int rank)
    {
        amrex::Print() << "AmrCatalystDataAdaptor::BuildGrid" << std::endl;

        // get levels
        amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels = this->Internals->amr->getAmrLevels();

        unsigned int nLevels = numActiveLevels(levels);

        // initialize new vtk datasets
        this->Internals->amrMesh = vtkOverlappingAMR::New();

        // num levels and blocks per level
        std::vector<int> nBlocks(nLevels);
        for (unsigned int i = 0; i < nLevels; ++i)
            nBlocks[i] = levels[i]->boxArray().size();

        this->Internals->amrMesh->Initialize(nLevels, nBlocks.data());

        // origin
        const amrex::RealBox& pd = levels[0]->Geom().ProbDomain();
        double origin[3] = {AMREX_ARLIM(pd.lo())};

        this->Internals->amrMesh->SetOrigin(origin);

        long gid = 0;
        for (unsigned int i = 0; i < nLevels; ++i)
        {
            // domain decomp
            const amrex::DistributionMapping &dmap = levels[i]->DistributionMap();

            // ghost zones
            amrex::MultiFab &state = levels[i]->get_new_data(0);
            unsigned int ng = state.nGrow();

            // spacing
            const amrex::Geometry &geom = levels[i]->Geom();
            double spacing [3] = {AMREX_ARLIM(geom.CellSize())};
            this->Internals->amrMesh->SetSpacing(i, spacing);

            // refinement ratio
            this->Internals->amrMesh->SetRefinementRatio(i, levels[i]->fineRatio()[0]);

            // loop over boxes
            const amrex::BoxArray& ba = levels[i]->boxArray();
            unsigned int nBoxes = ba.size();

            for (unsigned int j = 0; j < nBoxes; ++j)
            {
                // cell centered box
                amrex::Box cbox = ba[j];

                // cell centered dimensions
                int cboxLo[3] = {AMREX_ARLIM(cbox.loVect())};
                int cboxHi[3] = {AMREX_ARLIM(cbox.hiVect())};

                // vtk's representation of box metadata
                vtkAMRBox block(cboxLo, cboxHi);
                this->Internals->amrMesh->SetAMRBox(i, j, block);
                this->Internals->amrMesh->SetAMRBlockSourceIndex(i, j, gid++);

                // skip building a vtk amrMesh for the non local boxes
                if (dmap[j] != rank)
                    continue;

                // add ghost zones
                for (int q = 0; q < AMREX_SPACEDIM; ++q)
                    cbox.grow(q, ng);

                // node centered box
                amrex::Box nbox = surroundingNodes(cbox);

                // node centered dimensions
                int nboxLo[3] = {AMREX_ARLIM(nbox.loVect())};
                int nboxHi[3] = {AMREX_ARLIM(nbox.hiVect())};

                // new vtk uniform amrMesh, node centered
                vtkUniformGrid *ug = vtkUniformGrid::New();
                ug->SetOrigin(origin);
                ug->SetSpacing(spacing);
                ug->SetExtent(nboxLo[0], nboxHi[0],
                              nboxLo[1], nboxHi[1],
                              nboxLo[2], nboxHi[2]);

                // pass the block into vtk
                this->Internals->amrMesh->SetDataSet(i, j, ug);
                ug->Delete();
            }
        }

        return 0;
    }

//-----------------------------------------------------------------------------
    int AmrCatalystDataAdaptor::AddGhostCellsArray(int rank)
    {
        amrex::Print() << "AmrCatalystDataAdaptor::AddGhostCellsArray" << std::endl;

        // loop over levels
        amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels = this->Internals->amr->getAmrLevels();

        unsigned int nLevels = numActiveLevels(levels);

        std::vector<std::vector<unsigned char*>> masks(nLevels);
        for (unsigned int i = 0; i < nLevels; ++i)
        {
            // allocate mask arrays
            const amrex::BoxArray &boxes = levels[i]->boxArray();
            const amrex::DistributionMapping &dmap = levels[i]->DistributionMap();
            const amrex::Box &pdom = levels[i]->Domain();

            amrex::MultiFab& state = levels[i]->get_new_data(0);
            unsigned int ng = state.nGrow();

            std::vector<unsigned char*> mask;
            VTKGhostUtils::AllocateBoxArray<unsigned char>(rank, pdom, boxes, dmap, ng, mask);

            // mask ghost cells
            VTKGhostUtils::MaskGhostCells<unsigned char>(rank, pdom, boxes, dmap, ng, mask);

            // store mask array
            masks[i] = mask;
        }

        // loop over coarse levels
        unsigned int nCoarseLevels = nLevels - 1;
        for (unsigned int i = 0; i < nCoarseLevels; ++i)
        {
            int ii = i + 1;

            // mask regions covered by refinement
            amrex::MultiFab& state = levels[i]->get_new_data(0);
            unsigned int ng = state.nGrow();

            const amrex::Box &pdom = levels[i]->Domain();
            const amrex::BoxArray &cBoxes = levels[i]->boxArray();
            const amrex::DistributionMapping &cMap = levels[i]->DistributionMap();
            const amrex::BoxArray &fBoxes = levels[ii]->boxArray();
            amrex::IntVect fRefRatio = levels[i]->fineRatio();

            VTKGhostUtils::MaskCoveredCells<unsigned char>(rank, pdom, cBoxes, cMap, fBoxes, fRefRatio, ng, masks[i]);
        }

        // loop over levels
        for (unsigned int i = 0; i < nLevels; ++i)
        {
            const amrex::DistributionMapping &dMap = levels[i]->DistributionMap();

            // mask arrays for this level
            std::vector<unsigned char*> &mask = masks[i];

            // loop over boxes
            const amrex::BoxArray& ba = levels[i]->boxArray();
            unsigned int nBoxes = ba.size();

            for (unsigned int j = 0; j < nBoxes; ++j)
            {
                // skip non-local blocks
                if (dMap[j] != rank)
                    continue;

                vtkUniformGrid *blockMesh = this->Internals->amrMesh->GetDataSet(i, j);

                if (!blockMesh)
                {
                    amrex::Print() << "Error @ AmrCatalystDataAdaptor::AddGhostCellsArray: "
                            <<" Empty block " << i << ", " << j
                            << std::endl;
                    return -1;
                }

                long nCells = blockMesh->GetNumberOfCells();

                // transfer mask array into vtk
                vtkUnsignedCharArray *ga = vtkUnsignedCharArray::New();
                ga->SetName("vtkGhostType");
                ga->SetArray(mask[j], nCells, 0);
                blockMesh->GetCellData()->AddArray(ga);
                ga->Delete();
            }
        }

        return 0;
    }

//-----------------------------------------------------------------------------
    int AmrCatalystDataAdaptor::AddArray(int rank, int association, const std::string &arrayName)
    {
        amrex::Print() << "AmrCatalystDataAdaptor::AddArray" << std::endl;

        if ((association != vtkDataObject::CELL) &&
            (association != vtkDataObject::POINT))
        {
            amrex::Print() << "Error @ AmrCatalystDataAdaptor::AddArray: "
                    << "Invalid association " << association
                    << std::endl;
            return -1;
        }

        amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels = this->Internals->amr->getAmrLevels();

        // find the indices of the multifab and component within for
        // the named array
        int fab = 0;
        int comp = 0;
        if (this->Internals->descriptorMap.GetIndex(arrayName, association, fab, comp))
        {
            amrex::Print() << "Error @ AmrCatalystDataAdaptor::AddArray: "
                    << "Failed to locate descriptor for "
                    << this->Internals->descriptorMap.GetAttributesName(association)
                    << " data array " << arrayName
                    << std::endl;
            return -1;
        }

        // loop over levels
        unsigned int nLevels = numActiveLevels(levels);
        for (unsigned int i = 0; i < nLevels; ++i)
        {
            // domain decomp
            const amrex::DistributionMapping &dmap = levels[i]->DistributionMap();

            // ghost zones
            amrex::MultiFab& state = levels[i]->get_new_data(fab);
            unsigned int ng = state.nGrow();

            if (!((association == vtkDataObject::CELL) && state.is_cell_centered()) &&
                !((association == vtkDataObject::POINT) && state.is_nodal()))
            {
                amrex::Print() << "Error @ AmrCatalystDataAdaptor::AddArray: "
                        << "association does not match MultiFAB centering"
                        << std::endl;
                return -1;
            }

            // check component id
            int nComp = state.nComp();
            if (comp >= nComp)
            {
                amrex::Print() << "Error @ AmrCatalystDataAdaptor::AddArray: "
                        << "Component " << comp << " out of bounds"
                        << std::endl;
                return -1;
            }

            // loop over boxes
            const amrex::BoxArray& ba = levels[i]->boxArray();
            unsigned int nBoxes = ba.size();

            for (unsigned int j = 0; j < nBoxes; ++j)
            {
                // cell centered box
                amrex::Box cbox = ba[j];

                // add ghost zones
                for (int q = 0; q < AMREX_SPACEDIM; ++q)
                    cbox.grow(q, ng);

                // cell centered dimensions
                int cboxLo[3] = {AMREX_ARLIM(cbox.loVect())};
                int cboxHi[3] = {AMREX_ARLIM(cbox.hiVect())};

                // skip building a vtk mesh for the non local boxes
                if (dmap[j] != rank)
                    continue;

                // node centered box
                amrex::Box nbox = surroundingNodes(cbox);

                // node centered dimensions
                int nboxLo[3] = {AMREX_ARLIM(nbox.loVect())};
                int nboxHi[3] = {AMREX_ARLIM(nbox.hiVect())};

                // get the block mesh
                vtkUniformGrid *ug = this->Internals->amrMesh->GetDataSet(i, j);

                // node centered size
                long nlen = 1;
                for (int p = 0; p < 3; ++p)
                    nlen *= nboxHi[p] - nboxLo[p] + 1;

                // cell centered size
                long clen = 1;
                for (int p = 0; p < 3; ++p)
                    clen *= cboxHi[p] - cboxLo[p] + 1;

                // pointer to the data
                amrex_real *pcd = state[j].dataPtr(comp);

                // allocate vtk array
                VTKGhostUtils::amrex_tt<amrex_real>::vtk_type *da =
                        VTKGhostUtils::amrex_tt<amrex_real>::vtk_type::New();

                // set component name
                da->SetName(arrayName.c_str());

                if (state[j].box().ixType() == amrex::IndexType::TheCellType())
                {
                    // zero copy cell centered
                    da->SetArray(pcd, clen, 1);
                    ug->GetCellData()->AddArray(da);
                }
                else if (state[j].box().ixType() == amrex::IndexType::TheNodeType())
                {
                    // zero copy point centered
                    da->SetArray(pcd, nlen, 1);
                    ug->GetPointData()->AddArray(da);
                }
                else
                {
                    amrex::Print() << "Warning @ AmrCatalystDataAdaptor::AddArray: "
                                   << "Face or edge centered component " << comp << " skipped"
                                   << std::endl;
                }

                da->Delete();

            }
        }

        return 0;
    }

//-----------------------------------------------------------------------------
    int AmrCatalystDataAdaptor::AddArrays(int rank, const DescriptorList &descriptors) {
        amrex::Print() << "AmrCatalystDataAdaptor::AddArrays" << std::endl;

        int ndesc = descriptors.size();
        for (int i = 0; i < ndesc; ++i) {

            const StateDescriptor &desc = descriptors[i];
            int ncomp = desc.nComp();
            IndexType itype = desc.getType();

            for (int j = 0; j < ncomp; ++j) {

                std::string arrayName = desc.name(j);

                if (itype.cellCentered()) {

                    this->AddArray(rank, vtkDataObject::CELL, arrayName);

                } else if (itype.nodeCentered()) {

                    this->AddArray(rank, vtkDataObject::POINT, arrayName);

                }
            }
        }

        return 0;
    }

}
