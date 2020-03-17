#include "AMReX_AmrMeshCatalystDataAdaptor.H"

#include <chrono>
#include <timer/Timer.h>

#include <vtkObjectFactory.h>
#include <vtkOverlappingAMR.h>
#include <vtkAMRBox.h>
#include <vtkUniformGrid.h>
#include <vtkDataSetAttributes.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include <AMReX_ParmParse.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_Box.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_MFIter.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IndexType.H>
#include <AMReX_VTKGhostUtils.H>
#include <AMReX_Print.H>

namespace amrex {

    // helper to track names and centerings of the avaliable arrays
    class DescriptorMap : public amrex::StateMap
    {
    public:
        int Initialize(
                const std::vector<amrex::Vector<amrex::MultiFab> *> &states,
                const std::vector<std::vector<std::string>> &names);
    };

// --------------------------------------------------------------------------
    int DescriptorMap::Initialize(
            const std::vector<amrex::Vector<amrex::MultiFab> *> &states,
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

// data adaptor's internal data
    struct AmrMeshCatalystDataAdaptor::InternalsType
    {
        InternalsType() : Mesh(nullptr), States(), descriptorMap() {}

        amrex::AmrMesh *Mesh;
        vtkOverlappingAMR *amrMesh;
        std::vector<amrex::Vector<amrex::MultiFab> *> States;
        amrex::DescriptorMap descriptorMap;
    };

//-----------------------------------------------------------------------------
    AmrMeshCatalystDataAdaptor::AmrMeshCatalystDataAdaptor() :
            Internals(new AmrMeshCatalystDataAdaptor::InternalsType())
    {
    }

//-----------------------------------------------------------------------------
    AmrMeshCatalystDataAdaptor::~AmrMeshCatalystDataAdaptor()
    {
        delete this->Internals;
    }

    int AmrMeshCatalystDataAdaptor::CoProcess(unsigned int step, double time, amrex::AmrMesh *mesh,
                                              const std::vector<amrex::Vector<amrex::MultiFab>*> &states,
                                              const std::vector<std::vector<std::string>> &names)
    {
        int ret = 0;
        if (doCoProcess())
        {
            amrex::Print() << "Catalyst Begin CoProcess..." << std::endl;
            auto t0 = std::chrono::high_resolution_clock::now();

            timer::MarkEvent event("AmrMeshCatalystDataAdaptor::CoProcess");

            vtkNew<vtkCPDataDescription> dataDescription;
            dataDescription->AddInput("amrmesh-input");
            dataDescription->SetTimeData(time, step);
            if (this->Processor->RequestDataDescription(dataDescription) != 0)
            {
                vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
                int numberOfProcesses = controller->GetNumberOfProcesses();
                int processId = controller->GetLocalProcessId();

                this->Internals->Mesh = mesh;
                this->Internals->States = states;
                this->Internals->descriptorMap.Initialize(states, names);

                this->BuildGrid();
                this->AddGhostCellsArray();
                this->AddArrays(states, names);

                vtkCPInputDataDescription* inputDataDescription = dataDescription->GetInputDescriptionByName("amrex-input");
                inputDataDescription->SetGrid(this->Internals->amrMesh);

                // Set whole extent - ???
                int wholeExtent[6] = { 0, numberOfProcesses, 0, 1, 0, 1 };
                inputDataDescription->SetWholeExtent(wholeExtent);

                // CoProcess
                this->Processor->CoProcess(dataDescription);

                // Release data
                this->Internals->Mesh = nullptr;
                this->Internals->amrMesh = nullptr;
                this->Internals->States.clear();
                this->Internals->descriptorMap.Clear();
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            auto dt = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
            amrex::Print() << "Catalyst CoProcess complete (" << dt.count() << " sec)" << std::endl;
        }
        return ret;
    }

//-----------------------------------------------------------------------------
    int AmrMeshCatalystDataAdaptor::BuildGrid()
    {
        timer::MarkEvent("AmrMeshCatalystDataAdaptor::BuildGrid");

        unsigned int nLevels = this->Internals->Mesh->finestLevel() + 1;

        // initialize new vtk datasets
        this->Internals->amrMesh = vtkOverlappingAMR::New();

        // num levels and blocks per level
        std::vector<int> nBlocks(nLevels);
        for (unsigned int i = 0; i < nLevels; ++i)
            nBlocks[i] = this->Internals->Mesh->boxArray(i).size();

        this->Internals->amrMesh->Initialize(nLevels, nBlocks.data());

        // origin
        const amrex::RealBox& pd = this->Internals->Mesh->Geom(0).ProbDomain();
        double origin[3] = {AMREX_ARLIM(pd.lo())};

        this->Internals->amrMesh->SetOrigin(origin);

        long gid = 0;
        for (unsigned int i = 0; i < nLevels; ++i)
        {
            // domain decomp
            const amrex::DistributionMapping &dmap = this->Internals->Mesh->DistributionMap(i);

            // ghost zones
            amrex::MultiFab& state = this->Internals->States[0]->at(i);
            unsigned int ng = state.nGrow();

            // spacing
            const amrex::Geometry &geom = this->Internals->Mesh->Geom(i);
            double spacing [3] = {AMREX_ARLIM(geom.CellSize())};
            this->Internals->amrMesh->SetSpacing(i, spacing);

            // refinement ratio
            int cRefRatio = nLevels > 1 ? this->Internals->Mesh->refRatio(i)[0] : 1;
            this->Internals->amrMesh->SetRefinementRatio(i, cRefRatio);

            // loop over boxes
            const amrex::BoxArray& ba = this->Internals->Mesh->boxArray(i);
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
    int AmrMeshCatalystDataAdaptor::AddGhostCellsArray()
    {
        timer::MarkEvent("AmrMeshCatalystDataAdaptor::AddGhostCellsArray");

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // loop over levels
        unsigned int nLevels = this->Internals->Mesh->finestLevel() + 1;

        std::vector<std::vector<unsigned char*>> masks(nLevels);
        for (unsigned int i = 0; i < nLevels; ++i)
        {
            // allocate mask arrays
            const amrex::BoxArray &boxes = this->Internals->Mesh->boxArray(i);
            const amrex::DistributionMapping &dmap = this->Internals->Mesh->DistributionMap(i);
            const amrex::Box &pdom = this->Internals->Mesh->Geom(i).Domain();

            amrex::MultiFab& state = this->Internals->States[0]->at(i);
            unsigned int ng = state.nGrow();

            std::vector<unsigned char*> mask;
            VTKGhostUtils::AllocateBoxArray<unsigned char>(
                    pdom, boxes, dmap, ng, mask);

            // mask ghost cells
            VTKGhostUtils::MaskGhostCells<unsigned char>(
                    pdom, boxes, dmap, ng, mask);

            // store mask array
            masks[i] = mask;
        }

        // loop over coarse levels
        unsigned int nCoarseLevels = nLevels - 1;
        for (unsigned int i = 0; i < nCoarseLevels; ++i)
        {
            int ii = i + 1;

            // mask regions covered by refinement
            amrex::MultiFab& state = this->Internals->States[0]->at(i);
            unsigned int ng = state.nGrow();

            const amrex::Box &pdom = this->Internals->Mesh->Geom(i).Domain();
            const amrex::BoxArray &cBoxes = this->Internals->Mesh->boxArray(i);
            const amrex::DistributionMapping &cMap = this->Internals->Mesh->DistributionMap(i);
            const amrex::BoxArray &fBoxes = this->Internals->Mesh->boxArray(ii);
            amrex::IntVect cRefRatio = this->Internals->Mesh->refRatio(i);

            VTKGhostUtils::MaskCoveredCells<unsigned char>(
                    pdom, cBoxes, cMap, fBoxes, cRefRatio, ng, masks[i]);
        }

        // loop over levels
        for (unsigned int i = 0; i < nLevels; ++i)
        {
            const amrex::DistributionMapping &dmap = this->Internals->Mesh->DistributionMap(i);

            // mask arrays for this level
            std::vector<unsigned char*> &mask = masks[i];

            // loop over boxes
            const amrex::BoxArray& ba = this->Internals->Mesh->boxArray(i);
            unsigned int nBoxes = ba.size();

            for (unsigned int j = 0; j < nBoxes; ++j)
            {
                if (dmap[j] != rank)
                    continue;

                vtkUniformGrid *blockMesh = this->Internals->amrMesh->GetDataSet(i, j);

                if (!blockMesh)
                {
                    amrex::Print() << "Error @ AmrMeshCatalystDataAdaptor::AddGhostCellsArray: "
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
    int AmrMeshCatalystDataAdaptor::AddArray(int association, const std::string &arrayName)
    {
        timer::MarkEvent("AmrMeshCatalystDataAdaptor::AddArray");

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if ((association != vtkDataObject::CELL) &&
            (association != vtkDataObject::POINT))
        {
            amrex::Print() << "Error @ AmrMeshCatalystDataAdaptor::AddArray: "
                           << "Invalid association " << association
                           << std::endl;
            return -1;
        }

        // find the indices of the multifab and component within for
        // the named array
        int fab = 0;
        int comp = 0;
        if (this->Internals->descriptorMap.GetIndex(arrayName, association, fab, comp))
        {
            amrex::Print() << "Error @ AmrMeshCatalystDataAdaptor::AddArray: "
                           << "Failed to locate descriptor for "
                           << amrex::StateMap::GetAttributesName(association)
                           << " data array " << arrayName
                           << std::endl;
            return -1;
        }

        // loop over levels
        unsigned int nLevels = this->Internals->Mesh->finestLevel() + 1;

        for (unsigned int i = 0; i < nLevels; ++i)
        {
            // domain decomp
            const amrex::DistributionMapping &dmap = this->Internals->Mesh->DistributionMap(i);

            // ghost zones
            amrex::MultiFab& state = this->Internals->States[fab]->at(i);
            unsigned int ng = state.nGrow();

            // check centering
            if (!((association == vtkDataObject::CELL) && state.is_cell_centered()) &&
                !((association == vtkDataObject::POINT) && state.is_nodal()))
            {
                amrex::Print() << "Error @ AmrMeshCatalystDataAdaptor::AddArray: "
                               << "association does not match MultiFAB centering"
                               << std::endl;
                return -1;
            }

            // check component id
            int nComp = state.nComp();
            if (comp >= nComp)
            {
                amrex::Print() << "Error @ AmrMeshCatalystDataAdaptor::AddArray: "
                               << "Component " << comp << " out of bounds"
                               << std::endl;
                return -1;
            }

            // loop over boxes
            const amrex::BoxArray& ba = this->Internals->Mesh->boxArray(i);
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
                    amrex::Print() << "Warning @ AmrMeshCatalystDataAdaptor::AddArray: "
                                   << "Face or edge centered component " << comp << " skipped"
                                   << std::endl;
                }

                da->Delete();
            }
        }

        return 0;
    }

    //-----------------------------------------------------------------------------
    int AmrMeshCatalystDataAdaptor::AddArrays(const std::vector<amrex::Vector<amrex::MultiFab> *> &states,
                                              const std::vector<std::vector<std::string>> &names)
    {
        timer::MarkEvent("AmrMeshCatalystDataAdaptor::AddArrays");

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
                    this->AddArray(vtkDataObject::CELL, arrayName);
                }
                else if (state.is_nodal())
                {
                    this->AddArray(vtkDataObject::POINT, arrayName);
                }
            }
        }

        return 0;
    }
}
