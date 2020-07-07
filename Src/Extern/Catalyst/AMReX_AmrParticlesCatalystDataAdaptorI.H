
#include <chrono>
#include <map>
#include <utility>

#include <vtkAMRBox.h>                 // VTK::CommonDataModel
#include <vtkCellData.h>               // VTK::CommonDataModel
#include <vtkCellArray.h>              // VTK::CommonDataModel
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
#include <vtkPolyData.h>               // VTK::CommonDataModel
#include <vtkCPProcessor.h>            // ParaView::Catalyst
#include <vtkUniformGrid.h>            // VTK::CommonDataModel
#include <vtkUnsignedCharArray.h>      // VTK::CommonCore

#include <AMReX_Amr.H>
#include <AMReX_AmrDescriptorMap.H>
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
#include <AMReX_Print.H>

namespace amrex {
    //-----------------------------------------------------------------------------
     template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
     int AmrParticlesCatalystDataAdaptor::CoProcess(Amr *aamr,
             ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer,
             std::vector<std::string> real_comp_names,
             std::vector<std::string> int_comp_names)
    {
        int ret = 0;
        if (doCoProcess())
        {
            amrex::Print() << "Catalyst Begin CoProcess..." << std::endl;
            auto t0 = std::chrono::high_resolution_clock::now();

            this->amr = aamr;

            vtkNew<vtkCPDataDescription> dataDescription;
            dataDescription->AddInput("inputmesh");
            dataDescription->AddInput("inputparticles");
            dataDescription->SetTimeData(this->amr->cumTime(), this->amr->levelSteps(0));
            if (this->Processor->RequestDataDescription(dataDescription) != 0)
            {
                vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
                int numberOfProcesses = controller->GetNumberOfProcesses();
                int processId = controller->GetLocalProcessId();

                amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels = this->amr->getAmrLevels();

                this->descriptorMap = new AmrDescriptorMap;
                this->descriptorMap->Initialize(levels[0]->get_desc_lst());

                this->amrMesh = vtkOverlappingAMR::New();
                this->BuildGrid(processId);
                this->AddGhostCellsArray(processId);
                this->AddArrays(processId, levels[0]->get_desc_lst());

                vtkCPInputDataDescription* inputMeshDataDescription = dataDescription->GetInputDescriptionByName("inputmesh");
                inputMeshDataDescription->SetGrid(this->amrMesh);

                int wholeExtent[6] = { 0, numberOfProcesses, 0, 1, 0, 1 };
                inputMeshDataDescription->SetWholeExtent(wholeExtent);

                this->particleData = vtkPolyData::New();
                this->numberOfParticles = this->CountParticles(processId, particleContainer);
                this->BuildParticles(processId, particleContainer);
                this->AddParticlesArrays(processId, particleContainer, real_comp_names, int_comp_names);

                vtkCPInputDataDescription* inputParticlesDataDescription = dataDescription->GetInputDescriptionByName("inputparticles");
                inputParticlesDataDescription->SetGrid(this->particleData);

                //int particlesWholeExtent[6] = { 0, numberOfProcesses, 0, 1, 0, 1 };
                //inputParticlesDataDescription->SetWholeExtent(particlesWholeExtent);

                // CoProcess
                this->Processor->CoProcess(dataDescription);

                // Release data
                this->amr = nullptr;
                this->amrMesh->Delete();
                this->amrMesh = nullptr;
                this->descriptorMap->Clear();
                this->descriptorMap = nullptr;
                this->particleData->Delete();
                this->particleData = nullptr;
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            auto dt = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
            amrex::Print() << "Catalyst CoProcess complete (" << dt.count() << " sec)" << std::endl;
        }

        return ret;
    }

    //-----------------------------------------------------------------------------
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    int AmrParticlesCatalystDataAdaptor::CountParticles(
            int rank,
            ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer)
    {
        const auto& particles = particleContainer.GetParticles();
        int count = 0;
        for (int lev = 0; lev < particles.size();  lev++)
        {
            auto& pmap = particles[lev];
            for (const auto& kv : pmap)
            {
                const auto& aos = kv.second.GetArrayOfStructs();
                for (auto it = aos.cbegin(); it != aos.cend(); ++it)
                {
                    if (it->m_idata.id >= 0)
                    {
                        count++;
                    }
                }
            }
        }
        return count;
    }

    //-----------------------------------------------------------------------------
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    int AmrParticlesCatalystDataAdaptor::BuildParticles(
            int rank,
            ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer)
    {
        const auto& particles = particleContainer.GetParticles();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(this->numberOfParticles);
        int count = 0;
        for (int lev = 0; lev < particles.size();  lev++)
        {
            auto& pmap = particles[lev];
            for (const auto& kv : pmap)
            {
                const auto& aos = kv.second.GetArrayOfStructs();
                for (auto it = aos.cbegin(); it != aos.cend(); ++it)
                {
                    if (it->m_idata.id > 0)
                    {
                        if (AMREX_SPACEDIM == 2)
                            points->SetPoint(count, it->m_rdata.pos[0], it->m_rdata.pos[1], 0.0);
                        if (AMREX_SPACEDIM == 3)
                            points->SetPoint(count, it->m_rdata.pos[0], it->m_rdata.pos[1], it->m_rdata.pos[2]);
                        count++;
                    }
                }
            }
        }
        this->particleData->SetPoints(points);
        vtkSmartPointer<vtkCellArray> vertex = vtkSmartPointer<vtkCellArray>::New();
        vertex->Allocate(this->numberOfParticles * 2);
        for (int i = 0; i < this->numberOfParticles; i++)
        {
            vertex->InsertNextCell(1);
            vertex->InsertCellPoint(i);
        }
        this->particleData->SetVerts(vertex);

        return 0;
    }

    //-----------------------------------------------------------------------------
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    int AmrParticlesCatalystDataAdaptor::AddParticlesIDArray(
            int rank,
            ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer)
    {
        const auto& particles = particleContainer.GetParticles();
        vtkSmartPointer<vtkIntArray> intArray = vtkSmartPointer<vtkIntArray>::New();
        intArray->SetName("id");
        intArray->SetNumberOfComponents(1);
        intArray->SetNumberOfValues(this->numberOfParticles);
        int count = 0;
        for (int level = 0; level<particles.size();level++)
        {
            auto &pmap = particles[level];
            for (const auto &kv : pmap)
            {
                const auto &aos = kv.second.GetArrayOfStructs();
                for (auto it = aos.cbegin();it != aos.cend();++it)
                {
                    if (it->m_idata.id > 0)
                    {
                        intArray->SetValue(count, it->m_idata.id);
                        count++;
                    }
                }
            }
        }
        this->particleData->GetPointData()->AddArray(intArray);

        return 0;
    }

    //-----------------------------------------------------------------------------
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    int AmrParticlesCatalystDataAdaptor::AddParticlesCPUArray(
            int rank,
            const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer)
    {
        const auto& particles = particleContainer.GetParticles();
        vtkSmartPointer<vtkIntArray> intArray = vtkSmartPointer<vtkIntArray>::New();
        intArray->SetName("cpu");
        intArray->SetNumberOfComponents(1);
        intArray->SetNumberOfValues(this->numberOfParticles);
        int count = 0;
        for (int level = 0; level<particles.size();level++)
        {
            auto &pmap = particles[level];
            for (const auto &kv : pmap)
            {
                const auto &aos = kv.second.GetArrayOfStructs();
                for (auto it = aos.cbegin();it != aos.cend();++it)
                {
                    if (it->m_idata.id > 0)
                    {
                        intArray->SetValue(count, it->m_idata.cpu);
                        count++;
                    }
                }
            }
        }
        this->particleData->GetPointData()->AddArray(intArray);

        return 0;
    }

    //-----------------------------------------------------------------------------
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    int AmrParticlesCatalystDataAdaptor::AddParticlesAOSRealArray(
            int rank,
            int index,
            const std::string &arrayName,
            const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer)
    {
        const auto& particles = particleContainer.GetParticles();
        vtkSmartPointer<vtkDoubleArray> doubleArray = vtkSmartPointer<vtkDoubleArray>::New();
        doubleArray->SetName(arrayName.c_str());
        doubleArray->SetNumberOfComponents(1);
        doubleArray->SetNumberOfValues(this->numberOfParticles);
        int count = 0;
        for (int level = 0; level<particles.size();level++)
        {
            auto &pmap = particles[level];
            for (const auto &kv : pmap)
            {
                const auto &aos = kv.second.GetArrayOfStructs();
                for (auto it = aos.cbegin();it != aos.cend();++it)
                {
                    if (it->m_idata.id > 0)
                    {
                        doubleArray->SetValue(count, it->m_rdata.arr[AMREX_SPACEDIM + index]);
                        count++;
                    }
                }
            }
        }
        this->particleData->GetPointData()->AddArray(doubleArray);

        return 0;
    }

    //-----------------------------------------------------------------------------
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    int AmrParticlesCatalystDataAdaptor::AddParticlesAOSIntArray(
            int rank,
            int index,
            const std::string &arrayName,
            const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer)
    {
        const auto& particles = particleContainer.GetParticles();
        vtkSmartPointer<vtkIntArray> intArray = vtkSmartPointer<vtkIntArray>::New();
        intArray->SetName(arrayName.c_str());
        intArray->SetNumberOfComponents(1);
        intArray->SetNumberOfValues(this->numberOfParticles);
        int count = 0;
        for (int level = 0; level<particles.size();level++)
        {
            auto &pmap = particles[level];
            for (const auto &kv : pmap)
            {
                const auto &aos = kv.second.GetArrayOfStructs();
                for (auto it = aos.cbegin();it != aos.cend();++it)
                {
                    if (it->m_idata.id > 0)
                    {
                        intArray->SetValue(count, it->m_idata.arr[AMREX_SPACEDIM + index]);
                        count++;
                    }
                }
            }
        }
        this->particleData->GetPointData()->AddArray(intArray);

        return 0;
    }

    //-----------------------------------------------------------------------------
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    int AmrParticlesCatalystDataAdaptor::AddParticlesSOARealArray(
            int rank,
            int index,
            const std::string &arrayName,
            const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer)
    {
        const auto& particles = particleContainer.GetParticles();
        vtkSmartPointer<vtkDoubleArray> doubleArray = vtkSmartPointer<vtkDoubleArray>::New();
        doubleArray->SetName(arrayName.c_str());
        doubleArray->SetNumberOfComponents(1);
        doubleArray->SetNumberOfValues(this->numberOfParticles);
        int count = 0;
        for (int level = 0; level<particles.size();level++)
        {
            auto &pmap = particles[level];
            for (const auto &kv : pmap)
            {
                const auto &aos = kv.second.GetArrayOfStructs();
                const auto &soa = kv.second.GetStructOfArrays();
                int element = 0;
                for (auto it = aos.cbegin();it != aos.cend();++it)
                {
                    if (it->m_idata.id > 0)
                    {
                        doubleArray->SetValue(count, soa.GetRealData(index)[element]);
                        element++;
                        count++;
                    }
                }
            }
        }
        this->particleData->GetPointData()->AddArray(doubleArray);

        return 0;
    }

    //-----------------------------------------------------------------------------
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    int AmrParticlesCatalystDataAdaptor::AddParticlesSOAIntArray(
            int rank,
            int index,
            const std::string &arrayName,
            const ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer)
    {
        const auto& particles = particleContainer.GetParticles();
        vtkSmartPointer<vtkIntArray> intArray = vtkSmartPointer<vtkIntArray>::New();
        intArray->SetName(arrayName.c_str());
        intArray->SetNumberOfComponents(1);
        intArray->SetNumberOfValues(this->numberOfParticles);
        for (int level = 0; level<particles.size();level++)
        {
            auto &pmap = particles[level];
            for (const auto &kv : pmap)
            {
                const auto &aos = kv.second.GetArrayOfStructs();
                const auto &soa = kv.second.GetStructOfArrays();
                int count = 0;
                for (auto it = aos.cbegin();it != aos.cend();++it)
                {
                    if (it->m_idata.id > 0)
                    {
                        intArray->SetValue(soa.GetIntData(index)[count]);
                        count++;
                    }
                }
            }
        }
        this->particleData->GetPointData()->AddArray(intArray);

        return 0;
    }

    //-----------------------------------------------------------------------------
    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    int AmrParticlesCatalystDataAdaptor::AddParticlesArrays(
            int rank,
            ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> &particleContainer,
            std::vector<std::string> real_comp_names,
            std::vector<std::string> int_comp_names)
    {
        AddParticlesIDArray(rank, particleContainer);
        AddParticlesCPUArray(rank, particleContainer);
        for (int i = 0; i < NStructReal; ++i )
        {
            if (real_comp_names.size() == 0)
            {
                std::stringstream ss;
                ss << "real_comp" << i;
                this->AddParticlesAOSRealArray(rank, i, ss.str(), particleContainer);
            }
            else
            {
                this->AddParticlesAOSRealArray(rank, i, real_comp_names[i], particleContainer);
            }
        }
        for (int i = NStructReal; i < (NStructReal + NArrayReal); ++i )
        {
            if (real_comp_names.size() == 0)
            {
                std::stringstream ss;
                ss << "real_comp" << i;
                this->AddParticlesSOARealArray(rank, i - NStructReal, ss.str(), particleContainer);
            }
            else
            {
                this->AddParticlesSOARealArray(rank, i - NStructReal, real_comp_names[i], particleContainer);
            }
        }
        for (int i = 0; i < NStructInt; ++i )
        {
            if (int_comp_names.size() == 0)
            {
                std::stringstream ss;
                ss << "int_comp" << i;
                this->AddParticlesAOSIntArray(rank, i, ss.str(), particleContainer);
            }
            else
            {
                this->AddParticlesAOSIntArray(rank, i, int_comp_names[i], particleContainer);
            }
        }
        for (int i = NStructInt; i < (NStructInt + NArrayInt); ++i )
        {
            if (int_comp_names.size() == 0)
            {
                std::stringstream ss;
                ss << "int_comp" << i;
                this->AddParticlesSOARealArray(rank, i - NStructInt, ss.str(), particleContainer);
            }
            else
            {
                this->AddParticlesSOARealArray(rank, i - NStructInt, int_comp_names[i], particleContainer);
            }
        }

        return 0;
    }

}
