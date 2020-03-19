find_package(ParaView 5.8 REQUIRED)

target_link_libraries(amrex
   PUBLIC
   ParaView::Catalyst
   ParaView::PythonCatalyst
   VTK::CommonCore
   VTK::CommonDataModel
   VTK::ParallelCore)

target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
target_compile_definitions(amrex PUBLIC $<BUILD_INTERFACE:BL_USE_CATALYST_INSITU>)

target_sources(amrex
   PRIVATE
   AMReX_AmrCatalystDataAdaptor.H
   AMReX_AmrMeshCatalystDataAdaptor.H
   AMReX_CatalystDataAdaptor.H
   AMReX_StateMap.H
   AMReX_VTKGhostUtils.H
   AMReX_AmrCatalystDataAdaptor.cpp
   AMReX_AmrMeshCatalystDataAdaptor.cpp
   AMReX_CatalystDataAdaptor.cpp
   AMReX_StateMap.cpp
   )
