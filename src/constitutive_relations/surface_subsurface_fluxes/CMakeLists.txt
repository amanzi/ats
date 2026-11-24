# -*- mode: cmake -*-

#
#  ATS
#


set(ats_surf_subsurf_src_files
  overland_source_from_subsurface_flux_evaluator.cc
  surface_top_cells_evaluator.cc
  top_cells_surface_evaluator.cc
  volumetric_darcy_flux_evaluator.cc
  flux_divergence_evaluator.cc
 )

file(GLOB ats_surf_subsurf_inc_files "*.hh")


set(ats_surf_subsurf_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  solvers
  state
  )


add_amanzi_library(ats_surf_subsurf
                   SOURCE ${ats_surf_subsurf_src_files}
                   HEADERS ${ats_surf_subsurf_inc_files}
		   LINK_LIBS ${ats_surf_subsurf_link_libs})


if (BUILD_TESTS)
  include_directories(${UnitTest_INCLUDE_DIRS})
  include_directories(${MESH_FACTORY_SOURCE_DIR})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  # FluxDivergenceEvaluator must not depend on the mesh partition: the np2 run
  # compares a distributed mesh against a serial copy cell by cell.
  add_amanzi_test(surf_subsurf_flux_divergence surf_subsurf_flux_divergence
    KIND int
    SOURCE test/Main.cc test/surf_subsurf_flux_divergence.cc
    LINK_LIBS ats_surf_subsurf mesh_factory ${ats_surf_subsurf_link_libs} ${UnitTest_LIBRARIES})
  add_amanzi_test(surf_subsurf_flux_divergence_np2 surf_subsurf_flux_divergence NPROCS 2 KIND int)
endif()


