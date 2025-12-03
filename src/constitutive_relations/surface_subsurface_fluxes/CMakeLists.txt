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


