# -*- mode: cmake -*-

#
#  Upwinding Evaluators
#
#  Generic (non-flow-specific) evaluators that upwind a cell-based field onto
#  faces for use as a div-grad operator's nonlinear coefficient.
#
set(ats_upwinding_evals_src_files
    evaluator_upwinding.cc
    evaluator_upwind_total_flux.cc
    evaluator_upwind_arithmetic_mean.cc
    evaluator_upwind_potential_difference.cc
   )

file(GLOB ats_upwinding_evals_inc_files "*.hh")

set(ats_upwinding_evals_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  state
  )

add_amanzi_library(ats_upwinding_evals
                   SOURCE ${ats_upwinding_evals_src_files}
                   HEADERS ${ats_upwinding_evals_inc_files}
                   LINK_LIBS ${ats_upwinding_evals_link_libs})

register_evaluator_with_factory(
  HEADERFILE evaluator_upwind_total_flux_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE evaluator_upwind_arithmetic_mean_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE evaluator_upwind_potential_difference_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )
