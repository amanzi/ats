# -*- mode: cmake -*-

#
#  Generic Evaluators 
#
set(ats_generic_evals_src_files
    MultiplicativeEvaluator.cc
    AdditiveEvaluator.cc
    SubgridDisaggregateEvaluator.cc
    SubgridAggregateEvaluator.cc
    ColumnSumEvaluator.cc	
   )

file(GLOB ats_generic_evals_inc_files "*.hh")

set(ats_generic_evals_link_libs
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

add_amanzi_library(ats_generic_evals
                   SOURCE ${ats_generic_evals_src_files}
                   HEADERS ${ats_generic_evals_inc_files}
		   LINK_LIBS ${ats_generic_evals_link_libs})


