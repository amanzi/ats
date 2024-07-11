# -*- mode: cmake -*-

#
#  Generic Evaluators 
#
set(ats_column_integrator_src_files
  ColumnSumEvaluator.cc
  activelayer_average_temp_evaluator.cc  
  water_table_depth_evaluator.cc
  perched_water_table_depth_evaluator.cc
  thaw_depth_evaluator.cc
  )


set(ats_column_integrator_inc_files
  EvaluatorColumnIntegrator.hh
  ColumnSumEvaluator.hh
  activelayer_average_temp_evaluator.hh
  water_table_depth_evaluator.hh
  perched_water_table_depth_evaluator.hh
  thaw_depth_evaluator.hh
  )

set(ats_column_integrator_link_libs
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

add_amanzi_library(ats_column_integrator
                   SOURCE ${ats_column_integrator_src_files}
                   HEADERS ${ats_column_integrator_inc_files}
		   LINK_LIBS ${ats_column_integrator_link_libs})


