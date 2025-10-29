# -*- mode: cmake -*-

# operators -- layer between discretization and PK
add_subdirectory(eos)
add_subdirectory(surface_subsurface_fluxes)
add_subdirectory(generic_evaluators)
add_subdirectory(column_integrators)

#================================================================
# register evaluators/factories/pks

# eos/

register_evaluator_with_factory(
  HEADERFILE eos/eos_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )


register_evaluator_with_factory(
  HEADERFILE column_integrators/column_integrators_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

# constitutive_relations/surface_subsurface_fluxes/

register_evaluator_with_factory(
  HEADERFILE surface_subsurface_fluxes/volumetric_darcy_flux_evaluator_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE surface_subsurface_fluxes/surface_top_cells_evaluator_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE surface_subsurface_fluxes/top_cells_surface_evaluator_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE surface_subsurface_fluxes/overland_source_from_subsurface_flux_evaluator_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE surface_subsurface_fluxes/flux_divergence_evaluator_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )


register_evaluator_with_factory(
  HEADERFILE generic_evaluators/MultiplicativeEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/AdditiveEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/ReciprocalEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/SubgridDisaggregateEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/SubgridAggregateEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/TimeMaxEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/InitialTimeEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/ExtractionEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

generate_evaluators_registration_header(
  HEADERFILE ats_relations_registration.hh
  LISTNAME   ATS_RELATIONS_REG
  INSTALL    True
  )
