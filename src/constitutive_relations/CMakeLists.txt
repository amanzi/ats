# -*- mode: cmake -*-

# operators -- layer between discretization and PK
add_subdirectory(eos)
add_subdirectory(surface_subsurface_fluxes)
add_subdirectory(generic_evaluators)

#================================================================
# register evaluators/factories/pks

# eos/

register_evaluator_with_factory(
  HEADERFILE eos/eos_constant_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/eos_linear_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/eos_evaluator_tp_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/eos_evaluator_ctp_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/eos_ice_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )
register_evaluator_with_factory(
  HEADERFILE eos/eos_ideal_gas_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/eos_vapor_in_gas_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/eos_sw_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/eos_water_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/isobaric_eos_evaluator_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/molar_fraction_gas_evaluator_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/vapor_pressure_water_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/viscosity_constant_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/viscosity_evaluator_reg.hh
  LISTNAME   ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE eos/viscosity_water_reg.hh
  LISTNAME   ATS_RELATIONS_REG
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
  HEADERFILE generic_evaluators/MultiplicativeEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/AdditiveEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/SubgridDisaggregateEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/ColumnSumEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

register_evaluator_with_factory(
  HEADERFILE generic_evaluators/SubgridAggregateEvaluator_reg.hh
  LISTNAME ATS_RELATIONS_REG
  )

generate_evaluators_registration_header(
  HEADERFILE ats_relations_registration.hh
  LISTNAME   ATS_RELATIONS_REG
  INSTALL    True
  )
