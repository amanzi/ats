# -*- mode: cmake -*-
# ATS Surface balance PKs describe Evaporation, energy fluxes from
#  long/showtwave radiation, precip, etc etc etc
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/constitutive_relations/surface_subsurface_fluxes)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/constitutive_relations/land_cover)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/constitutive_relations/litter)
include_directories(${CLM_INCLUDE_DIRS})
include_directories(${ATS_SOURCE_DIR}/src/pks/flow/constitutive_relations/elevation) # subgrid model

set(ats_surface_balance_src_files
  constitutive_relations/land_cover/seb_physics_defs.cc
  constitutive_relations/land_cover/seb_physics_funcs.cc
  constitutive_relations/land_cover/longwave_evaluator.cc
  constitutive_relations/land_cover/incident_shortwave_radiation_model.cc
  constitutive_relations/land_cover/incident_shortwave_radiation_evaluator.cc
  constitutive_relations/land_cover/evaporation_downregulation_model.cc
  constitutive_relations/land_cover/evaporation_downregulation_evaluator.cc
  constitutive_relations/land_cover/LandCover.cc
  constitutive_relations/land_cover/drainage_evaluator.cc
  constitutive_relations/land_cover/interception_fraction_evaluator.cc
  constitutive_relations/land_cover/interception_fraction_model.cc
  constitutive_relations/land_cover/pet_priestley_taylor_evaluator.cc
  constitutive_relations/land_cover/plant_wilting_factor_evaluator.cc
  constitutive_relations/land_cover/plant_wilting_factor_model.cc  
  constitutive_relations/land_cover/rooting_depth_fraction_evaluator.cc  
  constitutive_relations/land_cover/rooting_depth_fraction_model.cc
  constitutive_relations/land_cover/snow_meltrate_evaluator.cc
  constitutive_relations/land_cover/transpiration_distribution_evaluator.cc  
  constitutive_relations/land_cover/albedo_twocomponent_evaluator.cc
  constitutive_relations/land_cover/albedo_threecomponent_evaluator.cc
  constitutive_relations/land_cover/area_fractions_twocomponent_evaluator.cc
  constitutive_relations/land_cover/area_fractions_threecomponent_evaluator.cc
  constitutive_relations/land_cover/area_fractions_threecomponent_microtopography_evaluator.cc
  constitutive_relations/land_cover/radiation_balance_evaluator.cc
  constitutive_relations/land_cover/seb_twocomponent_evaluator.cc
  constitutive_relations/land_cover/seb_threecomponent_evaluator.cc
  constitutive_relations/litter/interception_evaluator.cc
  constitutive_relations/litter/latent_heat_evaluator.cc
  constitutive_relations/litter/latent_heat_model.cc       
  constitutive_relations/litter/evaporative_flux_relaxation_evaluator.cc       
  constitutive_relations/litter/evaporative_flux_relaxation_model.cc
  constitutive_relations/litter/micropore_macropore_flux_evaluator.cc       
  constitutive_relations/litter/micropore_macropore_flux_model.cc       
  constitutive_relations/litter/macropore_surface_flux_evaluator.cc       
  constitutive_relations/litter/macropore_surface_flux_model.cc       
  surface_balance_base.cc
  surface_balance_implicit_subgrid.cc
)

if (ENABLE_CLM)
  list(APPEND ats_surface_balance_src_files
    CLM/ats_clm.F90
    CLM/ats_clm_interface.cc
    CLM/surface_balance_CLM.cc
    )
endif()

set(ats_surface_balance_inc_files
  constitutive_relations/land_cover/seb_physics_defs.hh
  constitutive_relations/land_cover/seb_physics_funcs.hh
  constitutive_relations/land_cover/longwave_evaluator.hh
  constitutive_relations/land_cover/incident_shortwave_radiation_model.hh
  constitutive_relations/land_cover/incident_shortwave_radiation_evaluator.hh
  constitutive_relations/land_cover/evaporation_downregulation_model.hh
  constitutive_relations/land_cover/evaporation_downregulation_evaluator.hh
  constitutive_relations/land_cover/LandCover.hh
  constitutive_relations/land_cover/drainage_evaluator.hh
  constitutive_relations/land_cover/interception_fraction_evaluator.hh
  constitutive_relations/land_cover/interception_fraction_model.hh
  constitutive_relations/land_cover/pet_priestley_taylor_evaluator.hh
  constitutive_relations/land_cover/plant_wilting_factor_evaluator.hh
  constitutive_relations/land_cover/plant_wilting_factor_model.hh  
  constitutive_relations/land_cover/rooting_depth_fraction_evaluator.hh  
  constitutive_relations/land_cover/rooting_depth_fraction_model.hh
  constitutive_relations/land_cover/snow_meltrate_evaluator.hh
  constitutive_relations/land_cover/transpiration_distribution_evaluator.hh  
  constitutive_relations/land_cover/albedo_twocomponent_evaluator.hh
  constitutive_relations/land_cover/albedo_threecomponent_evaluator.hh
  constitutive_relations/land_cover/area_fractions_twocomponent_evaluator.hh
  constitutive_relations/land_cover/area_fractions_threecomponent_evaluator.hh
  constitutive_relations/land_cover/area_fractions_threecomponent_microtopography_evaluator.hh
  constitutive_relations/land_cover/radiation_balance_evaluator.hh
  constitutive_relations/land_cover/seb_twocomponent_evaluator.hh
  constitutive_relations/land_cover/seb_threecomponent_evaluator.hh
  constitutive_relations/litter/interception_evaluator.hh
  constitutive_relations/litter/latent_heat_evaluator.hh
  constitutive_relations/litter/latent_heat_model.hh       
  constitutive_relations/litter/evaporative_flux_relaxation_evaluator.hh       
  constitutive_relations/litter/evaporative_flux_relaxation_model.hh
  constitutive_relations/litter/micropore_macropore_flux_evaluator.hh       
  constitutive_relations/litter/micropore_macropore_flux_model.hh       
  constitutive_relations/litter/macropore_surface_flux_evaluator.hh       
  constitutive_relations/litter/macropore_surface_flux_model.hh       
  surface_balance_base.hh
  surface_balance_implicit_subgrid.hh
)  

if (ENABLE_CLM)
  list(APPEND ats_surface_balance_inc_files
    CLM/ats_clm.F90
    CLM/surface_balance_CLM.hh
    CLM/ats_clm_interface.hh
    CLM/ats_clm_interface_private.hh
    )
endif()

set(ats_surface_balance_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  ${CLM_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  operators
  solvers
  time_integration
  state
  pks
  ats_operators
  ats_pks
  )

add_amanzi_library(ats_surface_balance
                   SOURCE ${ats_surface_balance_src_files}
                   HEADERS ${ats_surface_balance_inc_files}
		   LINK_LIBS ${ats_surface_balance_link_libs})


#================================================
# register evaluators/factories/pks

# land_cover evalutors
register_evaluator_with_factory(
  HEADERFILE constitutive_relations/land_cover/land_cover_evaluator_reg.hh
  LISTNAME ATS_SURFACE_BALANCE_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/interception_evaluator_reg.hh
  LISTNAME ATS_SURFACE_BALANCE_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/latent_heat_evaluator_reg.hh
  LISTNAME ATS_SURFACE_BALANCE_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/micropore_macropore_flux_evaluator_reg.hh
  LISTNAME ATS_SURFACE_BALANCE_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/macropore_surface_flux_evaluator_reg.hh
  LISTNAME ATS_SURFACE_BALANCE_REG
)


register_evaluator_with_factory(
  HEADERFILE constitutive_relations/litter/evaporative_flux_relaxation_evaluator_reg.hh
  LISTNAME ATS_SURFACE_BALANCE_REG
)

if (ENABLE_CLM)
  register_evaluator_with_factory(
    HEADERFILE CLM/surface_balance_CLM_reg.hh
    LISTNAME ATS_SURFACE_BALANCE_REG
    )
endif()

register_evaluator_with_factory(
  HEADERFILE surface_balance_base_reg.hh
  LISTNAME ATS_SURFACE_BALANCE_REG
)
register_evaluator_with_factory(
  HEADERFILE surface_balance_implicit_subgrid_reg.hh
  LISTNAME ATS_SURFACE_BALANCE_REG
)
generate_evaluators_registration_header(
  HEADERFILE ats_surface_balance_registration.hh
  LISTNAME   ATS_SURFACE_BALANCE_REG
  INSTALL    True
  )
