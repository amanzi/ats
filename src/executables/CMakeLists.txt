# -*- mode: cmake -*-


# Need this define. Errors from MSTK include files 
# about MPI_COMM_WORLD. --lpritch
add_definitions("-DMSTK_HAVE_MPI")

if (WITH_MSTK_2_20rc1_OR_NEWER)
  add_definitions("-DMSTK_2_20rc1_OR_NEWER")
endif ()

if (WITH_MSTK_2_21rc1_OR_NEWER)
  add_definitions("-DMSTK_2_21rc1_OR_NEWER")
endif ()


#
#  ATS
#    Executable
#
include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/mpc)
include_directories(${ATS_SOURCE_DIR}/src/pks/coupled_transport)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad)
include_directories(${ATS_SOURCE_DIR}/src/operators/divgrad/upwind_scheme)
include_directories(${ATS_SOURCE_DIR}/src/operators/advection)
include_directories(${ATS_SOURCE_DIR}/src/operators/deformation)
include_directories(${ATS_SOURCE_DIR}/src/pks/energy/base)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport)
include_directories(${ATS_SOURCE_DIR}/src/pks/transport/transport_amanzi)
include_directories(${ATS_SOURCE_DIR}/src/pks/chemistry)
include_directories(${ATS_SOURCE_DIR}/src/pks/chemistry/chemistry_amanzi)

include_directories(${ATS_SOURCE_DIR}/src/coordinator)

#------------------------------------------------------------------------------#
# ATS F90 includes.
#------------------------------------------------------------------------------#

include_directories(${ATS_SOURCE_DIR}/src/xlib/c)
include_directories(${ATS_BINARY_DIR}/src/xlib/fc)
include_directories(${ATS_BINARY_DIR}/src/xlib)


include_evaluators_directories(LISTNAME ATS_POROSITY_EVAL_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_EOS_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_SURFACE_SUBSURFACE_FLUXES_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME CONSTITUTIVE_RELATIONS_GENERIC_EVALUATORS_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_RELATIONS_INCLUDES)
include_evaluators_directories(LISTNAME WRM_MODELS_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME WMR_EVALUATORS_FACTORIES_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_ICY_OVERLAND_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_OVERLAND_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_OVERLAND_HEAD_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME FLOW_OVERLAND_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME RICHARDS_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME PERMAFROST_EVAL_FAC_INCLUDES)
include_evaluators_directories(LISTNAME DEFORMATION_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_ADVECTION_DIFFUSION_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_CONSTANT_TEMPERATURE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_RELATIONS_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_SURFACE_ICE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_TWO_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME ENERGY_THREE_PHASE_REG_INCLUDES)
include_evaluators_directories(LISTNAME SURFACE_BALANCE_SEB_REG_INCLUDES)
include_evaluators_directories(LISTNAME BGC_REG_INCLUDES)
include_evaluators_directories(LISTNAME TEST_PKS_REG_INCLUDES)
include_evaluators_directories(LISTNAME MPC_PK_FAC_INCLUDES)
#include_evaluators_directories(LISTNAME TRANSPORT_PASSIVE_TRACER_REG_INCLUDES)
  


add_executable(ats ats_mesh_factory.cc simulation_driver.cc main.cc )

set(ATS_LIBS
		     coordinator
                     mpc
                     #  pk_test_pks
                      # pk_energy_constant_temperature
                      pk_energy_advection_diffusion
		      pk_energy_two_phase
                      pk_energy_three_phase
                      pk_energy_surface_ice
                      pk_energy_base
		      pk_flow_permafrost
                      pk_flow_richards
                      pk_flow_snow_distribution
                      pk_flow_overland
                      pk_flow_overland_pressure
                      pk_flow_icy_overland
                      pk_transport
                      pk_coupled_transport
	              # pk_prescribed_deformation
#                      pk_prescribed_volumetric_deformation
                      pk_volumetric_deformation
                      pk_surface_balance_SEB
                      pk_BGC
                      pk_bases
                      #deform_relations_porosity
                      energy_relations
                      flow_relations
                      flow_relations_surface_subsurface_fluxes
                      generic_evaluators
                      relations_eos
                      relations_ewc
                      bc_factory
                      advection
                      divgrad
                      deformation_operator
)

set(AMANZI_LIBS
                      amanzi_operators
                      amanzi_pks
                      amanzi_state
	              amanzi_whetstone
                      amanzi_time_integration
	              amanzi_solvers
                      amanzi_data_structures
                      amanzi_mesh
                      amanzi_mesh_audit
                      amanzi_mesh_functions
                      amanzi_functions
                      amanzi_geometry
                      amanzi_mesh_factory
                      amanzi_interface
                      amanzi_output
                      amanzi_mstk_mesh
		      amanzi_mesh_logical
                      amanzi_geochemrxns
		      amanzi_geochembase
		      amanzi_geochemutil
		      amanzi_geochemsolvers
		      amanzi_chemistry_pk
		      amanzi_mpc_tree
		      amanzi_transport
                      amanzi_dbg
)

# set(GEOCHEM_LIBS 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libchemistry_pk.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libgeochembase.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libgeochemsolvers.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libgeochemutil.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libmpc_tree.a 
#                  /home/scratch/dasvyat/Amanzi/build/amanzi_default/lib/libtransport.a
# )

target_link_libraries(ats ${ATS_LIBS} ${AMANZI_LIBS})# ${GEOCHEM_LIBS})


install(TARGETS ats DESTINATION bin)

#------------------------------------------------------------------------------#
# ATS F90 test program
#------------------------------------------------------------------------------#

if (Amanzi_DISABLE_MPI_CXX_BINDINGS)  

   add_executable(ats_f90 main.f90)

   target_link_libraries(
        # Target
	ats_f90

        # ATS cross lib
	atsx

    ${ATS_LIBS}
	${AMANZI_LIBS}
	${Amanzi_TPL_Trilinos_LIBRARIES}

        # MPI
	${Amanzi_MPI_Fortran_LIBRARIES}
	${Amanzi_MPI_C_LIBRARIES}
   )

   set_target_properties(ats_f90 PROPERTIES LINKER_LANGUAGE Fortran)

   #---------------------------------------------------------------------------#
   # Install rule
   #---------------------------------------------------------------------------#

   install(TARGETS ats_f90 DESTINATION bin)



   add_executable(ats_test_clm_driver test_clm_driver.cc)

   target_link_libraries(
        # Target
	ats_test_clm_driver

        # ATS cross lib
	atsx

    ${ATS_LIBS}
	${AMANZI_LIBS}
	${Amanzi_TPL_Trilinos_LIBRARIES}

        # MPI
	${Amanzi_MPI_C_LIBRARIES}
   )

   #---------------------------------------------------------------------------#
   # Install rule
   #---------------------------------------------------------------------------#

   install(TARGETS ats_test_clm_driver DESTINATION bin)

endif()
