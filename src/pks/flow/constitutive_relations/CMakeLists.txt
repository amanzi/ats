# -*- mode: cmake -*-

#
#  ATS
#    Constitutive relations for flow
#


include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
add_definitions("-DMSTK_HAVE_MPI")

list(APPEND subdirs wrm porosity overland_conductivity elevation water_content)

set(sources "")

# include all subdirectories
foreach(lcv IN LISTS subdirs)
  include_directories($lcv)
endforeach(lcv)

# collect all registrations and register
file(GLOB_RECURSE registrations "./*_reg.hh" )
message("My reg: ${registrations}")
foreach(reg_lcv IN LISTS registrations)
  message("REGISTERING: ${reg_lcv}")
  register_abs_evaluator_with_factory(HEADERFILE ${reg_lcv} LISTNAME FLOW_RELATIONS)
endforeach(reg_lcv)

generate_evaluators_registration_header(
  HEADERFILE flow_relations_registration.hh
  LISTNAME   FLOW_RELATIONS
  INSTALL    True
  )


# collect all sources
foreach(lcv IN LISTS subdirs)

file(GLOB subdir_sources "./${lcv}/*.cc")
  set(sources ${sources} ${subdir_sources})
endforeach(lcv)

# make the library
add_library(flow_relations ${sources})

    
