# -*- mode: cmake -*-

#
# ATS Version Information:
# 
# Information about the current source is extracted from the git repository and used to 
# create the version string (ATS_VERSION).  
#
# NOTE: this information won't be accessible without the full repository.
#       So for releases we need to extract this and set it as part of the tarball creation.
#
#   * if ats_version.hh does not exist create it
#       * if git is found
#            use git to create version strings 
#       * else
#            use statically defined version strings
#       * endif
#   * endif
#   install ats_version.hh
#

include(PrintVariable)
include(InstallManager)

message(STATUS "")
message(STATUS ">>>>>>>> ATSVersion.cmake")

message(STATUS ">>>> JDM: ${CMAKE_SOURCE_DIR}")
set(ATS_SUBMODULE_DIR ${CMAKE_SOURCE_DIR}/src/physics/ats)
message(STATUS ">>>> JDM: ${ATS_SUBMODULE_DIR}")

find_package(Git)

if ( (EXISTS ${ATS_SUBMODULE_DIR}/.git) AND (GIT_FOUND) ) 

  # Get the name of the current branch.
  set(GIT_ARGS status)
  execute_process(COMMAND ${GIT_EXECUTABLE} ${GIT_ARGS}
	          WORKING_DIRECTORY ${ATS_SUBMODULE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE ATS_GIT_STATUS
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if (err_occurred)
    message(WARNING "Error executing git:\n ${cmd}\n${err}")
    set(cmd_output cmd_output-NOTFOUND)
    exit()
  endif()


  # too noisy command...
  # message(STATUS ">>>> JDM: ATS_GIT_STATUS:      ${ATS_GIT_STATUS}")


  # Put the status in a list
  STRING(REPLACE "\n" ";" ATS_GIT_STATUS_LIST ${ATS_GIT_STATUS})
  # Extract the first entry - reuse the ATS_GIT_STATUS variable
  LIST(GET ATS_GIT_STATUS_LIST 0 ATS_GIT_STATUS)
  if (${ATS_GIT_STATUS} MATCHES "(D|d)etached") 
    # For now just set branch to detached - we could add a lookup for tags later
    set(ATS_GIT_BRANCH detached)
  elseif(${ATS_GIT_STATUS} MATCHES "On branch")
    # Extract the branch name
    STRING(REPLACE "On branch " "" ATS_GIT_BRANCH ${ATS_GIT_STATUS})
  endif()

  message(STATUS ">>>> JDM: ATS_GIT_BRANCH = ${ATS_GIT_BRANCH}")

  # Extract the lastest tag of the form ats-*

  # Get the hash of the current version
  set(GIT_ARGS rev-parse --short HEAD)
  execute_process(COMMAND ${GIT_EXECUTABLE} ${GIT_ARGS}
	          WORKING_DIRECTORY ${ATS_SUBMODULE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE ATS_GIT_GLOBAL_HASH
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if(err_occurred)
    message(WARNING "Error executing git:\n ${cmd}\n${err}")
    set(cmd_output cmd_output-NOTFOUND)
    exit()
  endif()

  message(STATUS ">>>> JDM: ATS_GIT_GLOBAL_HASH: ${ATS_GIT_GLOBAL_HASH}")

  # Get the latest ats-* version number tag
  set(GIT_ARGS tag -l ats-*)
  execute_process(COMMAND  ${GIT_EXECUTABLE} ${GIT_ARGS}
	          WORKING_DIRECTORY ${ATS_SUBMODULE_DIR}
                  RESULT_VARIABLE err_occurred 
                  OUTPUT_VARIABLE ATS_GIT_LATEST_TAG
                  ERROR_VARIABLE err
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)

   # Put the tags in a list
   STRING(REPLACE "\n" ";" ATS_GIT_LATEST_TAG_LIST ${ATS_GIT_LATEST_TAG})
   # Extract the lastest tag of the form ats-*
   IF ( ${ATS_GIT_BRANCH} MATCHES "master" ) 
     FOREACH(atag ${ATS_GIT_LATEST_TAG_LIST})
       IF ( ${atag} MATCHES "^ats-.*-dev" )
         set ( ATS_GIT_LATEST_TAG ${atag} )
       ENDIF()
     ENDFOREACH()
   ELSE()
     FOREACH(atag ${ATS_GIT_LATEST_TAG_LIST})
       IF ( ${atag} MATCHES "^ats-[0-9]\\.[0-9][0-9]?\\.[0-9]" )
         set ( ATS_GIT_LATEST_TAG ${atag} )
       ENDIF()
     ENDFOREACH()
   ENDIF()

   # message(STATUS ">>>> JDM: GIT_EXEC        = ${GIT_EXECUTABLE}")
   # message(STATUS ">>>> JDM: GIT_ARGS        = ${GIT_ARGS}")
   # message(STATUS ">>>> JDM: RESULT_VARIABLE = ${err_occurred}")
   # message(STATUS ">>>> JDM: ATS_GIT_LATEST_TAG = ${ATS_GIT_LATEST_TAG}")

   STRING(REGEX REPLACE "ats-" "" ATS_GIT_LATEST_TAG_VER ${ATS_GIT_LATEST_TAG})	
   STRING(REGEX REPLACE "\\..*" "" ATS_GIT_LATEST_TAG_MAJOR ${ATS_GIT_LATEST_TAG_VER})	
   STRING(REGEX MATCH "\\.[0-9][0-9]?[\\.,-]" ATS_GIT_LATEST_TAG_MINOR ${ATS_GIT_LATEST_TAG_VER})  	
   STRING(REGEX REPLACE "[\\.,-]" "" ATS_GIT_LATEST_TAG_MINOR ${ATS_GIT_LATEST_TAG_MINOR} )	

   set(ATS_VERSION_MAJOR ${ATS_GIT_LATEST_TAG_MAJOR})
   set(ATS_VERSION_MINOR ${ATS_GIT_LATEST_TAG_MINOR})

   #
   # ATS version
   #
   set(ATS_VERSION ${ATS_GIT_LATEST_TAG_VER}_${ATS_GIT_GLOBAL_HASH})

   STRING(REGEX REPLACE ".*\\.[0-9][0-9]?[\\.,-]" "" ATS_VERSION_PATCH ${ATS_VERSION})
   STRING(REGEX REPLACE ".*_" "" ATS_VERSION_HASH ${ATS_VERSION_PATCH})
   STRING(REGEX REPLACE "_.*" "" ATS_VERSION_PATCH ${ATS_VERSION_PATCH})

else()

  message(STATUS "  >>>>>>>> Using static version information to create ats_version.hh")
  if ( NOT GIT_FOUND ) 
    message(STATUS "    >>>>>> Could not locate Git executable.")
  endif()
  if ( NOT EXISTS ${CMAKE_SOURCE_DIR}/.git/ )
    message(STATUS "    >>>>>> Release or snapshot, no .git directory found.")
  endif()

  #
  # For releases and snapshots, set static information before creating the tarball.
  #
  set(ATS_GIT_BRANCH master )
  set(ATS_GIT_GLOBAL_HASH )

  set(ATS_VERSION_MAJOR 1)
  set(ATS_VERSION_MINOR 0)
  set(ATS_VERSION_PATCH 0)
  set(ATS_VERSION_HASH ${ATS_GIT_GLOBAL_HASH})

  #
  # ATS version
  #
  set(ATS_VERSION ${ATS_VERSION_MAJOR}.${ATS_VERSION_MINOR}-${ATS_VERSION_PATCH}_${ATS_VERSION_HASH})

endif()

# Write the version header file
                              set(version_template ${ATS_SUBMODULE_DIR}/tools/cmake/ats_version.hh.in)
configure_file(${version_template}
               ${CMAKE_CURRENT_BINARY_DIR}/ats_version.hh
               @ONLY)
configure_file(${version_template}
               ${CMAKE_CURRENT_BINARY_DIR}/extras/ats_version.hh
               @ONLY)

add_install_include_file(${CMAKE_CURRENT_BINARY_DIR}/ats_version.hh)             

message(STATUS "\t >>>>>  ATS Version: ${ATS_VERSION}")
message(STATUS "\t >>>>>  MAJOR ${ATS_VERSION_MAJOR}")
message(STATUS "\t >>>>>  MINOR ${ATS_VERSION_MINOR}")
message(STATUS "\t >>>>>  PATCH ${ATS_VERSION_PATCH}")
message(STATUS "\t >>>>>  HASH  ${ATS_VERSION_HASH}")

