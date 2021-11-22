#
# ATS
#
# NOTE: ATS is not standalone code, and this is not a stand-alone
#       CMakeLists.txt.  Instead, it must be built as a subrepo of
#       Amanzi.  See ATS/INSTALL.md

message(STATUS "\n-------------------------------------------")
message(STATUS "\n-- CMake: Configuring ATS build/install.\n--")
message(STATUS "----------------------------------------")
project(ATS)
add_subdirectory(src)
#add_subdirectory(testing)

# ATS Custom CMake Modules and Macros
set(ATS_MODULE_PATH "${AMANZI_SOURCE_DIR}/src/physics/ats/tools/cmake")
list(APPEND CMAKE_MODULE_PATH "${ATS_MODULE_PATH}")

# ATS code version defined
include(ATSVersion)
