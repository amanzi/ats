#!/usr/bin/env bash

# NOTE: see INSTALL.md before using this!

if [ ${ATS_BUILD_TYPE} == Debug ]
then
    dbg_option=--debug
elif [ ${ATS_BUILD_TYPE} == RelWithDebInfo ]
then
    dbg_option=--relwithdebinfo
else
    dbg_option=--opt
fi

echo "Running Amanzi Boostrap with:"
echo "  AMANZI_SRC_DIR: ${AMANZI_SRC_DIR}"
echo "  AMANZI_BUILD_DIR: ${AMANZI_BUILD_DIR}"
echo "  AMANZI_INSTALL_DIR: ${AMANZI_DIR}"
echo "  AMANZI_TPLS_BUILD_DIR: ${AMANZI_TPLS_BUILD_DIR}"
echo "  AMANZI_TPLS_DIR: ${AMANZI_TPLS_DIR}"
echo " with build type: ${ATS_BUILD_TYPE}"
echo ""


${AMANZI_SRC_DIR}/bootstrap.sh \
   ${dbg_option} \
   --with-mpi=${OPENMPI_DIR} \
   --enable-shared \
   --disable-structured  --enable-unstructured \
   --disable-stk_mesh --enable-mstk_mesh \
   --enable-hypre \
   --enable-silo \
   --disable-petsc \
   --disable-physics \
   --enable-ats \
   --amanzi-install-prefix=${AMANZI_DIR} \
   --amanzi-build-dir=${AMANZI_BUILD_DIR} \
   --tpl-install-prefix=${AMANZI_TPLS_DIR} \
   --tpl-build-dir=${AMANZI_TPLS_BUILD_DIR} \
   --tpl-download-dir=${ATS_BASE}/amanzi-tpls/Downloads \
   --tools-download-dir=${ATS_BASE}/amanzi-tpls/Downloads \
   --tools-build-dir=${ATS_BASE}/build \
   --tools-install-prefix=${ATS_BASE}/install \
   --with-cmake=`which cmake` \
   --with-ctest=`which ctest`
