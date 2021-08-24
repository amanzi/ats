ifeq (${ATS_BUILD_TYPE},Debug)
	CXX_FLAGS = -g -O0 -DENABLE_DBC
else
	CXX_FLAGS = -g -O3
endif

CXX_FLAGS +=  -std=c++11

TPLS_LIB = ${AMANZI_TPLS_DIR}/lib
TPLS_INCLUDE = ${AMANZI_TPLS_DIR}/include
TPLS_LIBS = -lexodus -lnetcdf -lhdf5_hl -lhdf5 -lz

all: extrude_one_layer extrude_uniform extrude_homogeneous_uniform extrude_variable extrude_homogeneous_variable extrude_homogeneous_variable_with0

extrude_one_layer: extrude_one_layer.o extrude.a
	mpicxx -std=c++11 ${CXX_FLAGS} extrude_one_layer.o src/extrude.a -L${TPLS_LIB} ${TPLS_LIBS} -o extrude_one_layer

extrude_uniform: extrude_uniform.o extrude.a
	mpicxx -std=c++11 ${CXX_FLAGS} extrude_uniform.o src/extrude.a -L${TPLS_LIB} ${TPLS_LIBS} -o extrude_uniform

extrude_homogeneous_uniform: extrude_homogeneous_uniform.o extrude.a
	mpicxx -std=c++11 ${CXX_FLAGS} extrude_homogeneous_uniform.o src/extrude.a -L${TPLS_LIB} ${TPLS_LIBS} -o extrude_homogeneous_uniform

extrude_homogeneous_variable: extrude_homogeneous_variable.o extrude.a
	mpicxx -std=c++11 ${CXX_FLAGS} extrude_homogeneous_variable.o src/extrude.a -L${TPLS_LIB} ${TPLS_LIBS} -o extrude_homogeneous_variable

extrude_homogeneous_variable_with0: extrude_homogeneous_variable_with0.o extrude.a
	mpicxx -std=c++11 ${CXX_FLAGS} extrude_homogeneous_variable_with0.o src/extrude.a -L${TPLS_LIB} ${TPLS_LIBS} -o extrude_homogeneous_variable_with0

extrude_variable: extrude.a extrude_variable.o
	mpicxx -std=c++11 ${CXX_FLAGS} extrude_variable.o src/extrude.a -L${TPLS_LIB} ${TPLS_LIBS} -o extrude_variable

extrude.a:
	make -C src extrude.a

clean:
	rm -f ./*.o
	rm -f extrude_homogeneous_uniform extrude_homogeneous_variable extrude_one_layer extrude_uniform extrude_variable extrude_homogeneous_variable_with0
	rm -f .depend
	rm -f ./*.d
	make -C src clean
%.o: %.cc
	mpicxx -std=c++11 $(CXX_FLAGS) -I src -I $(TPLS_INCLUDE) -c $< -o $@




