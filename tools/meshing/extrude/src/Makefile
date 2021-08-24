ifeq (${ATS_BUILD_TYPE},Debug)
	CXX_FLAGS = -g -O0 -DENABLE_DBC
else
	CXX_FLAGS = -g -O3
endif

CXX_FLAGS +=  -std=c++11

TPLS_LIB = ${AMANZI_TPLS_DIR}/lib
TPLS_INCLUDE = ${AMANZI_TPLS_DIR}/include
TPLS_LIBS = -lexodus -lnetcdf -lhdf5 -lhdf5_hl -lz


SRCS = dbc.cc exceptions.cc Mesh2D.cc Mesh3D.cc readMesh2D.cc writeMesh3D.cc

OBJS=$(SRCS:%.cc=%.o)
DEPS=$(OBJS:%.o=%.d)

extrude.a: ${DEPS} ${OBJS}
	ar rcs extrude.a ${OBJS}

clean:
	rm -f ./*.o
	rm -f extrude.a
	rm -f .depend
	rm -f ./*.d

%.d: %.cc
	mpicxx -std=c++11 $(CXX_FLAGS) -I $(TPLS_INCLUDE) -MM -MF $(patsubst %.o,%.d,$@) -o $@ $<

%.o: %.cc
	mpicxx -std=c++11 $(CXX_FLAGS) -I $(TPLS_INCLUDE) -c $< -o $@

-include $(DEPS)



