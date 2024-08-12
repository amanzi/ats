/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   Amanzi

   License: see $AMANZI_DIR/COPYRIGHT
   Author: Svetlana Tokareva (tokareva@lanl.gov)

   Helper class for lake unit tests
   ------------------------------------------------------------------------- */

#ifndef LAKE_THERMO_TEST_CLASS_HH_
#define LAKE_THERMO_TEST_CLASS_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "State.hh"
#include "lake_thermo_pk.hh"

using namespace Amanzi;

class LakeThermoTest {

public:
  // constuctor
  LakeThermoTest(Teuchos::ParameterList& plist,
                const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                int num_components);

  // init
  void initialize();

  // data
  Teuchos::ParameterList parameter_list;
  Teuchos::RCP<State> S0, S1;
  Teuchos::RCP<Lake_Thermo_PK> LTPK;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  int num_components;
  double constant_value;

  // helper analytic solutions
  virtual double my_f(const AmanziGeometry::Point& x, double t) = 0;
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t);

  // helper evaluation methods
  void commit_step();
  void initialize_owned();
  void initialize_mass_flux();
  void evaluate_error_temp(double t, double* L1, double* L2);
};

class LakeThermoTestTestOne : public LakeThermoTest {
public:
  LakeThermoTestTestOne(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class LakeThermoTestStep : public LakeThermoTest {
public:
  LakeThermoTestStep(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class LakeThermoTestSmooth : public LakeThermoTest {
public:
  LakeThermoTestSmooth(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class LakeThermoTestCubic : public LakeThermoTest {
public:
  LakeThermoTestCubic(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class LakeThermoTestTwoDOne : public LakeThermoTest {
public:
  LakeThermoTestTwoDOne(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t);
};

#endif

