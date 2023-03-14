/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   Amanzi

   Helper class for energy unit tests
   ------------------------------------------------------------------------- */

#ifndef ENERGY_TEST_CLASS_HH_
#define ENERGY_TEST_CLASS_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "State.hh"
#include "two_phase.hh"

using namespace Amanzi;

class EnergyTest {
 public:
  // constuctor
  EnergyTest(Teuchos::ParameterList& plist,
             const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
             int num_components);

  // init
  void initialize();

  // data
  Teuchos::ParameterList parameter_list;
  Teuchos::RCP<State> S0, S1;
  Teuchos::RCP<Energy::TwoPhase> EPK;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  int num_components;
  double constant_value;

  // helper analytic solutions
  virtual double my_f(const AmanziGeometry::Point& x, double t) = 0;
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t);

  // helper evaluation methods
  void commit_step();
  void initialize_owned();
  void initialize_water_flux();
  void evaluate_error_temp(double t, double* L1, double* L2);
};

class EnergyTestOne : public EnergyTest {
 public:
  EnergyTestOne(Teuchos::ParameterList& plist,
                const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class EnergyTestStep : public EnergyTest {
 public:
  EnergyTestStep(Teuchos::ParameterList& plist,
                 const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                 int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class EnergyTestSmooth : public EnergyTest {
 public:
  EnergyTestSmooth(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                   int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class EnergyTestCubic : public EnergyTest {
 public:
  EnergyTestCubic(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                  int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
};

class EnergyTestTwoDOne : public EnergyTest {
 public:
  EnergyTestTwoDOne(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                    int num_components);

  virtual double my_f(const AmanziGeometry::Point& x, double t);
  virtual AmanziGeometry::Point my_u(const AmanziGeometry::Point& x, double t);
};

#endif
