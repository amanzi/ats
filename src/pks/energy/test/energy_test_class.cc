/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "energy_test_class.hh"

EnergyTest::EnergyTest(Teuchos::ParameterList& plist_,
                       const Teuchos::RCP<AmanziMesh::Mesh>& mesh_,
                       int num_components_)
  : mesh(mesh_), getParameterList(plist_), num_components(num_components_)
{
  // create states
  Teuchos::ParameterList state_plist = parameter_list.get<Teuchos::ParameterList>("State");
  S0 = Teuchos::rcp(new State(state_plist));
  S0->RegisterDomainMesh(mesh);
  Teuchos::ParameterList energy_plist = parameter_list.get<Teuchos::ParameterList>("Energy");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector("solution"));

  // create the PK
  EPK = Teuchos::rcp(new Energy::TwoPhase(energy_plist, S0, soln));
}

void
EnergyTest::initialize()
{
  // initialize state, including darcy flux, sat, density, poro from parameter list
  S0->Initialize();
  S0->set_time(0.0);

  // initialize energy
  EPK->initialize(S0);
  initialize_owned();
  initialize_water_flux();

  // finish checking state and create the state at the new timestep
  if (!S0->CheckAllInitialized()) { std::cout << "DID NOT INITIALIZE THINGS!" << std::endl; }
  S1 = Teuchos::rcp(new State(*S0));
  *S1 = *S0;
  EPK->set_states(S0, S0, S1);
}

void
EnergyTest::commit_step()
{
  *S0 = *S1;
}

void
EnergyTest::initialize_owned()
{
  Teuchos::RCP<CompositeVector> temp = S0->GetFieldData("temperature", "energy");

  int c_owned = temp->size("cell");
  for (int c = 0; c != c_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    for (int lcv_comp = 0; lcv_comp != num_components; ++lcv_comp) {
      (*temp)("cell", lcv_comp, c) = my_f(xc, 0.0);
    }
  }
  S0->GetField("temperature", "energy")->set_initialized();
}

void
EnergyTest::initialize_water_flux()
{
  const Epetra_BlockMap& fmap = mesh->getMap(AmanziMesh::Entity_kind::FACE,true);
  Teuchos::RCP<CompositeVector> water_flux = S0->GetFieldData("water_flux", "state");

  for (int f = fmap.MinLID(); f <= fmap.MaxLID(); f++) {
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
    const AmanziGeometry::Point& fc = mesh->getFaceCentroid(f);
    (*water_flux)("face", f) = my_u(fc, 0.0) * normal;
  }
  S0->GetField("water_flux", "state")->set_initialized();
}

void
EnergyTest::evaluate_error_temp(double t, double* L1, double* L2)
{
  const Epetra_BlockMap& cmap = mesh->getMap(AmanziMesh::Entity_kind::CELL,true);
  Teuchos::RCP<const CompositeVector> temp = S1->GetFieldData("temperature");

  double d;
  *L1 = 0.0;
  *L2 = 0.0;
  for (int c = cmap.MinLID(); c <= cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    double volume = mesh->getCellVolume(c);

    for (int lcv_comp = 0; lcv_comp != num_components; ++lcv_comp) {
      d = (*temp)("cell", lcv_comp, c) - my_f(xc, t);

      *L1 += fabs(d) * volume;
      *L2 += d * d * volume;
    }
  }
  *L2 = sqrt(*L2);
}

AmanziGeometry::Point
EnergyTest::my_u(const AmanziGeometry::Point& x, double t)
{
  return AmanziGeometry::Point(1.0, 0.0, 0.0);
}

// test problem with constant solution of 1
EnergyTestOne::EnergyTestOne(Teuchos::ParameterList& plist_,
                             const Teuchos::RCP<AmanziMesh::Mesh>& mesh_,
                             int num_components_)
  : EnergyTest::EnergyTest(plist_, mesh_, num_components_)
{}

double
EnergyTestOne::my_f(const AmanziGeometry::Point& x, double t)
{
  return 293.15;
}

// test problem with step solution
EnergyTestStep::EnergyTestStep(Teuchos::ParameterList& plist_,
                               const Teuchos::RCP<AmanziMesh::Mesh>& mesh_,
                               int num_components_)
  : EnergyTest::EnergyTest(plist_, mesh_, num_components_)
{}

double
EnergyTestStep::my_f(const AmanziGeometry::Point& x, double t)
{
  if (x[0] <= 1 + t) return 293.15;
  return 0;
}

// // test problem with smooth solution
// EnergyTestSmooth::EnergyTestSmooth(Teuchos::ParameterList& plist_,
//         const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
//   EnergyTest::EnergyTest(plist_, mesh_, num_components_) {}

// double EnergyTestSmooth::my_f(const AmanziGeometry::Point& x, double t) {
//   return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
// }

// // test problem with cubic solution
// EnergyTestCubic::EnergyTestCubic(Teuchos::ParameterList& plist_,
//         const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
//   EnergyTest::EnergyTest(plist_, mesh_, num_components_) {}

// double EnergyTestCubic::my_f(const AmanziGeometry::Point& x, double t) {
//   if( x[0] < 1 + t ) return 1;
//   if( x[0] > 3 + t ) return 0;
//   double z = (x[0]-1-t) / 2;
//   return 2*z*z*z - 3*z*z + 1;
// }

// test problem on 2D square with velocity in the <1,1> x-y direction, still virtual
EnergyTestTwoDOne::EnergyTestTwoDOne(Teuchos::ParameterList& plist_,
                                     const Teuchos::RCP<AmanziMesh::Mesh>& mesh_,
                                     int num_components_)
  : EnergyTest::EnergyTest(plist_, mesh_, num_components_)
{}

AmanziGeometry::Point
EnergyTestTwoDOne::my_u(const AmanziGeometry::Point& x, double t)
{
  return AmanziGeometry::Point(1.0, 1.0);
}

double
EnergyTestTwoDOne::my_f(const AmanziGeometry::Point& x, double t)
{
  return 293.15;
}
