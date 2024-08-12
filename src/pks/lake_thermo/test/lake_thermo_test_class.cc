#include "lake_thermo_test_class.hh"

LakeThermoTest::LakeThermoTest(Teuchos::ParameterList& plist_,
                             const Teuchos::RCP<AmanziMesh::Mesh>& mesh_,
                             int num_components_) :
    mesh(mesh_), parameter_list(plist_), num_components(num_components_) {

  // create states
  Teuchos::ParameterList state_plist =
    parameter_list.get<Teuchos::ParameterList>("State");
  S0 = Teuchos::rcp(new State(state_plist));
  S0->RegisterDomainMesh(mesh);
  Teuchos::ParameterList energy_plist =
    parameter_list.get<Teuchos::ParameterList>("Energy");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector("solution"));

  // create the PK
  LTPK = Teuchos::rcp(new Lake_Thermo_PK(energy_plist, S0, soln));
}

void LakeThermoTest::initialize() {
  // initialize state, including darcy flux, sat, density, poro from parameter list
  S0->Initialize();
  S0->set_time(0.0);

  // initialize energy
  LTPK->initialize(S0);
  initialize_owned();
  initialize_mass_flux();

  // finish checking state and create the state at the new timestep
  if (!S0->CheckAllInitialized()) {
    std::cout << "DID NOT INITIALIZE THINGS!" << std::endl;
  }
  S1 = Teuchos::rcp(new State(*S0));
  *S1 = *S0;
  LTPK->set_states(S0,S0,S1);
}

void LakeThermoTest::commit_step() {
  *S0 = *S1;
}

void LakeThermoTest::initialize_owned() {
  Teuchos::RCP<CompositeVector> temp = S0->GetFieldData("temperature", "energy");

  int c_owned = temp->size("cell");
  for (int c=0; c != c_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    for (int lcv_comp=0; lcv_comp != num_components; ++lcv_comp) {
      (*temp)("cell",lcv_comp,c) = my_f(xc, 0.0);
    }
  }
  S0->GetField("temperature", "energy")->set_initialized();
}

void LakeThermoTest::initialize_mass_flux() {
  const Epetra_BlockMap& fmap = mesh->face_map(true);
  Teuchos::RCP<CompositeVector> mass_flux =
    S0->GetFieldData("mass_flux", "state");

  for (int f=fmap.MinLID(); f<=fmap.MaxLID(); f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    const AmanziGeometry::Point& fc = mesh->face_centroid(f);
    (*mass_flux)("face",f) = my_u(fc, 0.0) * normal;
  }
  S0->GetField("mass_flux", "state")->set_initialized();
}

void LakeThermoTest::evaluate_error_temp(double t, double* L1, double* L2) {
  const Epetra_BlockMap& cmap = mesh->cell_map(true);
  Teuchos::RCP<const CompositeVector> temp = S1->GetFieldData("temperature");

  double d;
  *L1 = 0.0;
  *L2 = 0.0;
  for (int c=cmap.MinLID(); c<=cmap.MaxLID(); c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    for (int lcv_comp=0; lcv_comp != num_components; ++lcv_comp) {
      d = (*temp)("cell",lcv_comp,c) - my_f(xc, t);

      *L1 += fabs(d) * volume;
      *L2 += d * d * volume;
    }
  }
  *L2 = sqrt(*L2);
}

AmanziGeometry::Point LakeThermoTest::my_u(const AmanziGeometry::Point& x, double t) {
  return AmanziGeometry::Point(1.0, 0.0, 0.0);
}

// test problem with constant solution of 1
LakeThermoTestOne::LakeThermoTestOne(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  LakeThermoTest::LakeThermoTest(plist_, mesh_, num_components_) {}

double LakeThermoTestOne::my_f(const AmanziGeometry::Point& x, double t) {
  return 293.15;
}

// test problem with step solution
LakeThermoTestStep::LakeThermoTestStep(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  LakeThermoTest::LakeThermoTest(plist_, mesh_, num_components_) {}

double LakeThermoTestStep::my_f(const AmanziGeometry::Point& x, double t) {
  if (x[0] <= 1 + t) return 293.15;
  return 0;
}

// // test problem with smooth solution
// LakeThermoTestSmooth::LakeThermoTestSmooth(Teuchos::ParameterList& plist_,
//         const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
//   LakeThermoTest::LakeThermoTest(plist_, mesh_, num_components_) {}

// double LakeThermoTestSmooth::my_f(const AmanziGeometry::Point& x, double t) {
//   return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
// }

// // test problem with cubic solution
// LakeThermoTestCubic::LakeThermoTestCubic(Teuchos::ParameterList& plist_,
//         const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
//   LakeThermoTest::LakeThermoTest(plist_, mesh_, num_components_) {}

// double LakeThermoTestCubic::my_f(const AmanziGeometry::Point& x, double t) {
//   if( x[0] < 1 + t ) return 1;
//   if( x[0] > 3 + t ) return 0;
//   double z = (x[0]-1-t) / 2;
//   return 2*z*z*z - 3*z*z + 1;
// }

// test problem on 2D square with velocity in the <1,1> x-y direction, still virtual
LakeThermoTestTwoDOne::LakeThermoTestTwoDOne(Teuchos::ParameterList& plist_,
        const Teuchos::RCP<AmanziMesh::Mesh>& mesh_, int num_components_) :
  LakeThermoTest::LakeThermoTest(plist_, mesh_, num_components_) {}

AmanziGeometry::Point LakeThermoTestTwoDOne::my_u(const AmanziGeometry::Point& x, double t) {
  return AmanziGeometry::Point(1.0, 1.0);
}

double LakeThermoTestTwoDOne::my_f(const AmanziGeometry::Point& x, double t) {
  return 293.15;
}


