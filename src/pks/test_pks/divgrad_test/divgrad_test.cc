/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

/* -------------------------------------------------------------------------
  A high level test class for the MatrixMFD operator.

------------------------------------------------------------------------- */

#include "test_pk_bc_factory.hh"

#include "divgrad_test.hh"

namespace Amanzi {
namespace TestPKs {

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void
DivGradTest::setup(const Teuchos::Ptr<State>& S)
{
  PKPhysicalBase::setup(S);

  // Require fields and evaluators for those fields.
  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2, 1);
  locations2[0] = AmanziMesh::Entity_kind::CELL;
  locations2[1] = AmanziMesh::Entity_kind::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->Require<CompositeVector, CompositeVectorSpace>(key_, Tags::NEXT, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponents(names2, locations2, num_dofs2);

  // Get data for non-field quanitites.
  S->RequireEvaluator("cell_volume");

  // Create the absolute permeability tensor.
  int c_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>(c_owned));
  for (int c = 0; c != c_owned; ++c) {
    (*K)[c].Init(mesh_->getSpaceDimension(), 1);
    (*K)[c](0, 0) = 1.0;
  }

  // Create the boundary condition data structures.
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  TestPKBCFactory bc_factory(mesh_, bc_plist);
  bc_dirichlet_ = bc_factory.CreateDirichlet();
  bc_neumann_ = bc_factory.CreateNeumann();

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_->sublist("diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, mesh_));
  matrix_->set_symmetric(true);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->InitPreconditioner();
  matrix_->CreateMFDmassMatrices(K.ptr());
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void
DivGradTest::initialize(const Teuchos::Ptr<State>& S)
{
  // Check for PK-specific initialization
  if (!plist_->isSublist("initial condition")) {
    std::stringstream messagestream;
    messagestream << name_ << " has no initial condition parameter list.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // make sure the initial condition doesn't set faces in another way
  Teuchos::ParameterList ic_plist = plist_->sublist("initial condition");
  ic_plist.set("initialize faces from cells", false);

  // initialize primary variable from the ic_plist condition
  PKPhysicalBase::initialize(S);

  S->Get<CompositeVector>(key_).ScatterMasterToGhosted("face");

  // initialize boundary conditions
  int nfaces = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  bc_markers_.resize(nfaces, Operators::OPERATOR_BC_NONE);
  bc_values_.resize(nfaces, 0.0);

  // assemble and set BCs
  bc_dirichlet_->Compute(0.0);
  bc_neumann_->Compute(0.0);
  UpdateBoundaryConditions_();

  // assemble matrix
  matrix_->CreateMFDstiffnessMatrices(Teuchos::null);
  matrix_->CreateMFDrhsVectors();
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  // derive consistent faces
  Teuchos::RCP<CompositeVector> soln = S->GetPtrW<CompositeVector>(key_, name_);

  bool fail = TestRegularFaceValues_(soln);
  if (fail) {
    std::cout << "Failed test, which is good" << std::endl;
  } else {
    std::cout << "Passed test, which is BAD!" << std::endl;
  }
  AMANZI_ASSERT(fail);

  matrix_->UpdateConsistentFaceConstraints(soln.ptr());

  // test for correctness -- this should only work for regular meshes
  fail = TestRegularFaceValues_(soln);
  if (fail) {
    std::cout << "Failed test, which is BAD!" << std::endl;
  } else {
    std::cout << "Passed test, which is good" << std::endl;
  }
  AMANZI_ASSERT(!fail);
};


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void
DivGradTest::UpdateBoundaryConditions_()
{
  for (unsigned int n = 0; n != bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::OPERATOR_BC_NONE;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  for (bc = bc_dirichlet_->begin(); bc != bc_dirichlet_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  for (bc = bc_neumann_->begin(); bc != bc_neumann_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_values_[f] = bc->second;
  }
};


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void
DivGradTest::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& pres)
{
  Epetra_MultiVector& pres_f = *pres->ViewComponent("face", true);
  int nfaces = pres->size("face");
  for (int f = 0; f != nfaces; ++f) {
    if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) { pres_f[0][f] = bc_values_[f]; }
  }
};


bool
DivGradTest::TestRegularFaceValues_(const Teuchos::RCP<CompositeVector>& pres)
{
  double eps(1.e-8);
  int nfail = 0;

  int nfaces = pres->size("face");
  for (int f = 0; f != nfaces; ++f) {
    AmanziMesh::Entity_ID_List cells;
    cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::OWNED);

    if (cells.size() == 1) {
      if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) {
        if (std::abs((*pres)("face", f) - bc_values_[f]) > eps) nfail++;
      } else {
        if (bc_markers_[f] == Operators::OPERATOR_BC_NONE) { bc_values_[f] = 0.0; }

        AmanziGeometry::Point fpoint = mesh_->getFaceCentroid(f);
        AmanziGeometry::Point cpoint = mesh_->getCellCentroid(cells[0]);
        double dx = std::sqrt((fpoint - cpoint) * (fpoint - cpoint));

        double dp = std::abs(bc_values_[f]) * dx;
        double p = (*pres)("cell", cells[0]);
        if (bc_values_[f] > 0) {
          p = p + dp;
        } else {
          p = p - dp;
        }
        if (std::abs((*pres)("face", f) - p) > eps) nfail++;
      }
    } else {
      double p = ((*pres)("cell", cells[0]) + (*pres)("cell", cells[1])) / 2.0;
      if (std::abs((*pres)("face", f) - p) > eps) nfail++;
    }
  }

  if (nfail > 0) {
    return true;
  } else {
    return false;
  }
}

} // namespace TestPKs
} // namespace Amanzi
