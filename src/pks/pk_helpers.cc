/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A set of helper functions for doing common things in PKs.
#include "Mesh_Algorithms.hh"
#include "Chemistry_PK.hh"
#include "pk_helpers.hh"

namespace Amanzi {
namespace PKHelpers {

bool
aliasVector(State& S, const Key& key, const Tag& target, const Tag& alias)
{
  if (S.HasEvaluator(key, target) && !S.HasEvaluator(key, alias)) {
    S.SetEvaluator(key, alias, S.GetEvaluatorPtr(key, target));
    S.GetRecordSetW(key).AliasRecord(target, alias);
    return true;
  }
  return false;
}


// -----------------------------------------------------------------------------
// Given a vector, apply the Dirichlet data to that vector.
// -----------------------------------------------------------------------------
void
applyDirichletBCs(const Operators::BCs& bcs, CompositeVector& u)
{
  if (u.hasComponent("face")) {
    auto u_f = u.viewComponent("face", false);
    auto bc_value = bcs.bc_value();
    auto bc_model = bcs.bc_model();
    Kokkos::parallel_for("applyDirichletBCs", u_f.extent(0),
                         KOKKOS_LAMBDA(const int& f) {
                           if (bc_model(f) == Operators::OPERATOR_BC_DIRICHLET) {
                             u_f(f,0) = bc_value(f);
                           }
                         });
  }

  if (u.hasComponent("boundary_face")) {
    auto u_f = u.viewComponent("boundary_face", false);
    auto bc_value = bcs.bc_value();
    auto bc_model = bcs.bc_model();

    const AmanziMesh::Mesh* mesh = u.getMesh().get();

    Kokkos::parallel_for("applyDirichletBCs", u_f.extent(0),
                         KOKKOS_LAMBDA(const int& bf) {
                           auto f = getBoundaryFaceFace(*mesh, bf);
                           if (bc_model(f) == Operators::OPERATOR_BC_DIRICHLET) {
                             u_f(bf,0) = bc_value(f);
                           }
                         });
  }
}


// // -----------------------------------------------------------------------------
// // Given a vector and a face ID, get the value at that location.
// //
// // Looks in the following order:
// //  -- face component
// //  -- boundary Dirichlet data
// //  -- boundary_face value
// //  -- internal cell
// // -----------------------------------------------------------------------------
// double
// getFaceOnBoundaryValue(AmanziMesh::Entity_ID f, const CompositeVector& u, const Operators::BCs& bcs)
// {
//   if (u.hasComponent("face")) {
//     return (*u.viewComponent("face", false))[0][f];
//   } else if (bcs.bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
//     return bcs.bc_value()[f];
//     // } else if (u.hasComponent("boundary_face")) {
//     //   AmanziMesh::Entity_ID bf = getFaceOnBoundaryBoundaryFace(*u.getMesh(), f);
//     //   return (*u.viewComponent("boundary_face",false))[0][bf];
//   } else {
//     auto c = getFaceOnBoundaryInternalCell(*u.getMesh(), f);
//     return (*u.viewComponent("cell", false))[0][c];
//   }
//   return -1;
// }


// // -----------------------------------------------------------------------------
// // Get the directional int for a face that is on the boundary.
// // -----------------------------------------------------------------------------
// int
// getBoundaryDirection(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f)
// {
//   AmanziMesh::Entity_ID_List cells;
//   cells = mesh.getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
//   AMANZI_ASSERT(cells.size() == 1);
//   AmanziMesh::Entity_ID_List faces;
//   std::vector<int> dirs;
//   mesh.getCellFacesAndDirections(cells[0], &faces, &dirs);
//   return dirs[std::find(faces.begin(), faces.end(), f) - faces.begin()];
// }


// -----------------------------------------------------------------------------
// Get a primary variable evaluator for a key at tag
// -----------------------------------------------------------------------------
Teuchos::RCP<EvaluatorPrimaryCV>
requireEvaluatorPrimary(const Key& key, const Tag& tag, State& S, bool or_die)
{
  // first check, is there one already
  if (S.HasEvaluator(key, tag)) {
    // if so, make sure it is primary
    Teuchos::RCP<Evaluator> eval = S.GetEvaluatorPtr(key, tag);
    Teuchos::RCP<EvaluatorPrimaryCV> eval_pv = Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(eval);
    if (or_die && eval_pv == Teuchos::null) {
      Errors::Message msg;
      msg << "Expected primary variable evaluator for " << key << " @ " << tag.get();
      Exceptions::amanzi_throw(msg);
    }
    return eval_pv;
  }

  // if not, create one, only at this tag, not to be shared across tags.  By
  // this, we mean we don't stick the "type" = "primary" back into the
  // evaluator list -- this allows "copy evaluators" e.g. "water content at the
  // old tag" to differ from the standard evalulator, e.g. "water content at
  // the new tag" which is likely a secondary variable evaluator.
  Teuchos::ParameterList plist(key);
  plist.set("evaluator type", "primary variable");
  plist.set("tag", tag.get());
  auto eval_pv = Teuchos::rcp(new EvaluatorPrimaryCV(plist));
  S.SetEvaluator(key, tag, eval_pv);
  return eval_pv;
}


// -----------------------------------------------------------------------------
// Marks a primary evaluator as changed.
// -----------------------------------------------------------------------------
bool
changedEvaluatorPrimary(const Key& key, const Tag& tag, State& S, bool or_die)
{
  bool changed = false;
  Teuchos::RCP<Evaluator> eval = S.GetEvaluatorPtr(key, tag);
  Teuchos::RCP<EvaluatorPrimaryCV> eval_pv = Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(eval);
  if (eval_pv == Teuchos::null) {
    if (or_die) {
      Errors::Message msg;
      msg << "Expected primary variable evaluator for " << key << " @ " << tag.get();
      Exceptions::amanzi_throw(msg);
    }
  } else {
    eval_pv->SetChanged();
    changed = true;
  }
  return changed;
}


// -----------------------------------------------------------------------------
// Require a vector and a primary variable evaluator at current tag(s).
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireAtCurrent(const Key& key, const Tag& tag, State& S, const Key& name, bool is_eval)
{
  CompositeVectorSpace& cvs = S.Require<CompositeVector, CompositeVectorSpace>(key, tag);
  if (!name.empty()) {
    Key owner = S.GetRecord(key, tag).owner();
    if (owner.empty()) {
      S.Require<CompositeVector, CompositeVectorSpace>(key, tag, name);
      if (is_eval) requireEvaluatorAssign(key, tag, S);
    }

    if (tag != Tags::CURRENT) {
      S.Require<CompositeVector, CompositeVectorSpace>(key, Tags::CURRENT);
      Key current_owner = S.GetRecord(key, Tags::CURRENT).owner();
      if (owner.empty()) {
        S.Require<CompositeVector, CompositeVectorSpace>(key, Tags::CURRENT, name);
        if (is_eval) requireEvaluatorAssign(key, Tags::CURRENT, S);
      }
    }
  } else {
    if (is_eval) S.RequireEvaluator(key, tag);
  }
  return cvs;
}


// -----------------------------------------------------------------------------
// Require a vector and a primary variable evaluator at next tag(s).
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireAtNext(const Key& key, const Tag& tag, State& S, const Key& name)
{
  CompositeVectorSpace& cvs = S.Require<CompositeVector, CompositeVectorSpace>(key, tag);
  if (!name.empty()) {
    S.Require<CompositeVector, CompositeVectorSpace>(key, tag, name);
    requireEvaluatorPrimary(key, tag, S);
  } else {
    S.RequireEvaluator(key, tag);
  }

  if (tag != Tags::NEXT) { aliasVector(S, key, tag, Tags::NEXT); }
  return cvs;
}


// -----------------------------------------------------------------------------
// Require assignment evaluator, which allows tracking old data.
// -----------------------------------------------------------------------------
Teuchos::RCP<EvaluatorPrimaryCV>
requireEvaluatorAssign(const Key& key, const Tag& tag, State& S)
{
  // in the future, this will likely derive from primary instead of just being
  // primary.  This will allow confirming that the times are the same.
  return requireEvaluatorPrimary(key, tag, S, false);
}

// -----------------------------------------------------------------------------
// Assign if it is an assignment evaluator.
// -----------------------------------------------------------------------------
void
assign(const Key& key, const Tag& tag_dest, const Tag& tag_source, State& S)
{
  S.GetEvaluator(key, tag_source).Update(S, Keys::getKey(key, tag_dest));
  bool changed = changedEvaluatorPrimary(key, tag_dest, S, false);
  if (changed) S.Assign(key, tag_dest, tag_source);
}


// // -----------------------------------------------------------------------------
// // Helper functions for working with Amanzi's Chemistry PK
// // -----------------------------------------------------------------------------
// void
// convertConcentrationToAmanzi(const Epetra_MultiVector& mol_dens,
//                              int num_aqueous,
//                              const Epetra_MultiVector& tcc_ats,
//                              Epetra_MultiVector& tcc_amanzi)
// {
//   // convert from mole fraction [mol C / mol H20] to [mol C / L]
//   for (int k = 0; k != num_aqueous; ++k) {
//     for (int c = 0; c != tcc_ats.MyLength(); ++c) {
//       // 1.e-3 converts L to m^3
//       tcc_amanzi[k][c] = tcc_ats[k][c] * mol_dens[0][c] * 1.e-3;
//     }
//   }
// }


// void
// convertConcentrationToATS(const Epetra_MultiVector& mol_dens,
//                           int num_aqueous,
//                           const Epetra_MultiVector& tcc_amanzi,
//                           Epetra_MultiVector& tcc_ats)
// {
//   // convert from [mol C / L] to mol fraction [mol C / mol H20]
//   for (int k = 0; k != num_aqueous; ++k) {
//     for (int c = 0; c != tcc_amanzi.MyLength(); ++c) {
//       tcc_ats[k][c] = tcc_amanzi[k][c] / (mol_dens[0][c] * 1.e-3);
//     }
//   }
// }


// bool
// advanceChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
//                  double t_old,
//                  double t_new,
//                  bool reinit,
//                  const Epetra_MultiVector& mol_dens,
//                  Teuchos::RCP<Epetra_MultiVector> tcc,
//                  Teuchos::Time& timer)
// {
//   bool fail = false;
//   int num_aqueous = chem_pk->num_aqueous_components();
//   convertConcentrationToAmanzi(mol_dens, num_aqueous, *tcc, *tcc);
//   chem_pk->set_aqueous_components(tcc);

//   {
//     auto monitor = Teuchos::rcp(new Teuchos::TimeMonitor(timer));
//     fail = chem_pk->AdvanceStep(t_old, t_new, reinit);
//   }
//   if (fail) return fail;

//   *tcc = *chem_pk->aqueous_components();
//   convertConcentrationToATS(mol_dens, num_aqueous, *tcc, *tcc);
//   return fail;
// }

// -------------------------------------------------------------
// Helper function and customization point for upwinded coefs.
// -------------------------------------------------------------
void
requireNonlinearDiffusionCoefficient(const Key& key, const Tag& tag,
                            const std::string& coef_location, State& S)
{
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S.GetMesh(Keys::getDomain(key));
  if (coef_location == "upwind: face") {
    S.Require<CompositeVector, CompositeVectorSpace>(key, tag, key)
      .SetMesh(mesh)
      ->SetGhosted()
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  } else if (coef_location == "standard: cell") {
    S.Require<CompositeVector, CompositeVectorSpace>(key, tag, key)
      .SetMesh(mesh)
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  } else {
    Errors::Message msg;
    msg << "Unknown upwind coefficient location for \"" << key << "\"";
    Exceptions::amanzi_throw(msg);
  }
  S.GetRecordW(key, tag, key).set_io_vis(false);
}





} // namespace PKHelpers
} // namespace Amanzi
