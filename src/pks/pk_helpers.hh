/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A set of helper functions for doing common things in PKs.
#pragma once

#include "Teuchos_TimeMonitor.hpp"

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "BCs.hh"

#include "EvaluatorPrimary.hh"
#include "State.hh"

namespace Amanzi {
namespace PKHelpers {

bool
aliasVector(State& S, const Key& key, const Tag& target, const Tag& alias);

// -----------------------------------------------------------------------------
// Given a vector, apply the Dirichlet data to that vector's boundary_face
// component.
// -----------------------------------------------------------------------------
void
applyDirichletBCs(const Operators::BCs& bcs, CompositeVector& u);


// // -----------------------------------------------------------------------------
// // Given a vector and a face ID, get the value at that location.
// //
// // Looks in the following order:
// //  -- face component
// //  -- boundary Dirichlet data
// //  -- boundary_face value (currently not used -- fix me --etc)
// //  -- internal cell
// // -----------------------------------------------------------------------------
// double
// getFaceOnBoundaryValue(AmanziMesh::Entity_ID f,
//                        const CompositeVector& u,
//                        const Operators::BCs& bcs);


// // -----------------------------------------------------------------------------
// // Get the directional int for a face that is on the boundary.
// // -----------------------------------------------------------------------------
// int
// getBoundaryDirection(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f);


// -----------------------------------------------------------------------------
// Get a primary variable evaluator for a key at tag
// -----------------------------------------------------------------------------
Teuchos::RCP<EvaluatorPrimaryCV>
requireEvaluatorPrimary(const Key& key, const Tag& tag, State& S, bool or_die = true);


// -----------------------------------------------------------------------------
// Mark primary variable evaluator as changed.
// -----------------------------------------------------------------------------
bool
changedEvaluatorPrimary(const Key& key, const Tag& tag, State& S, bool or_die = true);


// -----------------------------------------------------------------------------
// Require a vector and a primary variable evaluator at current tag(s).
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireAtCurrent(const Key& key,
                 const Tag& tag,
                 State& S,
                 const Key& name = "",
                 bool is_eval = true);


// -----------------------------------------------------------------------------
// Require a vector and a primary variable evaluator at next tag(s).
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireAtNext(const Key& key, const Tag& tag, State& S, const Key& name = "");

// -----------------------------------------------------------------------------
// Require assignment evaluator, which allows tracking old data.
// -----------------------------------------------------------------------------
Teuchos::RCP<EvaluatorPrimaryCV>
requireEvaluatorAssign(const Key& key, const Tag& tag, State& S);

// -----------------------------------------------------------------------------
// Assign if it is an assignment evaluator, used in CommitStep to advance the
// interval
// -----------------------------------------------------------------------------
void
assign(const Key& key, const Tag& tag_dest, const Tag& tag_source, State& S);


// -----------------------------------------------------------------------------
// Helper functions for working with Amanzi's Chemistry PK
// -----------------------------------------------------------------------------
// void
// convertConcentrationToAmanzi(const Epetra_MultiVector& mol_den,
//                              int num_aqueous,
//                              const Epetra_MultiVector& tcc_ats,
//                              Epetra_MultiVector& tcc_amanzi);

// void
// convertConcentrationToATS(const Epetra_MultiVector& mol_den,
//                           int num_aqueous,
//                           const Epetra_MultiVector& tcc_ats,
//                           Epetra_MultiVector& tcc_amanzi);

// bool
// advanceChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
//                  double t_old,
//                  double t_new,
//                  bool reinit,
//                  const Epetra_MultiVector& mol_dens,
//                  Teuchos::RCP<Epetra_MultiVector> tcc,
//                  Teuchos::Time& timer);


// -----------------------------------------------------------------------------
// Require an upwinded coeficient for diffusion operators, given a location
// defined by the upwinding scheme.
// -----------------------------------------------------------------------------
void
requireNonlinearDiffusionCoefficient(const Key& key,
                                     const Tag& tag,
                                     const std::string& coef_location,
                                     State& S);

} // namespace PKHelpers
} // namespace Amanzi
