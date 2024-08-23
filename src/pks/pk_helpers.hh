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
#include "Chemistry_PK.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Given a vector, apply the Dirichlet data to that vector's boundary_face
// component.
// -----------------------------------------------------------------------------
void
applyDirichletBCs(const Operators::BCs& bcs, CompositeVector& u);


// -----------------------------------------------------------------------------
// Given a vector and a face ID, get the value at that location.
//
// Looks in the following order:
//  -- face component
//  -- boundary Dirichlet data
//  -- boundary_face value (currently not used -- fix me --etc)
//  -- internal cell
// -----------------------------------------------------------------------------
double
getFaceOnBoundaryValue(AmanziMesh::Entity_ID f,
                       const CompositeVector& u,
                       const Operators::BCs& bcs);


// -----------------------------------------------------------------------------
// Get the directional int for a face that is on the boundary.
// -----------------------------------------------------------------------------
int
getBoundaryDirection(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f);


// -----------------------------------------------------------------------------
// Create an alias -- a pointer is copied such that Tag alias's evaluator and
// data both point to Tag target's.
// -----------------------------------------------------------------------------
bool
aliasVector(State& S, const Key& key, const Tag& target, const Tag& alias);

// -----------------------------------------------------------------------------
// Get a primary variable evaluator for a key at tag
//
// If or_die is true, throw an error if the xml has a sublist for this key.
// -----------------------------------------------------------------------------
Teuchos::RCP<EvaluatorPrimaryCV>
requireEvaluatorPrimary(const Key& key, const Tag& tag, State& S, bool or_die = true);


// -----------------------------------------------------------------------------
// Mark primary variable evaluator as changed.
//
// If or_die is true, throw an error if GetEvaluator(key,tag) is not castable
// to EvaluatorPrimary
// -----------------------------------------------------------------------------
bool
changedEvaluatorPrimary(const Key& key, const Tag& tag, State& S, bool or_die = true);


// -----------------------------------------------------------------------------
// Require a vector and an evaluator at current tag(s).
//
// If owner is not empty, this variable is claimed and a EvaluatorPrimary is
// created.  Also, a CURRENT variable is also required, independent of tag, to
// ensure that we have a way to recover from failed timesteps.  This is a
// separate copy of data, and an "assigment" evaluator.
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireAtCurrent(const Key& key,
                 const Tag& tag,
                 State& S,
                 const Key& owner = "");

// -----------------------------------------------------------------------------
// Require a vector and a primary variable evaluator at next tag(s).
//
// If managed_here is true, this indicates that this variable is a quantity
// that will definitely be computed at tag, and the calling PK will take
// responsibility for making sure that key@tag is also copied to key@NEXT,
// thereby allowing other PKs to use it.  Effectively this creates an aliased
// evaluator from key@NEXT --> key@tag.
//
// If owner is not empty, key@tag is claimed by owner, and a EvaluatorPrimary
// is created.  This also implies managed_here.
// -----------------------------------------------------------------------------
CompositeVectorSpace&
requireAtNext(const Key& key,
              const Tag& tag,
              State& S,
              bool managed_here,
              const Key& owner = "");

inline CompositeVectorSpace&
requireAtNext(const Key& key,
              const Tag& tag,
              State& S,
              const Key& owner = "") {
  return requireAtNext(key, tag, S, false, owner);
}




// -----------------------------------------------------------------------------
// Require assignment evaluator, which allows tracking old data.
// -----------------------------------------------------------------------------
Teuchos::RCP<EvaluatorPrimaryCV>
requireEvaluatorAssign(const Key& key, const Tag& tag, State& S);

// -----------------------------------------------------------------------------
// Assign if it is an assignment evaluator.
// -----------------------------------------------------------------------------
void
assign(const Key& key, const Tag& tag_dest, const Tag& tag_source, State& S);


// -----------------------------------------------------------------------------
// Helper functions for working with Amanzi's Chemistry PK
// -----------------------------------------------------------------------------
void
convertConcentrationToAmanzi(const Epetra_MultiVector& mol_den,
                             int num_aqueous,
                             const Epetra_MultiVector& tcc_ats,
                             Epetra_MultiVector& tcc_amanzi);

void
convertConcentrationToATS(const Epetra_MultiVector& mol_den,
                          int num_aqueous,
                          const Epetra_MultiVector& tcc_ats,
                          Epetra_MultiVector& tcc_amanzi);

bool
advanceChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                 double t_old,
                 double t_new,
                 bool reinit,
                 const Epetra_MultiVector& mol_dens,
                 Teuchos::RCP<Epetra_MultiVector> tcc,
                 Teuchos::Time& timer);


void
copyMeshCoordinatesToVector(const AmanziMesh::Mesh& mesh, CompositeVector& vec);
void
copyVectorToMeshCoordinates(const CompositeVector& vec, AmanziMesh::Mesh& mesh);


// Compute pairs of value + location
typedef struct ValLoc {
  double value;
  AmanziMesh::Entity_ID gid;
} ENorm_t;

int
commMaxValLoc(const Comm_type& comm, const ValLoc& local, ValLoc& global);
ValLoc
maxValLoc(const Epetra_Vector& vec);

} // namespace Amanzi
