/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@ornl.gov)
*/

//! A set of helper functions for doing common things in PKs.

#pragma once

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "BCs.hh"

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
getFaceOnBoundaryValue(AmanziMesh::Entity_ID f, const CompositeVector& u, const Operators::BCs& bcs);


// -----------------------------------------------------------------------------
// Get the directional int for a face that is on the boundary.
// -----------------------------------------------------------------------------
int
getBoundaryDirection(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f);


} // namespace Amanzi
