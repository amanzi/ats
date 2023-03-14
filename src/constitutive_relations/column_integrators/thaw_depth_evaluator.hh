/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Computes the depth to a saturated water table.

Evaluator name: `"water table depth`"

.. _water-table-depth-spec:
.. admonition:: water-table-depth-spec

    KEYS:
      `"saturation_gas`"

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorColumnIntegrator.hh"

namespace Amanzi {
namespace Relations {

struct ParserThawDepth {
  ParserThawDepth(Teuchos::ParameterList& plist, const KeyTag& key_tag);
  KeyTagSet dependencies;
};


class IntegratorThawDepth {
 public:
  IntegratorThawDepth(Teuchos::ParameterList& plist,
                      std::vector<const Epetra_MultiVector*>& deps,
                      const AmanziMesh::Mesh* mesh);
  int scan(AmanziMesh::Entity_ID col, AmanziMesh::Entity_ID c, AmanziGeometry::Point& p);
  double coefficient(AmanziMesh::Entity_ID col);

 private:
  const Epetra_MultiVector* temp_;
  const Epetra_MultiVector* cv_;
  const Epetra_MultiVector* surf_cv_;
  double trans_temp_;
};

using ThawDepthEvaluator = EvaluatorColumnIntegrator<ParserThawDepth, IntegratorThawDepth>;

} //namespace Relations
} //namespace Amanzi
