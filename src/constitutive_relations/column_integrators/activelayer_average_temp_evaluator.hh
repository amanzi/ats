/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Computes the average temperature in the active layer.

Evaluator name: `"active layer average temperature`"

.. _activelayer-average_evaluator-temp-spec
.. admonition:: activelayer-average_evaluator-temp-spec

   * `"transition width [K]`" ``[double]`` **0.2**

   KEYS:
   - `"temperature`"

*/

#ifndef AMANZI_FLOWRELATIONS_ALTTEMP_EVALUATOR_
#define AMANZI_FLOWRELATIONS_ALTTEMP_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorColumnIntegrator.hh"

namespace Amanzi {
namespace Relations {

struct ParserActiveLayerAverageTemp {
  ParserActiveLayerAverageTemp(Teuchos::ParameterList& plist, const KeyTag& key_tag);
  KeyTagSet dependencies;
};


class IntegratorActiveLayerAverageTemp {
 public:
  IntegratorActiveLayerAverageTemp(Teuchos::ParameterList& plist,
                                   std::vector<const Epetra_MultiVector*>& deps,
                                   const AmanziMesh::Mesh* mesh);
  int scan(AmanziMesh::Entity_ID col, AmanziMesh::Entity_ID c, AmanziGeometry::Point& p);
  double coefficient(AmanziMesh::Entity_ID col) { return 1; }

 private:
  double trans_temp_;
  const Epetra_MultiVector* temp_;
  const Epetra_MultiVector* cv_;
  const AmanziMesh::Mesh* mesh_;
};

using ActiveLayerAverageTempEvaluator =
  EvaluatorColumnIntegrator<ParserActiveLayerAverageTemp, IntegratorActiveLayerAverageTemp>;

} // namespace Relations
} // namespace Amanzi

#endif
