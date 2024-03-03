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
#include "BlockVector_decl.hh"
#include "EvaluatorColumnIntegrator.hh"

namespace Amanzi {
namespace Relations {

struct ParserWaterTableDepth {
  ParserWaterTableDepth(Teuchos::ParameterList& plist, const KeyTag& key_tag);
  KeyTagSet dependencies;
};


class IntegratorWaterTableDepth {
 public:
  using cView_type = BlockVector<double>::cMultiVectorView_type<Amanzi::DefaultDevice>;

  IntegratorWaterTableDepth(Teuchos::ParameterList& plist,
                            std::vector<cView_type>& deps,
                            const AmanziMesh::Mesh& mesh);

  KOKKOS_INLINE_FUNCTION
  int scan(const AmanziMesh::Entity_ID col, const AmanziMesh::Entity_ID c, AmanziGeometry::Point& p) const;

  KOKKOS_INLINE_FUNCTION
  double coefficient(const AmanziMesh::Entity_ID col) const;

 private:
  cView_type sat_;
  cView_type cv_;
  cView_type surf_cv_;
};

using WaterTableDepthEvaluator =
  EvaluatorColumnIntegrator<ParserWaterTableDepth, IntegratorWaterTableDepth>;


KOKKOS_INLINE_FUNCTION
int
IntegratorWaterTableDepth::scan(const AmanziMesh::Entity_ID col,
                                const AmanziMesh::Entity_ID c,
                                AmanziGeometry::Point& p) const
{
  if (sat_(c, 0) > 0.0) {
    p[0] += cv_(c, 0);
    return false;
  }
  return true;
}

KOKKOS_INLINE_FUNCTION
double
IntegratorWaterTableDepth::coefficient(const AmanziMesh::Entity_ID col) const
{
  return 1. / surf_cv_(col, 0);
}



} //namespace Relations
} //namespace Amanzi
