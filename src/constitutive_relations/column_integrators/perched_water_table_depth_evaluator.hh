/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

Computes the depth to the perched water table.

`"evaluator type`" = `"perched water table depth`"

.. _perched_water_table_depth-spec:
.. admonition:: perched_water_table_depth-spec

    * `"interpolate depth from pressure`" ``[bool]`` **false** Default to calculate
      perched water table depth by locating the bottom face of the last continuously
      unsaturated cell from top downward. If true, use the height and pressure at the
      centroids of the last continuously unsaturated cell and its adjacent saturated
      cell to determine the perched water table depth through interpolation.

    KEYS:

    - `"saturation of gas`" **SUBSURFACE_DOMAIN-saturation-gas**
    - `"pressure`" **SUBSURFACE_DOMAIN-pressure**
    - `"subsurface cell volume`" **SUBSURFACE_DOMAIN-cell_volume**
    - `"surface cell volume`" **DOMAIN-cell_volume**


*/

#pragma once

#include "Factory.hh"
#include "PerchedWaterTableColumnIntegrator.hh"

namespace Amanzi {
namespace Relations {

struct ParserPerchedWaterTableDepth {
  ParserPerchedWaterTableDepth(Teuchos::ParameterList& plist, const KeyTag& key_tag);
  KeyTagSet dependencies;
};


class IntegratorPerchedWaterTableDepth {
 public:
  IntegratorPerchedWaterTableDepth(Teuchos::ParameterList& plist,
                                   std::vector<const Epetra_MultiVector*>& deps,
                                   const AmanziMesh::Mesh* mesh);
  int scan(AmanziMesh::Entity_ID col, AmanziMesh::Entity_ID c, AmanziGeometry::Point& p);
  double coefficient(AmanziMesh::Entity_ID col);

 private:
  const Epetra_MultiVector* sat_;
  const Epetra_MultiVector* pres_;
  const Epetra_MultiVector* cv_;
  const Epetra_MultiVector* surf_cv_;
  const AmanziMesh::Mesh* mesh_;
  bool is_interp_;
};

using PerchedWaterTableDepthEvaluator =
  PerchedWaterTableColumnIntegrator<ParserPerchedWaterTableDepth, IntegratorPerchedWaterTableDepth>;

} //namespace Relations
} //namespace Amanzi
