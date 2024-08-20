/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

//! Calculates the water table depth by locating the end of continuously unsaturated cells from top downwards.
#include "perched_water_table_depth_evaluator.hh"

namespace Amanzi {
namespace Relations {

ParserPerchedWaterTableDepth::ParserPerchedWaterTableDepth(Teuchos::ParameterList& plist, const KeyTag& key_tag)
{
  Key domain = Keys::getDomain(key_tag.first);
  Tag tag = key_tag.second;

  Key domain_ss = Keys::readDomainHint(plist, domain, "surface", "subsurface");
  Key sat_key = Keys::readKey(plist, domain_ss, "saturation of gas", "saturation_gas");
  dependencies.insert(KeyTag{ sat_key, tag });

  Key pres_key = Keys::readKey(plist, domain_ss, "pressure", "pressure");
  dependencies.insert(KeyTag{ pres_key, tag });

  Key cv_key = Keys::readKey(plist, domain_ss, "subsurface cell volume", "cell_volume");
  dependencies.insert(KeyTag{ cv_key, tag });

  Key cv_surf_key = Keys::readKey(plist, domain, "surface cell volume", "cell_volume");
  dependencies.insert(KeyTag{ cv_surf_key, tag });
}


IntegratorPerchedWaterTableDepth::IntegratorPerchedWaterTableDepth(Teuchos::ParameterList& plist,
                                                     std::vector<const Epetra_MultiVector*>& deps,
                                                     const AmanziMesh::Mesh* mesh) : mesh_(mesh)
{
  AMANZI_ASSERT(deps.size() == 4);
  sat_ = deps[0];
  pres_ = deps[1];
  cv_ = deps[2];
  surf_cv_ = deps[3];
  is_interp_ = plist.get<bool>("interpolate depth from pressure", false);
}

int
IntegratorPerchedWaterTableDepth::scan(AmanziMesh::Entity_ID col,
                                       AmanziMesh::Entity_ID c,
                                       AmanziGeometry::Point& p)
{
  if ((*sat_)[0][c] > 0.0) {
    p[2] = mesh_->getCellCentroid(c)[2]; // last unsaturated cell centroid
    if (is_interp_) {
      p[0] = (*pres_)[0][c]; // last unsaturated cell pressure
    } else {
      p[0] += (*cv_)[0][c]; // cumulative unsaturated cell volume
    }
    return false;
  }
  if (is_interp_) {
    p[1] = (*pres_)[0][c]; // first saturated cell pressure
  }
  return true;
}

double
IntegratorPerchedWaterTableDepth::coefficient(AmanziMesh::Entity_ID col)
{
  return 1. / (*surf_cv_)[0][col];
}

} // namespace Relations
} // namespace Amanzi
