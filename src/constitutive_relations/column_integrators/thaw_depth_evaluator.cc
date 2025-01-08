/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! An evaluator for calculating the water table height, relative to the surface.
#include "thaw_depth_evaluator.hh"

namespace Amanzi {
namespace Relations {

ParserThawDepth::ParserThawDepth(Teuchos::ParameterList& plist, const KeyTag& key_tag)
{
  Key domain = Keys::getDomain(key_tag.first);
  Tag tag = key_tag.second;

  Key domain_ss = Keys::readDomainHint(plist, domain, "surface", "subsurface");
  Key temp_key = Keys::readKey(plist, domain_ss, "temperature", "temperature");
  dependencies.insert(KeyTag{ temp_key, tag });

  Key cv_key = Keys::readKey(plist, domain_ss, "subsurface cell volume", "cell_volume");
  dependencies.insert(KeyTag{ cv_key, tag });

  Key cv_surf_key = Keys::readKey(plist, domain, "surface cell volume", "cell_volume");
  dependencies.insert(KeyTag{ cv_surf_key, tag });
}


IntegratorThawDepth::IntegratorThawDepth(Teuchos::ParameterList& plist,
                                         std::vector<const Epetra_MultiVector*>& deps,
                                         const AmanziMesh::Mesh* mesh)
{
  temp_ = deps[0];
  cv_ = deps[1];
  surf_cv_ = deps[2];

  double trans_width = plist.get<double>("transition width [K]", 0.2);
  trans_temp_ = 273.15 + 0.5 * trans_width;
}

int
IntegratorThawDepth::scan(AmanziMesh::Entity_ID col,
                          AmanziMesh::Entity_ID c,
                          AmanziGeometry::Point& p)
{
  if ((*temp_)[0][c] > trans_temp_) {
    p[0] += (*cv_)[0][c];
    return false;
  }
  return true;
}

double
IntegratorThawDepth::coefficient(AmanziMesh::Entity_ID col)
{
  return 1. / (*surf_cv_)[0][col];
}

} // namespace Relations
} // namespace Amanzi
