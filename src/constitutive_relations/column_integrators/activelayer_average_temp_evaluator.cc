/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "activelayer_average_temp_evaluator.hh"

namespace Amanzi {
namespace Relations {

ParserActiveLayerAverageTemp::ParserActiveLayerAverageTemp(Teuchos::ParameterList& plist,
                                                           const KeyTag& key_tag)
{
  Key domain = Keys::getDomain(key_tag.first);
  Tag tag = key_tag.second;

  Key domain_ss = Keys::readDomainHint(plist, domain, "surface", "subsurface");
  Key temp_key = Keys::readKey(plist, domain_ss, "temperature", "temperature");
  dependencies.insert(KeyTag{ temp_key, tag });
}

IntegratorActiveLayerAverageTemp::IntegratorActiveLayerAverageTemp(
  Teuchos::ParameterList& plist,
  std::vector<const Epetra_MultiVector*>& deps,
  const AmanziMesh::Mesh* mesh)
  : mesh_(mesh)
{
  temp_ = deps[0];
  double trans_width = plist.get<double>("transition wdith [K]", 0.2);
  trans_temp_ = 273.15 + 0.5 * trans_width;
}

int
IntegratorActiveLayerAverageTemp::scan(AmanziMesh::Entity_ID col,
                                       AmanziMesh::Entity_ID c,
                                       AmanziGeometry::Point& p)
{
  if ((*temp_)[0][c] >= trans_temp_) {
    double cv = mesh_->getCellVolume(c);
    p[0] += (*temp_)[0][c] * cv;
    p[1] += cv;
    return false;
  }
  return true;
}

} // namespace Relations
} // namespace Amanzi
