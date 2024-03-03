/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Calculates the water table depth where the water table is defined as the top cell with 0 gas saturation.
#include "water_table_depth_evaluator.hh"

namespace Amanzi {
namespace Relations {

ParserWaterTableDepth::ParserWaterTableDepth(Teuchos::ParameterList& plist, const KeyTag& key_tag)
{
  Key domain = Keys::getDomain(key_tag.first);
  Tag tag = key_tag.second;

  Key domain_ss = Keys::readDomainHint(plist, domain, "surface", "subsurface");
  Key sat_key = Keys::readKey(plist, domain_ss, "saturation of gas", "saturation_gas");
  dependencies.insert(KeyTag{ sat_key, tag });

  Key cv_key = Keys::readKey(plist, domain_ss, "subsurface cell volume", "cell_volume");
  dependencies.insert(KeyTag{ cv_key, tag });

  Key cv_surf_key = Keys::readKey(plist, domain, "surface cell volume", "cell_volume");
  dependencies.insert(KeyTag{ cv_surf_key, tag });
}


IntegratorWaterTableDepth::IntegratorWaterTableDepth(Teuchos::ParameterList& plist,
                                                     std::vector<cView_type>& deps,
                                                     const AmanziMesh::Mesh& mesh)
{
  assert(deps.size() == 3);
  sat_ = deps[0];
  cv_ = deps[1];
  surf_cv_ = deps[2];
}

} // namespace Relations
} // namespace Amanzi
