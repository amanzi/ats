/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

//! Sums a subsurface field vertically only a surface field.
#include "ColumnSumEvaluator.hh"

namespace Amanzi {
namespace Relations {
namespace Impl {

ParserColumnSum::ParserColumnSum(Teuchos::ParameterList& plist, const KeyTag& key_tag)
{
  Key surf_domain = Keys::getDomain(key_tag.first);
  Key domain = Keys::readDomainHint(plist, surf_domain, "surface", "subsurface");
  Key dep_key = Keys::readKey(plist, domain, "summed", Keys::getVarName(key_tag.first));
  dependencies.insert(KeyTag{ dep_key, key_tag.second });

  // dependency: cell volume, surface cell volume
  bool include_vol_factor = plist.get<bool>("include volume to surface area factor", false);
  if (include_vol_factor) {
    Key cv_key = Keys::readKey(plist, domain, "cell volume", "cell_volume");
    dependencies.insert(KeyTag{ cv_key, key_tag.second });

    Key surf_cv_key = Keys::readKey(plist, surf_domain, "surface cell volume", "cell_volume");
    dependencies.insert(KeyTag{ surf_cv_key, key_tag.second });
  }

  if (plist.get<bool>("volume averaged", false)) {
    if (include_vol_factor) {
      Errors::Message msg;
      msg << "ColumnSumEvaluator for " << key_tag.first
          << ": cannot use both options \"include volume to surface area factor\""
          << " and \"volume averaged\"";
      Exceptions::amanzi_throw(msg);
    }
    Key cv_key = Keys::readKey(plist, domain, "cell volume", "cell_volume");
    dependencies.insert(KeyTag{ cv_key, key_tag.second });
  }

  if (plist.get<bool>("divide by density", false)) {
    Key molar_dens_key = Keys::readKey(plist, domain, "molar density", "molar_density_liquid");
    dependencies.insert(KeyTag{ molar_dens_key, key_tag.second });
  }
}


IntegratorColumnSum::IntegratorColumnSum(Teuchos::ParameterList& plist,
                                         std::vector<const Epetra_MultiVector*>& deps,
                                         const AmanziMesh::Mesh* mesh)
{
  int i_dep(0);
  integrand_ = deps[i_dep++];

  volume_factor_ = plist.get<bool>("include volume to surface area factor");
  volume_average_ = plist.get<bool>("volume averaged");
  divide_by_density_ = plist.get<bool>("divide by density");

  if (volume_factor_) {
    AMANZI_ASSERT(deps.size() >= 3);
    cv_ = deps[i_dep++];
    surf_cv_ = deps[i_dep++];
  }
  if (volume_average_) {
    AMANZI_ASSERT(deps.size() >= 2);
    cv_ = deps[i_dep++];
  }

  if (divide_by_density_) {
    AMANZI_ASSERT(deps.size() > i_dep);
    dens_ = deps[i_dep];
  }

  coef_ = plist.get<double>("coefficient", 1.0);
}


int
IntegratorColumnSum::scan(AmanziMesh::Entity_ID col,
                          AmanziMesh::Entity_ID c,
                          AmanziGeometry::Point& p)
{
  double contrib = (*integrand_)[0][c];
  if (volume_average_ || volume_factor_) { contrib *= (*cv_)[0][c]; }
  if (divide_by_density_) { contrib /= (*dens_)[0][c]; }
  p[0] += contrib;

  if (volume_average_) p[1] += (*cv_)[0][c];
  return false;
}


double
IntegratorColumnSum::coefficient(AmanziMesh::Entity_ID col)
{
  if (volume_factor_) {
    return coef_ / (*surf_cv_)[0][col];
  } else {
    return coef_;
  }
}


} //namespace Impl
} //namespace Relations
} //namespace Amanzi
