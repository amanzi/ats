/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the Kr associated with the unfrozen fraction of water.

*/

#include <cmath>

#include "dbc.hh"
#include "errors.hh"
#include "one_uf_relperm_model.hh"

namespace Amanzi {
namespace Flow {

OneUFRelPermModel::OneUFRelPermModel(Teuchos::ParameterList& plist) : plist_(plist), pi_(M_PI)
{
  alpha_ = plist_.get<int>("unfrozen rel perm alpha", 4);
  if (alpha_ % 2 != 0) {
    Errors::Message message("Unfrozen Fraction Rel Perm: alpha must be an even integer");
    Exceptions::amanzi_throw(message);
  }

  h_cutoff_up_ = plist_.get<double>("unfrozen rel perm cutoff pressure [Pa]", 10.);
  h_cutoff_dn_ = plist_.get<double>("unfrozen rel perm cutoff pressure, below [Pa]", 0.);

  if (plist_.isParameter("unfrozen rel perm cutoff pressure") ||
      plist_.isParameter("unfrozen rel perm cutoff height")) {
    Errors::Message message(
      "surface-relative_permeability: old-style parameters detected, please use \"surface rel perm "
      "cutoff pressure [Pa]\" and \"surface rel perm cutoff pressure, below [Pa]\", not \"surface "
      "rel perm cutoff height\" or \"surface rel perm cutoff pressure\"");
    Exceptions::amanzi_throw(message);
  }
}

double
OneUFRelPermModel::SurfaceRelPerm(double uf, double h)
{
  double kr;

  if (h >= 101325. + h_cutoff_up_) {
    kr = std::pow(std::sin(pi_ * uf / 2.), alpha_);
  } else if (h <= 101325 + h_cutoff_dn_) {
    kr = 1.;
  } else {
    double fac = (h - (101325 + h_cutoff_dn_)) / (h_cutoff_up_ - h_cutoff_dn_);
    kr = (1 - fac) + fac * std::pow(std::sin(pi_ * uf / 2.), alpha_);
  }
  return kr;
}


} // namespace Flow
} // namespace Amanzi
