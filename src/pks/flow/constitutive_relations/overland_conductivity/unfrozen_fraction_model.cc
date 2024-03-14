/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the unfrozen fraction of water.

*/

#include <cmath>

#include "errors.hh"
#include "unfrozen_fraction_model.hh"

namespace Amanzi {
namespace Flow {

UnfrozenFractionModel::UnfrozenFractionModel(Teuchos::ParameterList& plist)
  : plist_(plist), pi_(M_PI)
{
  if (plist_.isParameter("transition width")) {
    Errors::Message message("Unfrozen Fraction Evaluator: parameter changed from \"transition "
                            "width\" to \"transition width [K]\"");
    Exceptions::amanzi_throw(message);
  }
  halfwidth_ = plist_.get<double>("transition width [K]", 0.2) / 2.;

  if (plist_.isParameter("freezing point")) {
    Errors::Message message("Unfrozen Fraction Evaluator: parameter changed from \"freezing "
                            "point\" to \"freezing point [K]\"");
    Exceptions::amanzi_throw(message);
  }
  T0_ = plist_.get<double>("freezing point [K]", 273.15);

  min_uf_ = plist_.get<double>("minimum unfrozen fraction [-]", 0.);
}

double
UnfrozenFractionModel::UnfrozenFraction(double temp) const
{
  double adj_temp = temp - T0_;
  double uf;
  if (adj_temp > halfwidth_) {
    uf = 1.;
  } else if (adj_temp < -halfwidth_) {
    uf = min_uf_;
  } else {
    uf = std::max(min_uf_, (std::sin(pi_ / 2. * adj_temp / halfwidth_) + 1.) / 2.);
  }
  return uf;
}

double
UnfrozenFractionModel::DUnfrozenFractionDT(double temp) const
{
  double adj_temp = temp - T0_;
  double duf;
  if (adj_temp > halfwidth_) {
    duf = 0.;
  } else if (adj_temp < -halfwidth_) {
    duf = 0.;
  } else {
    duf = std::cos(pi_ / 2. * adj_temp / halfwidth_) / 2. * pi_ / 2. / halfwidth_;
  }
  return duf;
}

} // namespace Flow
} // namespace Amanzi
