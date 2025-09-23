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
#include "unfrozen_fraction_relperm_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

UnfrozenFractionRelPermModel::UnfrozenFractionRelPermModel(Teuchos::ParameterList& plist)
  : plist_(plist), pi_(M_PI)
{
  alpha_ = plist_.get<int>("unfrozen rel perm alpha", 4);
  if (alpha_ % 2 != 0) {
    Errors::Message message("Unfrozen Fraction Rel Perm: alpha must be an even integer");
    Exceptions::amanzi_throw(message);
  }
}

double
UnfrozenFractionRelPermModel::SurfaceRelPerm(double uf, double h)
{
  return std::pow(std::sin(pi_ * uf / 2.), alpha_);
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
