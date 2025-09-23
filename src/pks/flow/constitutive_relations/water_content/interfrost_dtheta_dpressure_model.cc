/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The interfrost dtheta_dpressure model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Interfrost water content portion sl.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "interfrost_dtheta_dpressure_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
InterfrostDthetaDpressureModel::InterfrostDthetaDpressureModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
InterfrostDthetaDpressureModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  beta_ = plist.get<double>("compressibility [1/Pa]");
}


// main method
double
InterfrostDthetaDpressureModel::DThetaDpCoef(double nl, double sl, double phi) const
{
  return beta_ * nl * phi * sl;
}

double
InterfrostDthetaDpressureModel::DDThetaDpCoefDMolarDensityLiquid(double nl,
                                                                 double sl,
                                                                 double phi) const
{
  return beta_ * phi * sl;
}

double
InterfrostDthetaDpressureModel::DDThetaDpCoefDSaturationLiquid(double nl,
                                                               double sl,
                                                               double phi) const
{
  return beta_ * nl * phi;
}

double
InterfrostDthetaDpressureModel::DDThetaDpCoefDPorosity(double nl, double sl, double phi) const
{
  return beta_ * nl * sl;
}

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
