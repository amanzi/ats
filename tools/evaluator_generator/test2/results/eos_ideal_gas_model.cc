/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The ideal gas equation of state model is an algebraic model with dependencies.

  Generated via evaluator_generator with:


*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "eos_ideal_gas_model.hh"

namespace Amanzi {
namespace General {
namespace Relations {

// Constructor from ParameterList
EosIdealGasModel::EosIdealGasModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
EosIdealGasModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  cv_ = plist.get<double>("heat capacity");
  T0_ = plist.get<double>("reference temperature [K]");
}


// main method
double
EosIdealGasModel::Density(double temp, double pres) const
{
  return cv_ * (-T0_ + temp);
}

double
EosIdealGasModel::DDensityDTemperature(double temp, double pres) const
{
  return cv_;
}

double
EosIdealGasModel::DDensityDPressure(double temp, double pres) const
{
  return 0;
}

} // namespace Relations
} // namespace General
} // namespace Amanzi
