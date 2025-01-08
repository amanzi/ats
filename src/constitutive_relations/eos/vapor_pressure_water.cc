/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  ATS

  Saturated Vapor Pressure for vapor over water or ice, Sonntag (1990)

*/

#include <cmath>
#include "errors.hh"
#include "vapor_pressure_water.hh"

namespace Amanzi {
namespace Relations {

VaporPressureWater::VaporPressureWater(Teuchos::ParameterList& plist)
  : plist_(plist),
    ka0_(16.635764),
    ka_(-6096.9385),
    kb_(-2.7111933e-2),
    kc_(1.673952e-5),
    kd_(2.433502)
{}

double
VaporPressureWater::SaturatedVaporPressure(double T)
{
  if (T < 100. || T > 373.0) {
    std::cout << "Invalid temperature, T = " << T << std::endl;
    Exceptions::amanzi_throw(Errors::CutTimestep());
  }
  return 100.0 * exp(ka0_ + ka_ / T + (kb_ + kc_ * T) * T + kd_ * log(T));
};

double
VaporPressureWater::DSaturatedVaporPressureDT(double T)
{
  if (T < 100. || T > 373.0) {
    std::cout << "Invalid temperature, T = " << T << std::endl;
    Errors::Message m("Cut timestep");
    Exceptions::amanzi_throw(m);
  }
  return SaturatedVaporPressure(T) * (-ka_ / (T * T) + kb_ + 2.0 * kc_ * T + kd_ / T);
};


} // namespace Relations
} // namespace Amanzi
