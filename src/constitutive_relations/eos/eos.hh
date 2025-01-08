/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  EOS -- purely virtual base class for an EOS.
  std::vector<double> params contains parameters which define EOS.

*/

#ifndef AMANZI_RELATIONS_EOS_HH_
#define AMANZI_RELATIONS_EOS_HH_

#include <vector>

namespace Amanzi {
namespace Relations {

class EOS {
 public:
  virtual ~EOS(){};

  // Virtual methods that form the EOS
  virtual double MassDensity(std::vector<double>& params) = 0;
  virtual double DMassDensityDT(std::vector<double>& params) { return 0.; }
  virtual double DMassDensityDp(std::vector<double>& params) { return 0.; }
  virtual double DMassDensityDC(std::vector<double>& params) { return 0.; }

  virtual double MolarDensity(std::vector<double>& params) = 0;
  virtual double DMolarDensityDT(std::vector<double>& params) { return 0.; }
  virtual double DMolarDensityDp(std::vector<double>& params) { return 0.; }
  virtual double DMolarDensityDC(std::vector<double>& params) { return 0.; }

  // If molar mass is constant, we can take some shortcuts if we need both
  // molar and mass densities.  MolarMass() is undefined if
  // !IsConstantMolarMass()
  virtual bool IsConstantMolarMass() = 0;
  virtual double MolarMass() = 0;
  virtual bool IsTemperature() = 0;
  virtual bool IsPressure() = 0;
  virtual bool IsConcentration() = 0;
};

} // namespace Relations
} // namespace Amanzi

#endif
