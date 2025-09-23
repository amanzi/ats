/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOSConstantMolarMass -- intermediate class for default implementations of
  EOS with constant molar masses.

  Note, while instantiating this class does work, it should not be done.
  Instead, this class is intended to be inherited and either the Molar or Mass
  methods replaced.  Use as it stands results in an infinite recursion...

*/

#ifndef AMANZI_RELATIONS_EOS_CONSTANT_MM_HH_
#define AMANZI_RELATIONS_EOS_CONSTANT_MM_HH_

#include "eos.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

class EOSConstantMolarMass : public EOS {
 public:
  EOSConstantMolarMass()
    : M_(0.0)
  {}
  explicit EOSConstantMolarMass(double M)
    : M_(M)
  {}

  virtual double MolarDensity(std::vector<double>& params) { return MassDensity(params) / M_; }

  virtual double DMolarDensityDT(std::vector<double>& params)
  {
    return DMassDensityDT(params) / M_;
  }

  virtual double DMolarDensityDp(std::vector<double>& params)
  {
    return DMassDensityDp(params) / M_;
  }

  virtual double DMolarDensityDMoleFraction(std::vector<double>& params)
  {
    return DMassDensityDMoleFraction(params) / M_;
  }

  virtual double MassDensity(std::vector<double>& params) { return MolarDensity(params) * M_; }

  virtual double DMassDensityDT(std::vector<double>& params)
  {
    return DMolarDensityDT(params) * M_;
  }

  virtual double DMassDensityDp(std::vector<double>& params)
  {
    return DMolarDensityDp(params) * M_;
  }

  virtual double DMassDensityDMoleFraction(std::vector<double>& params)
  {
    return DMolarDensityDMoleFraction(params) * M_;
  }

  virtual bool IsConstantMolarMass() { return true; }
  virtual double MolarMass() { return M_; }

 protected:
  double M_;
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi

#endif
