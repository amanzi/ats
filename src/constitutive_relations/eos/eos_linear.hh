/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  Simple EOS for compressibility with pressure.

  http://software.lanl.gov/ats/trac

*/

#ifndef AMANZI_RELATIONS_EOS_LINEAR_HH_
#define AMANZI_RELATIONS_EOS_LINEAR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSLinear : public EOSConstantMolarMass {
 public:
  explicit EOSLinear(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(std::vector<double>& params) override
  {
    return rho_ * (1 + beta_ * std::max(params[0] - 101325., 0.));
  }
  virtual double DMassDensityDp(std::vector<double>& params) override
  {
    return params[0] > 101325. ? rho_ * beta_ : 0.;
  }
  virtual double DMassDensityDT(std::vector<double>& params) override { return 0.; }
  virtual double DMassDensityDMolarRatio(std::vector<double>& params) override { return 0.; }

  virtual bool IsTemperature() override { return false; }
  virtual bool IsPressure() override { return true; }
  virtual bool IsMolarRatio() override { return false; }

 private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double rho_;
  double beta_;

  static Utils::RegisteredFactory<EOS, EOSLinear> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
