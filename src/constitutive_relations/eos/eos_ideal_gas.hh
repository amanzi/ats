/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  EOS for an ideal gas.

*/

#ifndef AMANZI_RELATIONS_EOS_IDEAL_GAS_HH_
#define AMANZI_RELATIONS_EOS_IDEAL_GAS_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSIdealGas : public EOSConstantMolarMass {
 public:
  explicit EOSIdealGas(Teuchos::ParameterList& eos_plist);

  virtual double MolarDensity(std::vector<double>& params) override;
  virtual double DMolarDensityDT(std::vector<double>& params) override;
  virtual double DMolarDensityDp(std::vector<double>& params) override;
  virtual double DMolarDensityDMolarRatio(std::vector<double>& params) override { return 0.; }

  virtual bool IsTemperature() override { return true; }
  virtual bool IsPressure() override { return true; }
  virtual bool IsMolarRatio() override { return false; }

 protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double R_;

 private:
  static Utils::RegisteredFactory<EOS, EOSIdealGas> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
