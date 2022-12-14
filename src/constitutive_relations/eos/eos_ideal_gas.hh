/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas.

  Authors: Ethan Coon (ecoon@lanl.gov)
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
  virtual double DMolarDensityDC(std::vector<double>& params) override { return 0.; }

  virtual bool IsTemperature() override { return true; }
  virtual bool IsPressure() override { return true; }
  virtual bool IsConcentration() override { return false; }

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
