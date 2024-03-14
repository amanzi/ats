/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  Simple EOS for constant density.
  Defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

*/

#ifndef AMANZI_RELATIONS_EOS_CONSTANT_HH_
#define AMANZI_RELATIONS_EOS_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSConstant : public EOSConstantMolarMass {
 public:
  explicit EOSConstant(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(std::vector<double>& params) override { return rho_; }

  virtual double DMolarDensityDT(std::vector<double>& params) override { return 0.; }

  virtual double DMolarDensityDp(std::vector<double>& params) override { return 0.; }

  virtual double DMolarDensityDC(std::vector<double>& params) override { return 0.; }

  virtual bool IsTemperature() override { return false; }
  virtual bool IsPressure() override { return false; }
  virtual bool IsConcentration() override { return false; }

 private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double rho_;

  static Utils::RegisteredFactory<EOS, EOSConstant> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
