/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  ATS

  EOS for a combination of air and vapor pressure.  Mass density is not
  available, not because it can't be calculated, but because it depends upon
  omega.  It's not really needed, and if it were, would not fit the EOS
  interface without serious work.

*/

#ifndef AMANZI_RELATIONS_EOS_VAPOR_IN_GAS_HH_
#define AMANZI_RELATIONS_EOS_VAPOR_IN_GAS_HH_

#include "Teuchos_ParameterList.hpp"
#include "eos.hh"
#include "dbc.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSVaporInGas : public EOS {
 public:
  EOSVaporInGas(Teuchos::ParameterList& eos_plist);

  double MassDensity(std::vector<double>& params) override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }
  double DMassDensityDT(std::vector<double>& params) override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }
  double DMassDensityDp(std::vector<double>& params) override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }
  double DMassDensityDC(std::vector<double>& params) override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }

  double MolarDensity(std::vector<double>& params) override;
  double DMolarDensityDT(std::vector<double>& params) override;
  double DMolarDensityDp(std::vector<double>& params) override;
  double DMolarDensityDC(std::vector<double>& params) override { return 0.; }

  virtual bool IsTemperature() override { return gas_eos_->IsTemperature(); }
  virtual bool IsPressure() override { return gas_eos_->IsPressure(); }
  virtual bool IsConcentration() override { return gas_eos_->IsConcentration(); }

  bool IsConstantMolarMass() override { return false; }
  double MolarMass() override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }

 protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  Teuchos::RCP<EOS> gas_eos_;

 private:
  static Utils::RegisteredFactory<EOS, EOSVaporInGas> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
