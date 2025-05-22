/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

A molar and mass density of gas including water vapor.  This wraps another EOS
for providing molar density, but removes mass density functions as they are not
valid for this use case.  Typically this wraps the ideal gas law EOS.

`"EOS type`" = `"vapor in gas`"

.. _eos-vapor-in-gas-spec:
.. admonition:: eos-vapor-in-gas-spec

   * `"gas EOS parameters`" ``[eos-typedinline-spec]``

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
  double DMassDensityDMoleFraction(std::vector<double>& params) override
  {
    AMANZI_ASSERT(0);
    return 0.0;
  }

  double MolarDensity(std::vector<double>& params) override;
  double DMolarDensityDT(std::vector<double>& params) override;
  double DMolarDensityDp(std::vector<double>& params) override;
  double DMolarDensityDMoleFraction(std::vector<double>& params) override { return 0.; }

  virtual bool IsTemperature() override { return gas_eos_->IsTemperature(); }
  virtual bool IsPressure() override { return gas_eos_->IsPressure(); }
  virtual bool IsMoleFraction() override { return gas_eos_->IsMoleFraction(); }

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
