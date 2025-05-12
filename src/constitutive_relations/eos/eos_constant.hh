/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

A constant molar and mass density of water, independent of temperature or
concentration, but related to each other by a molar mass of water.

Note, users should prefer to use `"independent variable constant`" to this.

`"EOS type`" = `"constant`"

.. _eos-constant-spec:
.. admonition:: eos-constant-spec

   ONE OF

   * `"molar mass [kg mol^-1]`" ``[double]`` **0.0180153**

   OR

   * `"molar mass [g mol^-1]`" ``[double]`` **18.0153**

   END

   ONE OF

   * `"density [mol m^-3]`" ``[double]`` molar density of water

   OR

   * `"density [kg m^-3]`" ``[double]`` mass density of water

   END



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

  virtual double DMolarDensityDMolarRatio(std::vector<double>& params) override { return 0.; }

  virtual bool IsTemperature() override { return false; }
  virtual bool IsPressure() override { return false; }
  virtual bool IsMolarRatio() override { return false; }

 private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double rho_;

  static Utils::RegisteredFactory<EOS, EOSConstant> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
