/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*

A molar and mass density of water that are linear in pressure.

.. math::

   n &= n_0 + \beta (p - p_{atm})
   \rho &= M n

`"EOS type`" = `"linear`"

.. _eos-linear-spec:
.. admonition:: eos-linear-spec

   ONE OF

   * `"molar mass [kg mol^-1]`" ``[double]`` **0.0180153** :math:`M` above

   OR

   * `"molar mass [g mol^-1]`" ``[double]`` **18.0153**

   END

   ONE OF

   * `"density [mol m^-3]`" ``[double]`` molar density of water, :math:`n_0` above

   OR

   * `"density [kg m^-3]`" ``[double]`` mass density of water

   END

   * `"compressibility [Pa^-1]`" ``[double]`` :math:`\beta` above.

*/

#ifndef AMANZI_RELATIONS_EOS_LINEAR_HH_
#define AMANZI_RELATIONS_EOS_LINEAR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace ATS_Physics {
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
  virtual double DMassDensityDMoleFraction(std::vector<double>& params) override { return 0.; }

  virtual bool IsTemperature() override { return false; }
  virtual bool IsPressure() override { return true; }
  virtual bool IsMoleFraction() override { return false; }

 private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double rho_;
  double beta_;

  static Utils::RegisteredFactory<EOS, EOSLinear> factory_;
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi

#endif
