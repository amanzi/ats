/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Energy content evaluator for a standard Richards equation, including energy in the gas phase.
/*!

Calculates energy, in [KJ], via the equation:

.. math::
  E = V * ( \phi (u_l s_l n_l + u_g s_g n_g)  + (1-\phi_0) u_r \rho_r )

`"evaluator type`" = `"liquid+gas energy`"

Note this equation assumes that porosity is compressible, but is based on the
uncompressed rock grain density (not bulk density).  This means that porosity
is the base, uncompressible value when used with the energy in the grain, but
the larger, compressible value when used with the energy in the water.

.. _evaluator-liquid-gas-energy-spec:
.. admonition:: evaluator-liquid-gas-energy-spec

   DEPENDENCIES:

   - `"porosity`" The porosity, including any compressibility. [-]
   - `"base porosity`" The uncompressed porosity (note this may be the same as
     porosity for incompressible cases) [-]
   - `"molar density liquid`" [mol m^-3]
   - `"saturation liquid`" [-]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density gas`" [mol m^-3]
   - `"saturation gas`" [-]
   - `"internal energy gas`" [KJ mol^-1]
   - `"density rock`" Units may be either [kg m^-3] or [mol m^-3]
   - `"internal energy rock`" Units may be either [KJ kg^-1] or [KJ mol^-1],
     but must be consistent with the above density.
   - `"cell volume`" [m^3]

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {
namespace Relations {

class LiquidGasEnergyModel;

class LiquidGasEnergyEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit LiquidGasEnergyEvaluator(Teuchos::ParameterList& plist);
  LiquidGasEnergyEvaluator(const LiquidGasEnergyEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<LiquidGasEnergyModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

 protected:
  Key phi_key_;
  Key phi0_key_;
  Key sl_key_;
  Key nl_key_;
  Key ul_key_;
  Key sg_key_;
  Key ng_key_;
  Key ug_key_;
  Key rho_r_key_;
  Key ur_key_;
  Key cv_key_;

  Teuchos::RCP<LiquidGasEnergyModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, LiquidGasEnergyEvaluator> reg_;
};

} // namespace Relations
} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
