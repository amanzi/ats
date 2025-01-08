/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Energy content evaluator for a three-phase, gas, liquid, ice system including the surrounding soil.
/*!

Calculates energy, in [KJ], via the equation:

.. math::
  E = V * ( \phi (u_l s_l n_l + u_i s_i n_i + u_g s_g n_g)  + (1-\phi_0) u_r \rho_r )

Specified with evaluator type: `"three phase energy`"

Note this equation assumes that porosity is compressible, but is based on the
uncompressed rock grain density (not bulk density).  This means that porosity
is the base, uncompressible value when used with the energy in the grain, but
the larger, compressible value when used with the energy in the water.

.. _field_evaluator_type_three_phase_energy-spec:
.. admonition:: field_evaluator_type_three_phase_energy-spec

   DEPENDENCIES:

   - `"porosity`" The porosity, including any compressibility. [-]
   - `"base porosity`" The uncompressed porosity (note this may be the same as
     porosity for incompressible cases) [-]
   - `"molar density liquid`" [mol m^-3]
   - `"saturation liquid`" [-]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density ice`" [mol m^-3]
   - `"saturation ice`" [-]
   - `"internal energy ice`" [KJ mol^-1]
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
namespace Energy {
namespace Relations {

class ThreePhaseEnergyModel;

class ThreePhaseEnergyEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit ThreePhaseEnergyEvaluator(Teuchos::ParameterList& plist);
  ThreePhaseEnergyEvaluator(const ThreePhaseEnergyEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<ThreePhaseEnergyModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

  Key phi_key_;
  Key phi0_key_;
  Key sl_key_;
  Key nl_key_;
  Key ul_key_;
  Key si_key_;
  Key ni_key_;
  Key ui_key_;
  Key sg_key_;
  Key ng_key_;
  Key ug_key_;
  Key rho_r_key_;
  Key ur_key_;
  Key cv_key_;

  Teuchos::RCP<ThreePhaseEnergyModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, ThreePhaseEnergyEvaluator> reg_;
};

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
