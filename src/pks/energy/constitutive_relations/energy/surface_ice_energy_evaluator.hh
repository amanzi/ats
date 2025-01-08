/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Energy content for a surface water, partially frozen system.
/*!

The energy associated with ponded water, in [KJ], given by:

.. math::
  E = V * ( \eta h u_l n_l + (1 - \eta) h u_i n_i )

Specified with evaluator type: `"surface ice energy`"

.. _field_evaluator_type_surface_ice_energy-spec:
.. admonition:: field_evaluator_type_surface_ice_energy-spec

   DEPENDENCIES:

   - `"ponded depth`"  Height of water above the land surface [m]
   - `"unfrozen fraction`"  The fraction of unfrozen water ranges from 0 to 1. [-]
   - `"molar density liquid`" [mol m^-3]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density ice`" [mol m^-3]
   - `"internal energy ice`" [KJ mol^-1]
   - `"cell volume`" [m^2]

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

class SurfaceIceEnergyModel;

class SurfaceIceEnergyEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SurfaceIceEnergyEvaluator(Teuchos::ParameterList& plist);
  SurfaceIceEnergyEvaluator(const SurfaceIceEnergyEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<SurfaceIceEnergyModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

  Key h_key_;
  Key eta_key_;
  Key nl_key_;
  Key ul_key_;
  Key ni_key_;
  Key ui_key_;
  Key cv_key_;

  Teuchos::RCP<SurfaceIceEnergyModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, SurfaceIceEnergyEvaluator> reg_;
};

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
