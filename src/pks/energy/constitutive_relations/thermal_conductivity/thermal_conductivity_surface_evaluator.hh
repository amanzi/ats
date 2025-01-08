/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Interface for a thermal conductivity model with two phases.

*/

/*!

Thermal conductivity of surface water that can be either frozen or liquid phase.

`"evaluator type`" = `"surface thermal conductivity`"

.. _thermal_conductivity_surface_evaluator-spec:
.. admonition:: thermal_conductivity_surface_evaluator-spec

   * `"thermal conductivity parameters`" ``[thermal_conductivity_surface-spec]``

   KEYS:

   - `"unfrozen fraction`"
   - `"ponded depth`"

.. _thermal_conductivity_surface-spec:
.. admonition:: thermal_conductivity_surface-spec

   * `"thermal conductivity of water [W m^-1 K^-1]`" ``[double]`` **0.58**
   * `"thermal conductivity of ice [W m^-1 K^-1]`" ``[double]`` **2.18**
   * `"minimum thermal conductivity`" ``[double]`` **1.e-14**

*/

#ifndef AMANZI_ENERGY_RELATIONS_TC_SURFACE_EVALUATOR_HH_
#define AMANZI_ENERGY_RELATIONS_TC_SURFACE_EVALUATOR_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivitySurfaceEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  ThermalConductivitySurfaceEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivitySurfaceEvaluator(const ThermalConductivitySurfaceEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from SecondaryVariableFieldModel
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  // dependencies
  Key uf_key_;
  Key height_key_;

  double K_liq_;
  double K_ice_;
  double min_K_;

 private:
  static Utils::RegisteredFactory<Evaluator, ThermalConductivitySurfaceEvaluator> factory_;
};

} // namespace Energy
} // namespace Amanzi

#endif
