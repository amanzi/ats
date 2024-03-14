/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates the conductivity of surface flow.
/*!

This implements the conductivity of overland flow, which is the nonlinear
coefficient in the diffusion wave equation.  The term is given by:

.. math:
   k = \frac{\delta^{coef}}{n_{mann} \sqrt(| \nabla z |)}

Optionally, this may include a density factor, typically a molar density, which
converts the flow law to water flux rather than volumetric flux.

Also, this evaluator can be used in snow redistribution, and in that case needs
some extra factors (timestep size) to ensure the correct flow law in that case.

`"evaluator type`" = `"overland conductivity`"

.. _overland-conductivity-evaluator-spec:
.. admonition:: overland-conductivity-evaluator-spec

   * `"include density`" ``[bool]`` **true** Include the density prefactor,
     converting the flux from volumetric flux to water flux.
   * `"dt factor [s]`" ``[double]`` **-1** The artificial timestep size used in calculating
      snow redistribution, only used in that case.
   * `"swe density factor [-]`" ``[double]`` **10** Ratio of water to snow density.

   DEPENDENCIES:

   - `"mobile depth`" **DOMAIN-mobile_depth** Depth of the mobile water; delta
     in the above equation.
   - `"slope`" **DOMAIN-slope_magnitude** Magnitude of the bed surface driving
     flow; | \nabla z | above.
   - `"coefficient`" **DOMAIN-manning_coefficient** Surface roughness/shape
     coefficient; n_{mann} above.
   - `"molar density liquid`" **DOMAIN-molar_density_liquid** If `"include
     density`" is true, the density.

*/
#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel;

class OverlandConductivityEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  OverlandConductivityEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivityEvaluator(const OverlandConductivityEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<ManningConductivityModel> get_Model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 private:
  Teuchos::RCP<ManningConductivityModel> model_;

  Key mobile_depth_key_;
  Key slope_key_;
  Key coef_key_;
  Key dens_key_;
  double dt_swe_factor_;
  bool dens_;

 private:
  static Utils::RegisteredFactory<Evaluator, OverlandConductivityEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi
