/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Ahmad Jan
*/

//! Evaluates the conductivity of surface flow in the flow subgrid model.
/*!

This implements the conductivity of overland flow in the subgrid model case
from Jan et al WRR 2018.  This calculates the nonlinear coefficient in the
diffusion wave equation.  The term is given by:

.. math:
   k = n_l K^\beta \frac{\delta^{coef}}{n_{mann} \sqrt(| \nabla z |)}

This includes a density factor, typically a molar density, which
converts the flow law to water flux rather than volumetric flux.

.. _overland_conductivity_subgrid_evaluator-spec
.. admonition:: overland_conductivity_subgrid_evaluator-spec

   DEPENDENCIES:
   - `"mobile depth`" **DOMAIN-mobile_depth** Depth of the mobile water; delta
     in the above equation.
   - `"slope`" **DOMAIN-slope_magnitude** Magnitude of the bed surface driving
     flow; | \nabla z | above.
   - `"coefficient`" **DOMAIN-manning_coefficient** Surface roughness/shape
     coefficient; n_{mann} above.
   - `"molar density liquid`" **DOMAIN-molar_density_liquid** If `"include
     density`" is true, the density.
   - `"fractional conductance`" **DOMAIN-fractional_conductance** The leading
     conductance term, K in the above equation.
   - `"drag exponent`" **DOMAIN-drag_exponent** Power law for the fractional
     conductance, \beta above.

*/
#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel;

class OverlandConductivitySubgridEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  OverlandConductivitySubgridEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivitySubgridEvaluator(const OverlandConductivitySubgridEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<ManningConductivityModel> get_Model() { return model_; }

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 private:
  Teuchos::RCP<ManningConductivityModel> model_;

  Key slope_key_;
  Key coef_key_;
  Key dens_key_;
  Key mobile_depth_key_;
  Key drag_exp_key_;
  Key frac_cond_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, OverlandConductivitySubgridEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi
