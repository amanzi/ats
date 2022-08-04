/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
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

.. _overland-conductivity-subgrid-evaluator-spec
.. admonition:: overland-conductivity-subgrid-evaluator-spec

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
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel;

class OverlandConductivitySubgridEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  OverlandConductivitySubgridEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivitySubgridEvaluator(const OverlandConductivitySubgridEvaluator& other) = default;
  Teuchos::RCP<FieldEvaluator> Clone() const override;

  Teuchos::RCP<ManningConductivityModel> get_Model() { return model_; }
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 protected:

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

private:
  Teuchos::RCP<ManningConductivityModel> model_;

  Key slope_key_;
  Key coef_key_;
  Key dens_key_;
  Key mobile_depth_key_;
  Key drag_exp_key_;
  Key frac_cond_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,OverlandConductivitySubgridEvaluator> factory_;
};

} //namespace
} //namespace

