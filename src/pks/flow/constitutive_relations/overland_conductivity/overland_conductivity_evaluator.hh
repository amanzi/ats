/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
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
converts the flow law to mass flux rather than volumetric flux.

Also, this evaluator can be used in snow redistribution, and in that case needs
some extra factors (timestep size) to ensure the correct flow law in that case.

.. _overland-conductivity-evaluator-spec
.. admonition:: overland-conductivity-evaluator-spec

   * `"include density`" ``[bool]`` **true** Include the density prefactor,
     converting the flux from volumetric flux to mass flux.
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
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel;

class OverlandConductivityEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  OverlandConductivityEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivityEvaluator(const OverlandConductivityEvaluator& other) = default;
  Teuchos::RCP<FieldEvaluator> Clone() const override;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

  Teuchos::RCP<ManningConductivityModel> get_Model() { return model_; }

protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

private:
  Teuchos::RCP<ManningConductivityModel> model_;

  Key mobile_depth_key_;
  Key slope_key_;
  Key coef_key_;
  Key dens_key_;
  double dt_swe_factor_;
  bool dens_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,OverlandConductivityEvaluator> factory_;
};

} //namespace
} //namespace

