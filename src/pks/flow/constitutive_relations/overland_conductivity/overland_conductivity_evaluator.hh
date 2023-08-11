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
   k = n \delta K(\delta, |\nabla z|, n_{mann})

This writes the conductivity as a function of water density, the mobile depth,
or portion of the ponded depth that is mobile :math:`\delta`, and Manning's n,
and the slope magnitude.  Note that the flow law is specified as a model, which
may have multiple choices of implementation.

type : `"overland conductivity`"

.. _overland-conductivity-evaluator-spec
.. admonition:: overland-conductivity-evaluator-spec

   KEYS:
   - `"molar density liquid`" **DOMAIN-molar_density_liquid**
   - `"mobile depth`" **DOMAIN-mobile_depth** Depth of the mobile water; delta
     in the above equation.
   - `"Manning coefficient`" **DOMAIN-manning_coefficient** Surface roughness/shape
     coefficient; n_{mann} above.
   - `"slope magnitude`" **DOMAIN-slope_magnitude** Magnitude of the bed surface driving
     flow; | \nabla z | above.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ManningConductivityModel;

class OverlandConductivityEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  OverlandConductivityEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivityEvaluator(const OverlandConductivityEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  static const std::string name;
  virtual std::string getType() const override { return name; }

 private:
  Teuchos::RCP<ManningConductivityModel> model_;

  Key mobile_depth_key_;
  Key slope_key_;
  Key coef_key_;
  Key dens_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, OverlandConductivityEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
