/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  The trapping evaluator gets the trapping rates.
  we have expressed the trapping rate using the
  approach proposed by Palmer et al. [2004],

  ..math::
  Q_{tr} = \rho_s\,C\,U\,\epsilon d_s n_s min(h_s, p_D)

  ..math::
  \epsilon = \alpha_e\left(\frac{Ud_s}{\nu}\right)^{\beta_\epsilon}
             \left(\frac{d}{d_s}\right)^{\gamma_\epsilon}

  where math::'\rho_s' is the sediment density, math::'C' in a sediment molar ratio,
  math::'U' is the absolute value of the local flow speed, math::'\epsilon' is
  a capture efficiency which gives the rate at which
  transported sediment particles are captured by plant stems,
  math::'d' is a particle diameter, math::'\nu' is the kinematic viscosity,
  math::'d_s' is the stem diameter, math::'n_s' is the stem density,
  math::'h_s' is the average height of the stems, and math::'p_D" is the local flow
  ponded-depth

  `"evaluator type`" = `"trapping rate`"

   * `"kinematic viscosity`" ``[double]`` **1e-6**
   * `"particle diameter`" ``[double]`` **5.e-5**
   * `"alpha`" ``[double]`` **0.224**
   * `"beta`" ``[double]`` **0.718**
   * `"gamma`" ``[double]`` **2.08**

     DEPENDENCIES:

   - `"velocity`" **SURFACE-DOMAIN-velocity**
   - `"sediment`" **SURFACE-DOMAIN-sediment**
   - `"stem_diameter`" **SURFACE-DOMAIN-stem_diameter**
   - `"stem_height`" **SURFACE-DOMAIN-stem_height**
   - `"stem_density`" **SURFACE-DOMAIN-stem_density**
   - `"ponded depth`" **SURFACE-DOMAIN-ponded_depth**

*/
#pragma once

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class TrappingRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit TrappingRateEvaluator(Teuchos::ParameterList& plist);

  TrappingRateEvaluator(const TrappingRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  double visc_, d_p_, alpha_, beta_, gamma_;

  Key velocity_key_;
  Key sediment_key_;
  Key ponded_depth_key_;
  Key stem_diameter_key_, stem_height_key_, stem_density_key_;
  double sediment_density_;

  static Utils::RegisteredFactory<Evaluator, TrappingRateEvaluator> factory_;
};

} // namespace Amanzi
