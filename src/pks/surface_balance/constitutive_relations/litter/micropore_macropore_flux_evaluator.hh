/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Daniil Svyatsky (dasvyat@lanl.gov)
*/

//! Exchange flux between multiple continua.
/*!

Evaluates the following exchange flux model:

.. math::
   q_{exchange} = n_l k_r K \frac{\Gamma}{\delta} (p_M - p_m)

where :math:`p` is the pressure of the Macro and micro porespaces,
respectively, K is some measure of an absolute permeability, :math:`\Gamma [-]`
is the exchange coefficient, :math:`\delta [m]` is a unit of distance
characterizing the typical distance between pore, and :math:`k_r` is the
relative permeability, which is upwinded based on the larger of the two
pressures, 'n_l' is molar density liquid

Note that the expected domain for this is the micropore domain, but may be
changed on the input line.

.. _evaluator-micropore-macropore-flux-spec:
.. admonition:: evaluator-micropore-macropore-flux-spec

   * `"micropore domain`" ``[string]`` **DOMAIN** Defaults to the domain of the flux's
     variable name.
   * `"macropore domain`" ``[string]`` **MACROPORE_DOMAIN** Guesses based off of DOMAIN

   * `"micropore macropore flux model parameters`" ``[micropore-macropore-flux-model-spec]``

   KEYS:

   - `"micropore pressure`" **pressure**
   - `"macropore pressure`" **MACROPORE_DOMAIN-pressure**
   - `"micropore relative permeability`" **relative_permeability**
   - `"macropore relative permeability`" **MACROPORE_DOMAIN-relative_permeability**
   - `"permeability`" **permeability**
   - `"micropore molar density liquid`" **molar_density_liquid**

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class MicroporeMacroporeFluxModel;

class MicroporeMacroporeFluxEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit MicroporeMacroporeFluxEvaluator(Teuchos::ParameterList& plist);
  MicroporeMacroporeFluxEvaluator(const MicroporeMacroporeFluxEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<MicroporeMacroporeFluxModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

 protected:
  Key pm_key_;
  Key pM_key_;
  Key krM_key_;
  Key krm_key_;
  Key K_key_;
  Key den_key_;

  Teuchos::RCP<MicroporeMacroporeFluxModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, MicroporeMacroporeFluxEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
