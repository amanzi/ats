/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

This evaluator computes actual evaporation from potential evaporation due to a
dessicated zone via soil resistance.

.. math::

   E = \frac{E_{potential}}{1 + R_{soil}}

`"evaluator type`" = `"evaporation downregulation, soil resistance`"

.. _evaluator-evaporation-downregulation-soil-resistance:
.. admonition:: evaluator-evaporation-downregulation-soil-resistance

   DEPENDENCIES:

   - `"potential evaporation`"
   - `"soil resistance`"

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {


class EvaporationDownregulationEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit EvaporationDownregulationEvaluator(Teuchos::ParameterList& plist);
  EvaporationDownregulationEvaluator(const EvaporationDownregulationEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    // this will mostly be differentiated with respect to pressure for flow
    // Jacobians, but none of the terms that _really_ depend on p are actually
    // implemented.  That would require differentiating RSoil with respect to
    // s_l, s_g, etc.  But only derivatives wrt potential evaporation are
    // implemented.  That will rarely if ever be p-dependent.  Therefore, this
    // is just turned off to avoid lengthy calculations with 0.
    return false;
  }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key rsoil_key_;
  Key pot_evap_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaporationDownregulationEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
