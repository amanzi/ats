/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/* -----------------------------------------------------------------------------
ATS

Evaluator for internal energy.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

IE = phi * (s_liquid * n_liquid * u_liquid + s_gas * n_gas * u_gas
                + s_ice * n_ice * u_ice)
  + (1 - phi) * rho_rock * u_rock

This is simply the conserved quantity in the energy equation.
----------------------------------------------------------------------------- */


#ifndef AMANZI_INTERFROST_ENERGY_EVALUATOR_HH_
#define AMANZI_INTERFROST_ENERGY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class InterfrostEnergyEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit InterfrostEnergyEvaluator(Teuchos::ParameterList& energy_plist);

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  double beta_;

 private:
  static Utils::RegisteredFactory<Evaluator, InterfrostEnergyEvaluator> reg_;
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
