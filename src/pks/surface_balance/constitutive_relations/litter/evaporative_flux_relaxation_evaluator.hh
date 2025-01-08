/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The evaporative flux relaxation evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    myKeyFirst = evaporative
    evalNameString = evaporative flux relaxation
    evalNameCaps = EVAPORATIVE_FLUX_RELAXATION
    namespaceCaps = SURFACEBALANCE
    evalClassName = EvaporativeFluxRelaxation
    namespace = SurfaceBalance
    myMethodDeclarationArgs = double wc, double rho, double L
    myKey = evaporative_flux
    evalName = evaporative_flux_relaxation
    modelMethodDeclaration =   double EvaporativeFlux(double wc, double rho, double L) const;
    myKeyMethod = EvaporativeFlux
    myMethodArgs = wc_v[0][i], rho_v[0][i], L_v[0][i]

*/

#ifndef AMANZI_SURFACEBALANCE_EVAPORATIVE_FLUX_RELAXATION_EVALUATOR_HH_
#define AMANZI_SURFACEBALANCE_EVAPORATIVE_FLUX_RELAXATION_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class EvaporativeFluxRelaxationModel;

class EvaporativeFluxRelaxationEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit EvaporativeFluxRelaxationEvaluator(Teuchos::ParameterList& plist);
  EvaporativeFluxRelaxationEvaluator(const EvaporativeFluxRelaxationEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<EvaporativeFluxRelaxationModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

 protected:
  Key wc_key_;
  Key rho_key_;
  Key thickness_key_;
  Key cv_key_;

  Teuchos::RCP<EvaporativeFluxRelaxationModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaporativeFluxRelaxationEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi

#endif
