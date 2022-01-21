/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The carbon decompostion rate evaluator gets the subsurface temperature and pressure.
  Computes(integrates) CO2 decomposition rate.
  This is EvaluatorSecondaryMonotypeCV and depends on the subsurface temperature and pressure,

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_CARBONDECOM_EVALUATOR_
#define AMANZI_FLOWRELATIONS_CARBONDECOM_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class CarbonDecomposeRateEvaluator : public EvaluatorSecondaryMonotypeCV {

public:
  explicit
  CarbonDecomposeRateEvaluator(Teuchos::ParameterList& plist);
  CarbonDecomposeRateEvaluator(const CarbonDecomposeRateEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  virtual bool Update(State& S, const Key& request) override;
  virtual void EnsureCompatibility(State& S) override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
                         const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

  double Func_TempPres(double temp, double pres);

 protected:
  bool updated_once_;
  Key temp_key_, pres_key_, sat_key_, por_key_, cv_key_;
  Key domain_;
  double q10_;

 private:
  static Utils::RegisteredFactory<Evaluator,CarbonDecomposeRateEvaluator> reg_;

};

} //namespace
} //namespace

#endif
