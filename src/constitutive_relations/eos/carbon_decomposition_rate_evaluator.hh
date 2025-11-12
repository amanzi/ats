/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

/*
  The carbon decompostion rate evaluator gets the subsurface temperature and pressure.
  Computes(integrates) CO2 decomposition rate.
  This is EvaluatorSecondaryMonotypeCV and depends on the subsurface temperature and pressure,

*/

#ifndef AMANZI_FLOWRELATIONS_CARBONDECOM_EVALUATOR_
#define AMANZI_FLOWRELATIONS_CARBONDECOM_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class CarbonDecomposeRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit CarbonDecomposeRateEvaluator(Teuchos::ParameterList& plist);
  CarbonDecomposeRateEvaluator(const CarbonDecomposeRateEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  bool IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return false;
  }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override
  {}

  double Func_Pres(double pres) const;
  double Func_Temp(double temp, double q10) const;
  double Func_Depth(double depth) const;
  double Func_Satliq(double sat) const;
  double Func_Satgas(double sat) const;

 protected:
  Key temp_key_;
  Key pres_key_;
  Key sat_liq_key_;
  Key sat_gas_key_;
  Key por_key_;
  Key depth_key_;
  Key domain_;
  Key domain_surf_;

  double q10_;
  bool is_func_temp_, is_func_depth_, is_func_liq_, is_func_gas_;
  bool is_func_pres_, is_scaling_down_, is_threshold_temp_;

 private:
  static Utils::RegisteredFactory<Evaluator, CarbonDecomposeRateEvaluator> reg_;
};

} // namespace Relations
} // namespace Amanzi

#endif
