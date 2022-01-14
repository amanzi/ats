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
  CarbonDecomposeRateEvaluator(const CarbonDecomposeRateEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;
  
protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);
  
    
  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);
  
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  double  Func_TempPres(double temp, double pres);

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
