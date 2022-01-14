/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! FractionalConductanceEvaluator: an obstruction-drag factor.
/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan
*/

/*!

This implements the term,

.. math:
    \frac{\Phi(\delta) - \Phi(\delta_d)}{\delta - \delta_d}

from Jan et al WRR 2018.
 
*/

#ifndef AMANZI_FLOWRELATIONS_FRACTIONAL_CONDUCTANCE_EVALUATOR_
#define AMANZI_FLOWRELATIONS_FRACTIONAL_CONDUCTANCE_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class FractionalConductanceEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  FractionalConductanceEvaluator(Teuchos::ParameterList& plist);
  FractionalConductanceEvaluator(const FractionalConductanceEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const;
  
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 private:
  Key mobile_depth_key_, pd_key_, vpd_key_;
  Key depr_depth_key_, delta_ex_key_, delta_max_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,FractionalConductanceEvaluator> factory_;
};

} //namespace
} //namespace
} //namespace

#endif
