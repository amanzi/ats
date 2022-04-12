/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_OVERLAND_HEAD_NOT_FRESH_WATER_CONTENT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_OVERLAND_HEAD_NOT_FRESH_WATER_CONTENT_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class OverlandPressureUnfrozenWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  OverlandPressureUnfrozenWaterContentEvaluator(Teuchos::ParameterList& plist);
  OverlandPressureUnfrozenWaterContentEvaluator(const OverlandPressureUnfrozenWaterContentEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  void InitializeFromPlist_();

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  Key pres_key_, cv_key_, mass_dens_key_, molar_dens_key_;

  bool bar_;  // bar'd variable indicates this is potentially negative for
              // pressures less than atmospheric
  double rollover_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,OverlandPressureUnfrozenWaterContentEvaluator> reg_;

};

} //namespace
} //namespace

#endif
