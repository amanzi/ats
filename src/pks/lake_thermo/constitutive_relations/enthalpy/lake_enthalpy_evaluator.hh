/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for enthalpy.
----------------------------------------------------------------------------- */


#ifndef AMANZI_LAKE_ENTHALPY_EVALUATOR_HH_
#define AMANZI_LAKE_ENTHALPY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

class LakeEnthalpyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LakeEnthalpyEvaluator(Teuchos::ParameterList& plist);
  LakeEnthalpyEvaluator(const LakeEnthalpyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  Key pres_key_;
  Key dens_key_;
  Key ie_key_;
  bool include_work_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LakeEnthalpyEvaluator> factory_;

};

} // namespace
} // namespace

#endif
