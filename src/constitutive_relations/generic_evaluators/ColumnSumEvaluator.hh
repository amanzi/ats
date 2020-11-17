/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ColumnSumEvaluator is the generic evaluator for summuation of a column field.
  Return summation is put on the corresponding surface cell
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_COLUMNSUM_EVALUATOR_
#define AMANZI_FLOWRELATIONS_COLUMNSUM_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class ColumnSumEvaluator : public SecondaryVariableFieldEvaluator {

public:
  explicit
  ColumnSumEvaluator(Teuchos::ParameterList& plist);
  ColumnSumEvaluator(const ColumnSumEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S,
                                      Key request);
protected:
  double coef_;

  Key dep_key_;
  Key cv_key_;
  Key molar_dens_key_;
  Key surf_cv_key_;

  Key domain_;
  Key surf_domain_;

  bool updated_once_;
private:
  static Utils::RegisteredFactory<FieldEvaluator,ColumnSumEvaluator> factory_;

};

} //namespace
} //namespace

#endif
