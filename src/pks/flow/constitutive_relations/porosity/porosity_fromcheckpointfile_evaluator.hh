/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Takes porosity from checkpoint file.

  Compressible grains are both physically realistic (based on bulk modulus)
  and a simple way to provide a non-elliptic, diagonal term for helping
  solvers to converge.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_POROSITY_CHECKPOINT_FILE_EVALUATOR_HH_
#define AMANZI_FLOWRELATIONS_POROSITY_CHECKPOINT_FILE_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class PorosityFromCheckpointFileEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit
  PorosityFromCheckpointFileEvaluator(Teuchos::ParameterList& plist);
  PorosityFromCheckpointFileEvaluator(const PorosityFromCheckpointFileEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  //Teuchos::RCP<PorosityFromCheckpointFileModelPartition> get_Models() { return models_; }

protected:
  Key poro_key_;
  //Key pres_key_;

  //Teuchos::RCP<PorosityFromCheckpointFileModelPartition> models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,PorosityFromCheckpointFileEvaluator> fac_;



};

} // namespace
} // namespace

#endif
