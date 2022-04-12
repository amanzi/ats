/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the porosity, given a small compressibility of rock.

  Compressible grains are both physically realistic (based on bulk modulus)
  and a simple way to provide a non-elliptic, diagonal term for helping
  solvers to converge.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Compressible grains are both physically realistic (based on bulk modulus) and a
simple way to provide a non-elliptic, diagonal term for helping solvers to
converge.

* `"compressible porosity model parameters`" ``[compressible-porosity-model-spec-list]``

KEYS:
- `"pressure`" **DOMAIN-pressure**
- `"base porosity`" **DOMAIN-base_porosity**

*/


#ifndef AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_EVALUATOR_HH_
#define AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "compressible_porosity_model_partition.hh"

namespace Amanzi {
namespace Flow {

class CompressiblePorosityEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit
  CompressiblePorosityEvaluator(Teuchos::ParameterList& plist);
  CompressiblePorosityEvaluator(const CompressiblePorosityEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<CompressiblePorosityModelPartition> get_Models() { return models_; }

protected:
  Key poro_key_;
  Key pres_key_;

  Teuchos::RCP<CompressiblePorosityModelPartition> models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,CompressiblePorosityEvaluator> fac_;



};

} // namespace
} // namespace

#endif
