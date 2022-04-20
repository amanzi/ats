/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
//! Evaluates the porosity, given a small compressibility of rock.
/*!

Compressible grains are both physically realistic (based on bulk modulus) and a
simple way to provide a non-elliptic, diagonal term for helping solvers to
converge.

* `"compressible porosity model parameters`" ``[compressible-porosity-leijnse-model-spec]``

KEYS:
- `"pressure`" **DOMAIN-pressure**
- `"base porosity`" **DOMAIN-base_porosity**

*/

#ifndef AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_LEIJNSE_EVALUATOR_HH_
#define AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_LEIJNSE_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "compressible_porosity_leijnse_model_partition.hh"

namespace Amanzi {
namespace Flow {

class CompressiblePorosityLeijnseEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit
  CompressiblePorosityLeijnseEvaluator(Teuchos::ParameterList& plist);
  CompressiblePorosityLeijnseEvaluator(const CompressiblePorosityLeijnseEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<CompressiblePorosityLeijnseModelPartition> get_Models() { return models_; }

protected:
  Key poro_key_;
  Key pres_key_;

  Teuchos::RCP<CompressiblePorosityLeijnseModelPartition> models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,CompressiblePorosityLeijnseEvaluator> fac_;



};

} // namespace
} // namespace

#endif
