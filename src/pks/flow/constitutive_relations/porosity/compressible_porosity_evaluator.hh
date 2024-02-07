/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the porosity, given a small compressibility of rock.

  Compressible grains are both physically realistic (based on bulk modulus)
  and a simple way to provide a non-elliptic, diagonal term for helping
  solvers to converge.

*/

/*!

Compressible grains are both physically realistic (based on bulk modulus) and a
simple way to provide a non-elliptic, diagonal term for helping solvers to
converge.

`"evaluator type`" = `"compressible porosity`"

.. _compressible-porosity-evaluator-spec:
.. admonition:: compressible-porosity-evaluator-spec

   * `"compressible porosity model parameters`" ``[compressible-porosity-standard-model-spec-list]``

   KEYS:

   - `"pressure`" **DOMAIN-pressure**
   - `"base porosity`" **DOMAIN-base_porosity**

*/


#ifndef AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_EVALUATOR_HH_
#define AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "compressible_porosity_model_partition.hh"

namespace Amanzi {
namespace Flow {

class CompressiblePorosityEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit CompressiblePorosityEvaluator(Teuchos::ParameterList& plist);
  CompressiblePorosityEvaluator(const CompressiblePorosityEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<CompressiblePorosityModelPartition> get_Models() { return models_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key poro_key_;
  Key pres_key_;

  Teuchos::RCP<CompressiblePorosityModelPartition> models_;

 private:
  static Utils::RegisteredFactory<Evaluator, CompressiblePorosityEvaluator> fac_;
};

} // namespace Flow
} // namespace Amanzi

#endif
