/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Evaluates the porosity, given a small compressibility of rock.
/*!

Compressible grains are both physically realistic (based on bulk modulus) and a
simple way to provide a non-elliptic, diagonal term for helping solvers to
converge.  After Leijnse thesis, 1992.

`"evaluator type`" = `"compressible porosity leijnse`"

.. _evaluator-compressible-porosity-leijnse-spec
.. admonition:: evaluator-compressible-porosity-leijnse-spec

   * `"compressible porosity model parameters`" ``[compressible-porosity-leijnse-model-spec]``

   KEYS:

   - `"pressure`" **DOMAIN-pressure**
   - `"base porosity`" **DOMAIN-base_porosity**

*/

#ifndef AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_LEIJNSE_EVALUATOR_HH_
#define AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_LEIJNSE_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "compressible_porosity_leijnse_model_partition.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class CompressiblePorosityLeijnseEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit CompressiblePorosityLeijnseEvaluator(Teuchos::ParameterList& plist);
  CompressiblePorosityLeijnseEvaluator(const CompressiblePorosityLeijnseEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<CompressiblePorosityLeijnseModelPartition> get_Models() { return models_; }

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key poro_key_;
  Key pres_key_;

  Teuchos::RCP<CompressiblePorosityLeijnseModelPartition> models_;

 private:
  static Utils::RegisteredFactory<Evaluator, CompressiblePorosityLeijnseEvaluator> fac_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
