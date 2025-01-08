/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Evaluates carbon pool turnover.

*/

#ifndef AMANZI_BGCRELATIONS_POOL_DECOMP_HH_
#define AMANZI_BGCRELATIONS_POOL_DECOMP_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

class Epetra_SerialDenseVector;
class Epetra_SerialDenseMatrix;

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

class PoolTransferEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit PoolTransferEvaluator(Teuchos::ParameterList& plist);
  PoolTransferEvaluator(const PoolTransferEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const std::vector<Teuchos::Ptr<CompositeVector>>& results);
  virtual void
  EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                  Key wrt_key,
                                  const std::vector<Teuchos::Ptr<CompositeVector>>& results);


 protected:
  void InitModel_(const Teuchos::Ptr<State>& S, int npools);
  void InitCenturyModel_(double percent_sand);

 protected:
  Key carbon_key_;
  Key decay_key_;

  Key partition_key_;
  std::vector<Epetra_SerialDenseVector> resp_frac_;
  std::vector<Epetra_SerialDenseMatrix> transfer_frac_;
  bool init_model_;

 private:
  static Utils::RegisteredFactory<Evaluator, PoolTransferEvaluator> fac_;
};

} // namespace BGCRelations
} // namespace BGC
} // namespace Amanzi

#endif
