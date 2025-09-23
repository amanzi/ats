/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Pool transfer, transfer of carbon between pools
  Koven et al 13, eqn 1, sum_{i /= j} (1-r_i)T_ij k_j C_j

*/

#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "pool_transfer_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace BGC {
namespace BGCRelations {

PoolTransferEvaluator::PoolTransferEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  // dependencies
  carbon_key_ = plist_.get<std::string>("SOM key", "soil_organic_matter");
  dependencies_.insert(carbon_key_);
  decay_key_ = plist_.get<std::string>("pool decay rate key", "soil_carbon_decay_rate");
  dependencies_.insert(decay_key_);

  my_keys_.emplace_back((plist_.get<std::string>("soil carbon transfer key",
          "soil_carbon_transfer_rate"));
  my_keys_.emplace_back((plist_.get<std::string>("soil co2 production key",
          "soil_co2_production_rate"));
}

// partition key
partition_key_ = plist_.get<std::string>("partition key", "computational_domain");
init_model_ = false;
}


PoolTransferEvaluator::PoolTransferEvaluator(const PoolTransferEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other),
    carbon_key_(other.carbon_key_),
    decay_key_(other.decay_key_),
    partition_key_(other.partition_key_),
    resp_frac_(other.resp_frac_),
    transfer_frac_(other.transfer_frac_),
    init_model_(other.init_model_)
{}

Teuchos::RCP<Evaluator>
PoolTransferEvaluator::Clone() const
{
  return Teuchos::rcp(new PoolTransferEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
PoolTransferEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                                      const std::vector<Teuchos::Ptr<CompositeVector>>& results)
{
  if (!init_model_) InitModel_(S, result->NumVectors("cell"));

  Teuchos::RCP<const CompositeVector> carbon_cv = S->GetPtr<CompositeVector>(carbon_key_);
  const AmanziMesh::Mesh& mesh = *carbon_cv->Mesh();

  const Epetra_MultiVector& C = *carbon_cv->ViewComponent("cell", false);
  const Epetra_MultiVector& k =
    *S->GetPtr<CompositeVector>(decay_key_)->ViewComponent("cell", false);
  Epetra_MultiVector& transfer_c = *results[0]->ViewComponent("cell", false);
  Epetra_MultiVector& co2_c = *results[1]->ViewComponent("cell", false);
  transfer_c.PutScalar(0.);
  co2_c.PutScalar(0.);

  const MeshPartition& part = *S->GetMeshPartition(partition_key_);
  for (int c = 0; c != res_c.MyLength(); ++c) {
    const Epetra_SerialDenseMatrix& Tij = transfer_frac_[part[c]];
    const Epetra_SerialDenseVector& ri = resp_frac_[part[c]];

    for (int p = 0; p != res_c.NumVectors(); ++p) {
      double turnover = k[p][c] * C[p][c];

      // pool loss due to turnover
      transfer_c[p][c] -= turnover;

      for (int n = 0; n != res_c.NumVectors(); ++n) {
        double transfer = turnover * Tij[p][n];
        transfer_c[n][c] += transfer * (1 - ri[p]);
        co2_c[n][c] += transfer * ri[p];
      }
    }
  }
}


void
PoolTransferEvaluator::EvaluateFieldPartialDerivative_(
  const Teuchos::Ptr<State>& S,
  Key wrt_key,
  const std::vector<Teuchos::Ptr<CompositeVector>>& results)
{
  AMANZI_ASSERT(0);
}

void
PoolTransferEvaluator::InitModel_(const Teuchos::Ptr<State>& S, int npools)
{
  Teuchos::RCP<const MeshPartition> part = S->GetMeshPartition(partition_key_);
  Teuchos::ParameterList& models_list = plist_.sublist("models");

  const std::vector<std::string>& regions = part->regions();
  for (std::vector<std::string>::const_iterator r = regions.begin(); r != regions.end(); ++r) {
    // currently only handles century!
    Teuchos::ParameterList& model_list = models_list.sublist(*r);
    AMANZI_ASSERT(model_list.get<std::string>("model type", "century") == "century");
    AMANZI_ASSERT(model_list.get<int>("number of pools", 7) <= npools);
    AMANZI_ASSERT(npools == 7);
    double percent_sand = model_list.get<double>("percent sand");
    InitCentryModel_(percent_sand);
  }
}


void
PoolTransferEvaluator::InitCenturyModel_(double percent_sand)
{
  double tt = 0.85 - 0.68 * 0.01 * (100 - percent_sand);

  Epetra_SerialDenseVector RespF(7);
  // initialize the respiration fraction
  RespF[0] = 0.0;
  RespF[1] = 0.55;
  RespF[2] = 0.5;
  RespF[3] = 0.5;
  RespF[4] = tt;
  RespF[5] = 0.55;
  RespF[6] = 0.55;
  resp_frac_.push_back(RespF);

  Epetra_SerialDenseMatrix Tij(7, 7);
  // initialize conversion factors
  Tij[0][2] = 0.76;
  Tij[0][3] = 0.24;
  Tij[1][4] = 1.0;
  Tij[2][4] = 1.0;
  Tij[3][5] = 1.0;
  Tij[4][5] = 1.0 - 0.004 / (1.0 - tt);
  Tij[4][6] = 1.0 - Tij[4][5];
  Tij[5][4] = 0.93;
  Tij[5][6] = 0.07;
  Tij[6][4] = 1.0;
  transfer_frac_.push_back(Tij);
}

} //namespace
} //namespace
} //namespace
