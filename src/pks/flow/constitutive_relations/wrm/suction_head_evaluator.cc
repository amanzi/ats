/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  Suction head = \Psi( sat )

*/

#include "suction_head_evaluator.hh"

namespace Amanzi {
namespace Flow {

SuctionHeadEvaluator::SuctionHeadEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), min_val_(0.)
{
  std::string params_name = plist_.get<std::string>("model parameters", "WRM parameters");
  Teuchos::ParameterList& sublist = plist_.sublist(params_name);
  wrms_ = createWRMPartition(sublist);
  InitializeFromPlist_();
}

SuctionHeadEvaluator::SuctionHeadEvaluator(Teuchos::ParameterList& plist,
                                           const Teuchos::RCP<WRMPartition>& wrms)
  : EvaluatorSecondaryMonotypeCV(plist), wrms_(wrms), min_val_(0.)
{
  InitializeFromPlist_();
}


Teuchos::RCP<Evaluator>
SuctionHeadEvaluator::Clone() const
{
  return Teuchos::rcp(new SuctionHeadEvaluator(*this));
}


void
SuctionHeadEvaluator::InitializeFromPlist_()
{
  // dependencies
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- saturation liquid
  sat_key_ = Keys::readKey(plist_, domain_name, "saturation", "saturation_liquid");
  dependencies_.insert(KeyTag{ sat_key_, tag });

  // cutoff above 0?
  min_val_ = plist_.get<double>("minimum suction cutoff", 0.);
}


void
SuctionHeadEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result[0]->Mesh(), -1);
    wrms_->first->Verify();
  }

  Tag tag = my_keys_.front().second;
  // Evaluate suction.
  // -- Evaluate the model to calculate suction on cells.
  const Epetra_MultiVector& sat_c =
    *S.GetPtr<CompositeVector>(sat_key_, tag)->ViewComponent("cell", false);
  Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

  int ncells = res_c.MyLength();
  for (unsigned int c = 0; c != ncells; ++c) {
    int index = (*wrms_->first)[c];
    res_c[0][c] = wrms_->second[index]->suction_head(sat_c[0][c]);
  }
}


void
SuctionHeadEvaluator::EvaluatePartialDerivative_(const State& S,
                                                 const Key& wrt_key,
                                                 const Tag& wrt_tag,
                                                 const std::vector<CompositeVector*>& result)
{
  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result[0]->Mesh(), -1);
    wrms_->first->Verify();
  }

  Tag tag = my_keys_.front().second;
  if (wrt_key == sat_key_) {
    // d(psi) / dsl
    const Epetra_MultiVector& sat_c =
      *S.GetPtr<CompositeVector>(sat_key_, tag)->ViewComponent("cell", false);
    Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

    int ncells = res_c.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) {
      int index = (*wrms_->first)[c];
      res_c[0][c] = wrms_->second[index]->d_suction_head(sat_c[0][c]);
      // AMANZI_ASSERT(res_c[0][c] >= 0.);
    }
  }
}


} // namespace Flow
} // namespace Amanzi
