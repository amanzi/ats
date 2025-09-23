/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "EvaluatorSubgridReturn.hh"

namespace Amanzi {
namespace ATS_Physics {

EvaluatorSubgridReturn::EvaluatorSubgridReturn(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  const Tag& tag = my_keys_.front().second;

  domain_set_ = plist_.get<std::string>("subgrid domain set");
  if (Keys::isDomainSet(domain_set_)) { // strip the :*
    domain_set_ = Keys::getDomainSetName(domain_set_);
  }
  nonlocal_dependencies_ = true; // by definition!

  // first add the non-domain-set dependencies
  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert({ cv_key_, tag });

  lwc_key_ = Keys::readKey(plist, domain_, "liquid water content", "water_content");
  dependencies_.insert({ lwc_key_, tag });

  alpha_key_ = Keys::readKey(plist, domain_, "exchange coefficient", "exchange_coefficient_0");
  dependencies_.insert({ alpha_key_, tag });

  // now add the domain-set dependencies
  mf_suffix_ = plist.get<std::string>("subgrid field suffix", "mole_fraction");
}

Teuchos::RCP<Evaluator>
EvaluatorSubgridReturn::Clone() const
{
  return Teuchos::rcp(new EvaluatorSubgridReturn(*this));
}

void
EvaluatorSubgridReturn::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const Tag& tag = my_keys_.front().second;
  const Epetra_MultiVector& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& alpha =
    *S.Get<CompositeVector>(alpha_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& lwc =
    *S.Get<CompositeVector>(lwc_key_, tag).ViewComponent("cell", false);

  AMANZI_ASSERT(results.size() == 1);
  Epetra_MultiVector& res = *results[0]->ViewComponent("cell", false);

  // this should be n subdomains + lwc, cv, and alpha?
  AMANZI_ASSERT(res.MyLength() + 3 == dependencies_.size());
  auto ds = S.GetDomainSet(domain_set_);
  int c = 0;
  for (const auto& subdomain : *ds) {
    const Epetra_MultiVector& mf = *S.Get<CompositeVector>(Keys::getKey(subdomain, mf_suffix_), tag)
                                      .ViewComponent("cell", false);

    for (int k = 0; k != mf.NumVectors(); ++k) {
      double integral = 0.;
      for (int c_sg = 0; c_sg != mf.MyLength() ; ++c_sg) integral += mf[k][c_sg];

      res[k][c] = integral * alpha[k][c] * lwc[0][c] / cv[0][c] / mf.MyLength();
    }
    c++;
  }
  AMANZI_ASSERT(c == res.MyLength());
}


void
EvaluatorSubgridReturn::EnsureCompatibility_ToDeps_(State& S)
{
  auto akeytag = my_keys_.front();
  if (dependencies_.size() == 3) {
    auto ds = S.GetDomainSet(domain_set_);
    for (const auto& subdomain : *ds) {
      Key dep = Keys::getKey(subdomain, mf_suffix_);
      dependencies_.insert({ dep, akeytag.second });
    }
  }

  const auto& my_fac =
    S.Require<CompositeVector, CompositeVectorSpace>(akeytag.first, akeytag.second);

  if (my_fac.Mesh() != Teuchos::null && my_fac.HasComponent("cell")) {
    S.Require<CompositeVector, CompositeVectorSpace>(cv_key_, akeytag.second)
      .SetMesh(my_fac.Mesh())
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(lwc_key_, akeytag.second)
      .SetMesh(my_fac.Mesh())
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    S.Require<CompositeVector, CompositeVectorSpace>(alpha_key_, akeytag.second)
      .SetMesh(my_fac.Mesh())
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, my_fac.NumVectors("cell"));

    auto ds = S.GetDomainSet(domain_set_);
    for (const auto& subdomain : *ds) {
      S.Require<CompositeVector, CompositeVectorSpace>(Keys::getKey(subdomain, mf_suffix_),
                                                       akeytag.second)
        .SetMesh(S.GetMesh(subdomain))
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, my_fac.NumVectors("cell"));
    }
  }
}

} // namespace ATS_Physics
} // namespace Amanzi
