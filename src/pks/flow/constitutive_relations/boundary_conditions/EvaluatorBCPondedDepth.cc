/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Key.hh"
#include "MeshFunction.hh"
#include "BCs.hh"

#include "EvaluatorBCPondedDepth.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

const std::string EvaluatorBCPondedDepth::eval_type = "flow BC ponded depth";

EvaluatorBCPondedDepth::EvaluatorBCPondedDepth(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondary(plist)
{
  Key bc_key = Keys::cleanPListName(*plist);
  Key domain_name = Keys::getDomain(bc_key);
  Tag tag(plist->get<std::string>("tag"));

  Key elev_key = Keys::readKey(*plist_, domain_name, "elevation", "elevation");
  dependencies_.insert({elev_key, tag});
}

Teuchos::RCP<Evaluator>
EvaluatorBCPondedDepth::Clone() const
{
  return Teuchos::rcp(new EvaluatorBCPondedDepth(*this));
}


void
EvaluatorBCPondedDepth::EnsureCompatibility(State& S)
{
  const auto& keytag = my_keys_.front();
  S.CheckIsDebugEval(keytag.first, keytag.second, "ensure compatibilitied");

  // my data -- require and set to type DIRICHLET
  S.Require<MultiPatch<double>, MultiPatchSpace>(keytag.first, keytag.second, keytag.first)
    .set_flag(Operators::OPERATOR_BC_DIRICHLET);

  EnsureCompatibility_Flags_(S);

  // ensure a MPS that matches the functions that were set in the input spec
  auto& mps = S.Require<MultiPatch<double>, MultiPatchSpace>(keytag.first, keytag.second, keytag.first);
  if (func_ == Teuchos::null && mps.mesh != Teuchos::null) {
    AMANZI_ASSERT(mps.size() == 0);
    func_ = Teuchos::rcp(new Functions::MeshFunction(plist_->sublist("ponded depth"),
            mps.mesh,
            "boundary ponded depth [m]",
            AmanziMesh::Entity_kind::FACE,
            mps.flag_type));
    for (const auto& spec : *func_) mps.addPatch(std::get<1>(spec));
    AMANZI_ASSERT(mps.size() == func_->size());
    AMANZI_ASSERT(mps.size() > 0);
    AMANZI_ASSERT(mps.entity_kind == AmanziMesh::Entity_kind::FACE);

    // dependencies
    const auto& dep = dependencies_.front();
    S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
      .SetMesh(mps.mesh)->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

    EnsureCompatibility_DepEnsureCompatibility_(S);
  }
}


void
EvaluatorBCPondedDepth::Update_(State& S)
{
  auto& mp = S.GetW<MultiPatch<double>>(my_keys_.front().first, my_keys_.front().second, my_keys_.front().first);

  // NOTE: EvaluatorIndependentPatchFunctions own their own data.
  double time = S.get_time(my_keys_.front().second);
  func_->Compute(time, mp);

  // add on elevation
  const KeyTag& dep = dependencies_.front();
  const auto elev_f = S.Get<CompositeVector>(dep).viewComponent("face", false);

  for (auto& p : mp) {
    auto& data = p.data;
    const auto ids = p.space->getIDs();
    Kokkos::parallel_for("EvaluatorBCPondedDepth::Update", ids.size(),
                         KOKKOS_LAMBDA(const int i) {
                           data(i,0) += elev_f(ids(i),0);
                         });
  }
}


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
