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
#include "EvaluatorBCSeepageFaceHead.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

const std::string EvaluatorBCSeepageFaceHead::eval_type = "flow BC seepage face head";

EvaluatorBCSeepageFaceHead::EvaluatorBCSeepageFaceHead(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondary(plist)
{
  Key bc_key = Keys::cleanPListName(*plist);
  Key domain_name = Keys::getDomain(bc_key);
  Tag tag(plist->get<std::string>("tag"));

  elev_key_ = Keys::readKey(*plist_, domain_name, "elevation", "elevation");
  dependencies_.insert({ elev_key_, tag });
  pd_key_ = Keys::readKey(*plist_, domain_name, "ponded depth", "ponded_depth");
  dependencies_.insert({ pd_key_, tag });
}


Teuchos::RCP<Evaluator>
EvaluatorBCSeepageFaceHead::Clone() const
{
  return Teuchos::rcp(new EvaluatorBCSeepageFaceHead(*this));
}


void
EvaluatorBCSeepageFaceHead::EnsureCompatibility(State& S)
{
  const auto& keytag = my_keys_.front();
  S.CheckIsDebugEval(keytag.first, keytag.second, "ensure compatibilitied");

  // my data -- conditional: markers are provided per face in the "_flags" patch
  S.Require<MultiPatch<double>, MultiPatchSpace>(keytag.first, keytag.second, keytag.first)
    .set_flag(Operators::OPERATOR_BC_CONDITIONAL);
  S.Require<MultiPatch<int>, MultiPatchSpace>(keytag.first + "_flags", keytag.second, keytag.first);
  EnsureCompatibility_Flags_(S);

  auto& mps =
    S.Require<MultiPatch<double>, MultiPatchSpace>(keytag.first, keytag.second, keytag.first);
  if (func_ == Teuchos::null && mps.mesh != Teuchos::null) {
    AMANZI_ASSERT(mps.size() == 0);
    func_ = Teuchos::rcp(new Functions::MeshFunction(plist_->sublist("seepage face head"),
                                                     mps.mesh,
                                                     "boundary head",
                                                     AmanziMesh::Entity_kind::FACE,
                                                     mps.flag_type));
    for (const auto& spec : *func_) mps.addPatch(std::get<1>(spec));
    AMANZI_ASSERT(mps.size() == func_->size());
    AMANZI_ASSERT(mps.size() > 0);
    AMANZI_ASSERT(mps.entity_kind == AmanziMesh::Entity_kind::FACE);

    auto& flags_mps = S.Require<MultiPatch<int>, MultiPatchSpace>(
      keytag.first + "_flags", keytag.second, keytag.first);
    flags_mps = mps;

    // dependencies
    S.Require<CompositeVector, CompositeVectorSpace>(elev_key_, keytag.second)
      .SetMesh(mps.mesh)
      ->SetGhosted(true)
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(pd_key_, keytag.second)
      .SetMesh(mps.mesh)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    EnsureCompatibility_DepEnsureCompatibility_(S);
  }
}


void
EvaluatorBCSeepageFaceHead::Update_(State& S)
{
  const auto& keytag = my_keys_.front();
  auto& mp = S.GetW<MultiPatch<double>>(keytag.first, keytag.second, keytag.first);
  auto& flags = S.GetW<MultiPatch<int>>(keytag.first + "_flags", keytag.second, keytag.first);
  S.GetRecordW(keytag.first + "_flags", keytag.second, keytag.first).set_initialized();

  double time = S.get_time(keytag.second);
  func_->Compute(time, mp);

  const auto elev_f = S.Get<CompositeVector>(elev_key_, keytag.second).viewComponent("face", true);
  const auto elev_c = S.Get<CompositeVector>(elev_key_, keytag.second).viewComponent("cell", true);
  const auto pd_c = S.Get<CompositeVector>(pd_key_, keytag.second).viewComponent("cell", true);
  const AmanziMesh::MeshCache& m = mp.space->mesh->getCache();

  for (int i = 0; i != mp.size(); ++i) {
    auto& data = mp[i].data;
    auto& flag = flags[i].data;
    const auto ids = mp[i].space->getIDs();
    Kokkos::parallel_for(
      "EvaluatorBCSeepageFaceHead::Update", ids.size(), KOKKOS_LAMBDA(const int j) {
        auto f = ids(j);
        auto cells = m.getFaceCells(f);
        auto c = cells(0);
        double hz_f = data(j, 0) + elev_f(f, 0);
        double hz_c = pd_c(c, 0) + elev_c(c, 0);
        if (hz_f >= hz_c) {
          flag(j, 0) = Operators::OPERATOR_BC_NEUMANN;
          data(j, 0) = 0.;
        } else {
          flag(j, 0) = Operators::OPERATOR_BC_DIRICHLET;
          data(j, 0) = hz_f;
        }
      });
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
