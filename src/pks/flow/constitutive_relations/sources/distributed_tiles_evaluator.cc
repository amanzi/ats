/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*

  Evaluates water/solute sinks which represent distributed tiles in subsurface

  License: see $ATS_DIR/COPYRIGHT
  Author: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*!

Requires the following dependencies:


*/

#include "Key.hh"
#include "distributed_tiles_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


DistributedTilesRateEvaluator::DistributedTilesRateEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondary(plist),
    compatibility_checked_(false)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  subsurface_marks_key_ = Keys::readKey(plist, domain_, "catchments_id", "catchments_id");
  dist_sources_key_ = plist.get<std::string>("accumulated source key", "subdomain_sources");
  factor_key_ = plist.get<std::string>("factor field key", "");
  num_component_ = plist.get<int>("number of components", 1);

  num_ditches_ = plist.get<int>("number of ditches");
  p_enter_ = plist.get<double>("entering pressure", 101325);
  k_ = plist.get<double>("tile permeability");
  implicit_ = plist.get<bool>("implicit drainage", true);

  dependencies_.insert(KeyTag{subsurface_marks_key_, tag});
  pres_key_ = Keys::readKey(plist, domain_, "pressure", "pressure");
  dependencies_.insert(KeyTag{pres_key_, tag});
  mol_dens_key_ = Keys::readKey(plist, domain_, "molar density", "molar_density_liquid");
  dependencies_.insert(KeyTag{mol_dens_key_, tag});

  if (factor_key_ != "") dependencies_.insert(KeyTag{factor_key_, tag});
}

void
DistributedTilesRateEvaluator::Update_(State& S)
{
  auto key_tag = my_keys_.front();
  auto tag = key_tag.second;

  double dt = S.Get<double>("dt", tag);

  const auto& pres = *S.Get<CompositeVector>(pres_key_, tag).ViewComponent("cell", false);
  const auto& dens = *S.Get<CompositeVector>(mol_dens_key_, tag).ViewComponent("cell", false);
  const auto& cv = *S.Get<CompositeVector>(Keys::getKey(domain_,"cell_volume"), tag).ViewComponent("cell",false);
  const auto& sub_marks = *S.Get<CompositeVector>(subsurface_marks_key_, tag).ViewComponent("cell", false);

  auto& sub_sink = *S.GetW<CompositeVector>(key_tag.first, tag, key_tag.first).ViewComponent("cell",false);
  sub_sink.PutScalar(0);

  // NOTE: this is a hack for now -- in the future this should be made one of
  // my keys, but it would require that this eval was not a
  // EvaluatorSecondaryMonotype --ETC
  //State* S_nc = const_cast<State*>(&S);
  auto& dist_src_vec = S.GetW<Teuchos::Array<double>>(dist_sources_key_, tag, dist_sources_key_);
  for (int lcv=0; lcv!=num_ditches_; ++lcv) dist_src_vec[lcv] = 0.;

  AmanziMesh::Entity_ID ncells = sub_sink.MyLength();
  double total = 0.0;
  int num_vectors = 1;

  if (!factor_key_.empty()) {
    num_vectors = S.GetPtr<CompositeVector>(factor_key_, tag)->ViewComponent("cell",false)->NumVectors();
    AMANZI_ASSERT(num_vectors == sub_sink.NumVectors());
    AMANZI_ASSERT(num_vectors == num_component_);
  }

  if (std::abs(dt) > 1e-13) {
    const Epetra_MultiVector *factor = nullptr;
    if (!factor_key_.empty()) {
      factor = &(*(S.GetPtr<CompositeVector>(factor_key_, tag)->ViewComponent("cell", false)));
    }

    for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
      if (sub_marks[0][c] > 0) {
        double val = std::min(p_enter_ - pres[0][c], 0.0) * dens[0][c] * k_;

        if (!factor_key_.empty()) {
          for (int i=0; i<num_component_; ++i) {
            sub_sink[i][c] = (*factor)[i][c] * val;
            dist_src_vec[sub_marks[0][c] - 1 + i* num_ditches_] += (*factor)[i][c] * val * dt * cv[0][c];
            total = total + (*factor)[i][c] * val * dt * cv[0][c];
          }
        } else {
          sub_sink[0][c] = val;
          dist_src_vec[sub_marks[0][c] - 1] +=  val * dt * cv[0][c];
          total = total + val * dt * cv[0][c];
        }
      }
    }
  }

  total = 0.;
  for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
    total += sub_sink[0][c] * cv[0][c];
  }

  Teuchos::RCP<const Comm_type> comm_p = S.GetMesh(domain_)->get_comm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();
  MPI_Allreduce(MPI_IN_PLACE, &dist_src_vec.front(), num_ditches_, MPI_DOUBLE, MPI_SUM, comm);
}

void
DistributedTilesRateEvaluator::EnsureCompatibility(State& S)
{
  Key key = my_keys_.front().first;
  auto tag = my_keys_.front().second;
  if (!S.HasRecord(dist_sources_key_, tag)) {
    S.Require<Teuchos::Array<double>>(num_ditches_, dist_sources_key_, tag, dist_sources_key_);
  }

  S.Require<CompositeVector,CompositeVectorSpace>(key, tag, key)
    .SetMesh(S.GetMesh(domain_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  
  // For dependencies, all we really care is whether there is an evaluator or
  // not.  We do not use the data at all.
  for (const auto& dep : dependencies_) {
    S.RequireEvaluator(dep.first, dep.second);
  }
  
}



} //namespace
} //namespace
} //namespace

