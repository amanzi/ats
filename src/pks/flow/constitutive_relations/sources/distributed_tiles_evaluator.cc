/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#include "Key.hh"
#include "distributed_tiles_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


DistributedTilesRateEvaluator::DistributedTilesRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondary(plist)
{
  // my_keys are for distributed subsurface source, accumulated catchment-wide
  // source.  Order matters, so we must figure out which we were given and
  // which we need.
  Key akey = my_keys_.front().first;
  Key domain_name = Keys::getDomain(akey);
  Tag tag = my_keys_.front().second;
  my_keys_.clear();

  // there is no real good way to guess, so we just force the user to provide
  // if they want to override the defaults.
  acc_sources_key_ = Keys::readKey(plist, domain_, "accumulated source", "accumulated_source");
  dist_sources_key_ = Keys::readKey(plist, domain_, "distributed source", "distributed_source");
  if ((akey != acc_sources_key_) && (akey != dist_sources_key_)) {
    Errors::Message msg;
    msg << "DistributedTilesRateEvaluator: key requested \"" << akey
        << "\" not recognized as either the \"accumulated source key\""
        << "or the \"distributed source key\".  Please provide this name in the input spec.";
    Exceptions::amanzi_throw(msg);
  }
  my_keys_.emplace_back(KeyTag{ acc_sources_key_, tag });
  my_keys_.emplace_back(KeyTag{ dist_sources_key_, tag });

  // dependencies
  catch_id_key_ = Keys::readKey(plist, domain_, "catchment ID", "catchments_id");
  dependencies_.insert(KeyTag{ catch_id_key_, tag });

  factor_key_ = Keys::readKey(plist, domain_, "factor field", "NOT_PROVIDED");
  if (factor_key_ == Keys::getKey(domain_, "NOT_PROVIDED")) {
    factor_key_ = "";
  } else {
    dependencies_.insert(KeyTag{ factor_key_, tag });
  }

  pres_key_ = Keys::readKey(plist, domain_, "pressure", "pressure");
  dependencies_.insert(KeyTag{ pres_key_, tag });
  mol_dens_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ mol_dens_key_, tag });

  // other parameters
  num_ditches_ = plist.get<int>("number of ditches");
  k_ = plist.get<double>("tile permeability [m^2]");
  num_components_ = plist.get<int>("number of components", 1);
  p_enter_ = plist.get<double>("entering pressure [Pa]", 101325);
}

void
DistributedTilesRateEvaluator::Update_(State& S)
{
  auto tag = my_keys_.front().second;

  double dt = S.Get<double>("dt", tag);

  const auto& pres = *S.Get<CompositeVector>(pres_key_, tag).ViewComponent("cell", false);
  const auto& dens = *S.Get<CompositeVector>(mol_dens_key_, tag).ViewComponent("cell", false);
  const auto& cv =
    *S.Get<CompositeVector>(Keys::getKey(domain_, "cell_volume"), tag).ViewComponent("cell", false);
  const auto& sub_marks = *S.Get<CompositeVector>(catch_id_key_, tag).ViewComponent("cell", false);

  auto& sub_sink = *S.GetW<CompositeVector>(dist_sources_key_, tag, dist_sources_key_)
                      .ViewComponent("cell", false);
  sub_sink.PutScalar(0);

  auto& acc_src_vec = S.GetW<Teuchos::Array<double>>(acc_sources_key_, tag, acc_sources_key_);
  for (int lcv = 0; lcv != num_ditches_; ++lcv) acc_src_vec[lcv] = 0.;

  AmanziMesh::Entity_ID ncells = sub_sink.MyLength();
  double total = 0.0;
  int num_vectors = 1;

  if (!factor_key_.empty()) {
    num_vectors =
      S.GetPtr<CompositeVector>(factor_key_, tag)->ViewComponent("cell", false)->NumVectors();
    AMANZI_ASSERT(num_vectors == sub_sink.NumVectors());
    AMANZI_ASSERT(num_vectors == num_components_);
  }

  if (std::abs(dt) > 1e-13) {
    const Epetra_MultiVector* factor = nullptr;
    if (!factor_key_.empty()) {
      factor = &(*(S.GetPtr<CompositeVector>(factor_key_, tag)->ViewComponent("cell", false)));
    }

    for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
      if (sub_marks[0][c] > 0) {
        double val = std::min(p_enter_ - pres[0][c], 0.0) * dens[0][c] * k_;

        if (!factor_key_.empty()) {
          for (int i = 0; i < num_components_; ++i) {
            sub_sink[i][c] = (*factor)[i][c] * val;
            acc_src_vec[sub_marks[0][c] - 1 + i * num_ditches_] +=
              (*factor)[i][c] * val * dt * cv[0][c];
            total = total + (*factor)[i][c] * val * dt * cv[0][c];
          }
        } else {
          sub_sink[0][c] = val;
          acc_src_vec[sub_marks[0][c] - 1] += val * dt * cv[0][c];
          total = total + val * dt * cv[0][c];
        }
      }
    }
  }

  total = 0.;
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) { total += sub_sink[0][c] * cv[0][c]; }

  Teuchos::RCP<const Comm_type> comm_p = S.GetMesh(domain_)->getComm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();
  MPI_Allreduce(MPI_IN_PLACE, &acc_src_vec.front(), num_ditches_, MPI_DOUBLE, MPI_SUM, comm);
}

void
DistributedTilesRateEvaluator::EnsureCompatibility(State& S)
{
  auto tag = my_keys_.front().second;
  if (!S.HasRecord(acc_sources_key_, tag)) {
    S.Require<Teuchos::Array<double>>(num_ditches_, acc_sources_key_, tag, acc_sources_key_);
  }
  S.Require<CompositeVector, CompositeVectorSpace>(dist_sources_key_, tag, dist_sources_key_)
    .SetMesh(S.GetMesh(domain_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // Cop-out -- ensure not fully implemented for this evaluator.  FIXME --ETC
  for (const auto& dep : dependencies_) { S.RequireEvaluator(dep.first, dep.second); }
}


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
