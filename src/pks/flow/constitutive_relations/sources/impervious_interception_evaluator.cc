/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "impervious_interception_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


ImperviousInterceptionEvaluator::ImperviousInterceptionEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  auto domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // dependencies
  imp_frac_key_ = Keys::readKey(plist, domain, "impervious fraction", "impervious_fraction");
  dependencies_.insert(KeyTag{ imp_frac_key_, tag });

  imp_rec_id_key_ =
    Keys::readKey(plist, domain, "impervious runoff receiver", "impervious_runoff_receiver");
  dependencies_.insert(KeyTag{ imp_rec_id_key_, tag });

  src_key_ = Keys::readKey(plist, domain, "incoming water source", "precipitation_rain");
  dependencies_.insert(KeyTag{ src_key_, tag });

  cv_key_ = Keys::readKey(plist, domain, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  Qs_max_ = plist.get<double>("maximum specific diversion rate [m s^-1]", -1);
}


void
ImperviousInterceptionEvaluator::Evaluate_(const State& S,
                                           const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& src = *S.Get<CompositeVector>(src_key_, tag).ViewComponent("cell", false);
  const auto& imp_frac = *S.Get<CompositeVector>(imp_frac_key_, tag).ViewComponent("cell", false);
  const auto& imp_rec_id =
    *S.Get<CompositeVector>(imp_rec_id_key_, tag).ViewComponent("cell", false);

  Epetra_BlockMap natural_map(src.GlobalLength(), src.MyLength(), 1, 0, src.Comm());
  if (importer_ == Teuchos::null) {
    // target map is the destination/receiving ID.  Note if there is no
    // receiving ID, we send the water to ourselves.
    std::vector<AmanziMesh::Entity_ID> target_id(src.MyLength(), -1);
    for (int c = 0; c != src.MyLength(); ++c) {
      AmanziMesh::Entity_GID gid = natural_map.GID(c);
      AmanziMesh::Entity_GID target = (AmanziMesh::Entity_GID)(imp_rec_id[0][c]);
      target_id[c] = (target < 0) ? gid : target;
    }

    Epetra_BlockMap target_map(
      src.GlobalLength(), src.MyLength(), target_id.data(), 1, 0, src.Comm());
    importer_ = Teuchos::rcp(new Epetra_Import(target_map, natural_map));
  }

  Epetra_MultiVector& modified_src = *result[0]->ViewComponent("cell", false);
  Epetra_MultiVector diverted_water(natural_map, 1);

  int source_c = natural_map.LID(897);
  int dest_c = natural_map.LID(1466);

  // first split incoming water into diverted and not diverted
  for (int c = 0; c != diverted_water.MyLength(); ++c) {
    if (imp_rec_id[0][c] < 0) {
      diverted_water[0][c] = 0.;
      modified_src[0][c] = cv[0][c] * src[0][c];
    } else {
      double rate = Qs_max_ > 0 ? std::min(Qs_max_, src[0][c]) : src[0][c];
      diverted_water[0][c] = cv[0][c] * imp_frac[0][c] * rate;
      modified_src[0][c] = cv[0][c] * src[0][c] - diverted_water[0][c];
    }
  }

  // then sum diverted into modified
  modified_src.Export(diverted_water, *importer_, Epetra_AddLocalAlso);

  // lastly divide by the local cell volume
  modified_src.ReciprocalMultiply(1., cv, modified_src, 0.);
  diverted_water.PutScalar(0.);
  diverted_water.Update(1.0, modified_src, 1.0);
  diverted_water.Update(-1.0, src, 1.0);

  // double-check mass conservation
  double orig_mass = 0., mod_mass = 0.;
  int ierr = src.Dot(cv, &orig_mass);
  ierr |= modified_src.Dot(cv, &mod_mass);
  AMANZI_ASSERT(std::abs(orig_mass - mod_mass) < 1.e-10);
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
