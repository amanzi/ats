/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Epetra_Import.h"
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

  imp_rec_id_key_ = Keys::readKey(plist, domain, "impervious runoff receiver", "impervious_runoff_receiver");
  dependencies_.insert(KeyTag{ imp_rec_id_key_, tag });

  src_key_ = Keys::readKey(plist, domain, "incoming water source", "precipitation_rain");
  dependencies_.insert(KeyTag{ src_key_, tag });

  cv_key_ = Keys::readKey(plist, domain, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
}


void
ImperviousInterceptionEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& src = *S.Get<CompositeVector>(src_key_, tag).ViewComponent("cell", false);
  const auto& imp_frac = *S.Get<CompositeVector>(imp_frac_key_, tag).ViewComponent("cell", false);
  const auto& imp_rec_id = *S.Get<CompositeVector>(imp_rec_id_key_, tag).ViewComponent("cell", false);

  // if (importer_ == Teuchos::null) {
  //   // target map is the destination/receiving ID
  //   std::vector<AmanziMesh::Entity_ID> target_id(src.MyLength());
  //   for (int c = 0; c != src.MyLength(); ++c) {
  //     target_id[c] = (AmanziMesh::Entity_ID) imp_rec_id[0][c];
  //   }
  //   Epetra_BlockMap target_map(-1, target_id.size(), target_id.data(), 1, 0, src.Comm());
  //   importer_ = Teuchos::rcp(new Epetra_Import(target_map, src.Map()));
  // }

  Epetra_MultiVector& modified_src = *result[0]->ViewComponent("cell", false);
  Epetra_MultiVector diverted_water(modified_src.Map(), 1);
  diverted_water.PutScalar(0.);

  // first split incoming water into diverted and not diverted
  for (int c = 0; c != diverted_water.MyLength(); ++c) {
    diverted_water[0][c] = imp_rec_id[0][c] < 0 ? 0. : cv[0][c] * imp_frac[0][c] * src[0][c];
    modified_src[0][c] = imp_rec_id[0][c] < 0 ? cv[0][c] * src[0][c] : cv[0][c] * (1 - imp_frac[0][c]) * src[0][c];
  }

  // // then sum diverted into modified
  // modified_src.Import(diverted_water, *importer_, Add);

  for (int c = 0; c != diverted_water.MyLength(); ++c) {
      int rec_index = static_cast<int>(imp_rec_id[0][c]);  // Explicit cast to integer
      if (rec_index > -1) {
          modified_src[0][rec_index] += diverted_water[0][c];
      }
  }

  // lastly divide by the local cell volume
  modified_src.ReciprocalMultiply(1., cv, modified_src, 0.);

  // check if the total original and modified sources are the same
  double total_src = 0.0;
  double total_modified_src = 0.0;

  for (int c = 0; c != src.MyLength(); ++c) {
    total_src += cv[0][c] * src[0][c];
    total_modified_src += cv[0][c] * modified_src[0][c];
  }

  assert(std::abs(total_src - total_modified_src) < 1e-6 && "Total volumes differ between original and modified sources"); 
  // replace with ATS way of error message or exception

}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
