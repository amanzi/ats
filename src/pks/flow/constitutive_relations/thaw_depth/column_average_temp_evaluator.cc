/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The column average temperature evaluator gets the subsurface temperature and number of cells (related to depth), and returns the average column temperature.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "column_average_temp_evaluator.hh"

namespace Amanzi {
namespace Flow {

ColumnAverageTempEvaluator::ColumnAverageTempEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  domain_ = Keys::getDomain(my_keys_.front().first); //surface_column domain
  Key domain_ss = Keys::readDomainHint(plist_, domain_, "surface", "subsurface");

  temp_key_ = Keys::readKey(plist, domain_ss, "temperature", "temperature");
  dependencies_.insert(KeyTag{temp_key_, tag});

  depth_ = plist_.get<double>("depth from surface [m]", 0); // depth from the surface
  ncells_depth_ = plist_.get<int>("number of cells", -1); // or number of cells
}


Teuchos::RCP<Evaluator>
ColumnAverageTempEvaluator::Clone() const
{
  return Teuchos::rcp(new ColumnAverageTempEvaluator(*this));
}


void
ColumnAverageTempEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell",false);
  AMANZI_ASSERT(res_c.MyLength() == 1); // as implemented only valid for columns

  // search through the column and find the deepest unfrozen cell
  std::string domain_ss = Keys::getDomain(temp_key_);
  const auto& top_z_centroid = S.GetMesh(domain_ss)->face_centroid(0);
  AmanziGeometry::Point z_centroid(top_z_centroid);

  const auto& temp_c = *S.Get<CompositeVector>(temp_key_, tag).ViewComponent("cell", false);

  int col_cells = temp_c.MyLength();
  double temp_sum = 0;
  int count = 0 ;

  AMANZI_ASSERT (ncells_depth_ <= col_cells);
  for (int i=0; i!=col_cells; ++i) {
    if (depth_ > 0.0) {
      z_centroid = S.GetMesh(domain_ss)->face_centroid(i+1);
      double z_depth = top_z_centroid[2] - z_centroid[2];
      if (z_depth <= depth_) {
        temp_sum += temp_c[0][i];
        count += 1;
      }
      else {
        break;
      }
    } else if (ncells_depth_ > 0 && i <= ncells_depth_) {
      temp_sum += temp_c[0][i];
      count += 1;
    }
  }

  res_c[0][0] = count > 0 ? temp_sum / count : 0.0;
}

void
ColumnAverageTempEvaluator::EvaluatePartialDerivative_(const State& S,
               const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{}


// Custom EnsureCompatibility forces this to be updated once.
bool
ColumnAverageTempEvaluator::Update(State& S,
        const Key& request)
{
  bool changed = EvaluatorSecondaryMonotypeCV::Update(S, request);

  if (!updated_once_) {
    Update_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
ColumnAverageTempEvaluator::EnsureCompatibility(State& S)
{
  // note, no derivs are valid here
  EnsureCompatibility_ClaimOwnership_(S);
  EnsureCompatibility_Flags_(S);
  EnsureCompatibility_DepEvals_(S);

  CompositeVectorSpace fac;
  fac.AddComponent("cell", AmanziMesh::CELL, 1);
  EnsureCompatibility_DepsFromFac_(S, fac, true);

  EnsureCompatibility_DepEnsureCompatibility_(S);
}


} //namespace
} //namespace
