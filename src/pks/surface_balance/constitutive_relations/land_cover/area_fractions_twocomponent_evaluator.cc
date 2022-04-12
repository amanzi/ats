/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "area_fractions_twocomponent_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
AreaFractionsTwoComponentEvaluator::AreaFractionsTwoComponentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  //
  // NOTE: this evaluator simplifies the situation by assuming constant
  // density.  This make it so that ice and water see the same geometry per
  // unit pressure, which isn't quite true thanks to density differences.
  // However, we hypothesize that these differences, on the surface (unlike in
  // the subsurface) really don't matter much. --etc
  min_area_ = plist_.get<double>("minimum fractional area [-]", 1.e-5);
  if (min_area_ <= 0.) {
    Errors::Message message("AreaFractionsTwoComponentEvaluator: Minimum fractional area should be > 0.");
    Exceptions::amanzi_throw(message);
  }

  // get domain names
  domain_ = Keys::getDomain(my_key_); // surface
  domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");

  // get dependencies
  snow_depth_key_ = Keys::readKey(plist_, domain_snow_, "snow depth", "depth");
  dependencies_.insert(snow_depth_key_);
}


void
AreaFractionsTwoComponentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto mesh = result->Mesh();
  auto& res = *result->ViewComponent("cell",false);
  const auto& sd = *S->GetFieldData(snow_depth_key_)->ViewComponent("cell",false);

  for (const auto& lc : land_cover_) {
    AmanziMesh::Entity_ID_List lc_ids;
    mesh->get_set_entities(lc.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

    for (auto c : lc_ids) {
      // calculate area of land
      if (sd[0][c] >= lc.second.snow_transition_depth) {
        res[1][c] = 1.;
      } else if (sd[0][c] <= 0.) {
        res[1][c] = 0.;
      } else {
        res[1][c] = sd[0][c] / lc.second.snow_transition_depth;
      }

      // if any area is less than eps, give to other
      if (res[1][c] < min_area_) {
        res[1][c] = 0.;
      } else if (res[1][c] > (1-min_area_)) {
        res[1][c] = 1.;
      }
      res[0][c] = 1 - res[1][c];
    }
  }

  // debugging for bad input files
  int nerr = 0;
  for (int c=0; c!=res.MyLength(); ++c) {
    if (std::abs(1 - res[0][c] - res[1][c]) > 1e-10) nerr++;
  }
  int nerr_global = 0;
  mesh->get_comm()->SumAll(&nerr, &nerr_global, 1);
  if (nerr_global > 0) {
    Errors::Message msg("AreaFractionsTwoComponent: land cover types do not cover the mesh.");
    Exceptions::amanzi_throw(msg);
  }
}

void
AreaFractionsTwoComponentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(0.);
  // Errors::Message msg("NotImplemented: AreaFractionsTwoComponentEvaluator currently does not provide derivatives.");
  // Exceptions::amanzi_throw(msg);
}


void
AreaFractionsTwoComponentEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  if (land_cover_.size() == 0)
    land_cover_ = getLandCover(S->ICList().sublist("land cover types"),
            {"snow_transition_depth"});

  // see if we can find a master fac
  auto my_fac = S->RequireField(my_key_, my_key_);
  my_fac->SetMesh(S->GetMesh(domain_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 2);

  // Check plist for vis or checkpointing control.
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  for (auto dep_key : dependencies_) {
    auto fac = S->RequireField(dep_key);
    if (Keys::getDomain(dep_key) == domain_snow_) {
      fac->SetMesh(S->GetMesh(domain_snow_))
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
    } else {
      fac->SetMesh(S->GetMesh(domain_))
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);
    }

    // Recurse into the tree to propagate info to leaves.
    S->RequireFieldEvaluator(dep_key)->EnsureCompatibility(S);
  }
}


} //namespace
} //namespace
} //namespace

