/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Distributes and downregulates potential transpiration to the rooting zone.

#include "Function.hh"
#include "FunctionFactory.hh"
#include "transpiration_distribution_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
TranspirationDistributionEvaluator::TranspirationDistributionEvaluator(Teuchos::ParameterList& plist) :
  SecondaryVariableFieldEvaluator(plist)
{
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
TranspirationDistributionEvaluator::Clone() const
{
  return Teuchos::rcp(new TranspirationDistributionEvaluator(*this));
}


// Initialize by setting up dependencies
void
TranspirationDistributionEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  domain_sub_ = Keys::getDomain(my_key_);
  domain_surf_ = Keys::readDomainHint(plist_, domain_sub_, "domain", "surface");

  limiter_local_ = false;
  if (plist_.isSublist("water limiter function")) {
    Amanzi::FunctionFactory fac;
    limiter_ = Teuchos::rcp(fac.Create(plist_.sublist("water limiter function")));
  } else {
    limiter_local_ = plist_.get<bool>("water limiter local", true);
  }

  // - pull Keys from plist
  // dependency: pressure
  f_wp_key_ = Keys::readKey(plist_, domain_sub_, "plant wilting factor", "plant_wilting_factor");
  dependencies_.insert(f_wp_key_);

  // dependency: rooting_depth_fraction
  f_root_key_ = Keys::readKey(plist_, domain_sub_, "rooting depth fraction", "rooting_depth_fraction");
  dependencies_.insert(f_root_key_);

  // dependency: transpiration
  potential_trans_key_ = Keys::readKey(plist_, domain_surf_, "potential transpiration", "potential_transpiration");
  dependencies_.insert(potential_trans_key_);

  // dependency: cell volume, surface cell volume
  cv_key_ = Keys::readKey(plist_, domain_sub_, "cell volume", "cell_volume");
  dependencies_.insert(cv_key_);
  surf_cv_key_ = Keys::readKey(plist_, domain_surf_, "surface cell volume", "cell_volume");
  dependencies_.insert(surf_cv_key_);

  year_duration_ = plist_.get<double>("year duration", 1.0);
  std::string year_duration_units = plist_.get<std::string>("year duration units", "noleap");

  // deal with units
  Amanzi::Utils::Units units;
  bool flag;
  year_duration_ = units.ConvertTime(year_duration_, year_duration_units, "s", flag);
}


void
TranspirationDistributionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  // on the subsurface
  const Epetra_MultiVector& f_wp = *S->GetFieldData(f_wp_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& f_root = *S->GetFieldData(f_root_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell", false);

  // on the surface
  const Epetra_MultiVector& potential_trans = *S->GetFieldData(potential_trans_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv = *S->GetFieldData(surf_cv_key_)->ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell", false);

  double p_atm = *S->GetScalarData("atmospheric_pressure");

  auto& subsurf_mesh = *S->GetMesh(domain_sub_);
  auto& surf_mesh = *S->GetMesh(domain_surf_);

  result_v.PutScalar(0.);
  for (const auto& region_lc : land_cover_) {
    AmanziMesh::Entity_ID_List lc_ids;
    surf_mesh.get_set_entities(region_lc.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

    if (TranspirationPeriod_(S->time(), region_lc.second.leaf_on_doy, region_lc.second.leaf_off_doy)) {
      for (int sc : lc_ids) {
        double column_total = 0.;
        double f_root_total = 0.;
        double f_wp_total = 0.;
        double var_dz = 0.;
        for (auto c : subsurf_mesh.cells_of_column(sc)) {
          column_total += f_wp[0][c] * f_root[0][c] * cv[0][c];
          result_v[0][c] = f_wp[0][c] * f_root[0][c];
          if (f_wp[0][c] * f_root[0][c] > 0)
            var_dz += cv[0][c];
        }

        if (column_total > 0.) {
          double coef = potential_trans[0][sc] * surf_cv[0][sc] / column_total;
          if (limiter_.get()) {
            auto column_total_vector = std::vector<double>(1, column_total / surf_cv[0][sc]);
            double limiting_factor = (*limiter_)(column_total_vector);
            AMANZI_ASSERT(limiting_factor >= 0.);
            AMANZI_ASSERT(limiting_factor <= 1.);
            coef *= limiting_factor;
          }

          for (auto c : subsurf_mesh.cells_of_column(sc)) {
            result_v[0][c] *= coef;
            if (limiter_local_) {
              result_v[0][c] *= f_wp[0][c];
            }
          }
        }
      }
    }
  }
}


void
TranspirationDistributionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(0.); // this would be a nontrivial calculation, as it is technically nonlocal due to rescaling issues?
}


void
TranspirationDistributionEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // new state!
  if (land_cover_.size() == 0)
    land_cover_ = getLandCover(S->ICList().sublist("land cover types"),
            {"leaf_on_doy", "leaf_off_doy"});

  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(!my_key_.empty());
  auto my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  Key domain = Keys::getDomain(my_key_);
  my_fac->SetMesh(S->GetMesh(domain))
    ->SetComponent("cell", AmanziMesh::CELL, 1)
    ->SetGhosted(true);

  // Create an unowned factory to check my dependencies.
  // -- first those on the subsurface mesh
  CompositeVectorSpace dep_fac(*my_fac);
  dep_fac.SetOwned(false);
  S->RequireField(f_root_key_)->Update(dep_fac);

  CompositeVectorSpace dep_fac_one;
  dep_fac_one.SetMesh(my_fac->Mesh())
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(f_wp_key_)->Update(dep_fac_one);
  S->RequireField(cv_key_)->Update(dep_fac_one);

  // -- next those on the surface mesh
  CompositeVectorSpace surf_fac;
  surf_fac.SetMesh(S->GetMesh(Keys::getDomain(surf_cv_key_)))
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(potential_trans_key_)->Update(surf_fac);

  CompositeVectorSpace surf_fac_one;
  surf_fac_one.SetMesh(S->GetMesh(Keys::getDomain(surf_cv_key_)))
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(surf_cv_key_)->Update(surf_fac_one);

  // Recurse into the tree to propagate info to leaves.
  for (KeySet::const_iterator key=dependencies_.begin();
       key!=dependencies_.end(); ++key) {
    S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
  }
}

bool
TranspirationDistributionEvaluator::TranspirationPeriod_(double time, double leaf_on_doy, double leaf_off_doy)
{
  if (leaf_on_doy < 0 || leaf_off_doy < 0) {
    return true; // evergreen
  }

  double time_of_year = fmod(time, year_duration_);
  double leaf_on_time = leaf_on_doy * 86400;
  double leaf_off_time = leaf_off_doy * 86400;

  if (leaf_on_time < leaf_off_time) {
    // northern hemisphere
    if ((leaf_on_time <= time_of_year) && (time_of_year < leaf_off_time)) {
      //summer
      return true;
    } else {
      return false;
    }
  } else {
    // southern hemisphere
    if ((leaf_off_time <= time_of_year) && (time_of_year < leaf_on_time)) {
      // southern hemisphere summer
      return true;
    } else {
      return false;
    }
  }
}


} //namespace
} //namespace
} //namespace
