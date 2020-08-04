/*
  Distributes transpiration based upon a rooting depth and a wilting-point water-potential factor.
*/

#include "Function.hh"
#include "FunctionFactory.hh"
#include "transpiration_distribution_evaluator.hh"

namespace Amanzi {
namespace Flow {
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
  Key domain_name = Keys::getDomainPrefix(my_key_);
  Key surface_name;
  if (domain_name.empty() || domain_name == "domain") {
    surface_name = "surface";
  } else {
    surface_name = Key("surface_")+domain_name;
  }
  surface_name = plist_.get<std::string>("surface domain name", surface_name);
  if (surface_name == domain_name) {
    Errors::Message message("TranspirationDistributionEvaluator: subsurface domain name and surface domain name cannot be the same.");
    Exceptions::amanzi_throw(message);
  }

  limiter_local_ = false;
  if (plist_.isSublist("water limiter function")) {
    Amanzi::FunctionFactory fac;
    limiter_ = Teuchos::rcp(fac.Create(plist_.sublist("water limiter function")));
  } else {
    limiter_local_ = plist_.get<bool>("water limiter local", true);
  }
  
  // - pull Keys from plist
  // dependency: pressure
  f_wp_key_ = Keys::readKey(plist_, domain_name, "plant wilting factor", "plant_wilting_factor");
  dependencies_.insert(f_wp_key_);

  // dependency: rooting_depth_fraction
  f_root_key_ = Keys::readKey(plist_, domain_name, "rooting depth fraction", "rooting_depth_fraction");
  dependencies_.insert(f_root_key_);

  // dependency: transpiration
  potential_trans_key_ = Keys::readKey(plist_, surface_name, "potential transpiration", "potential_transpiration");
  dependencies_.insert(potential_trans_key_);

  // dependency: cell volume, surface cell volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(cv_key_);
  surf_cv_key_ = Keys::readKey(plist_, surface_name, "surface cell volume", "cell_volume");
  dependencies_.insert(surf_cv_key_);

  npfts_ = plist_.get<int>("number of PFTs", 1);
  if (npfts_ > 1) {
    Errors::Message message("TranspirationDistributionEvaluator: it is unclear whether this evaluator works for >1 PFT.  Contact developers if you need this.");
    Exceptions::amanzi_throw(message);
  }
    
  leaf_on_time_ = plist_.get<double>("leaf on time", 0);
  std::string leaf_on_units = plist_.get<std::string>("leaf on time units", "d");

  leaf_off_time_ = plist_.get<double>("leaf off time",366);
  std::string leaf_off_units = plist_.get<std::string>("leaf on time units", "d");

  year_duration_ = plist_.get<double>("year duration", 1.0);
  std::string year_duration_units = plist_.get<std::string>("year duration units", "noleap");

  // deal with units
  Amanzi::Utils::Units units;
  bool flag;
  leaf_on_time_ = units.ConvertTime(leaf_on_time_, leaf_on_units, "s", flag);
  leaf_off_time_ = units.ConvertTime(leaf_off_time_, leaf_off_units, "s", flag);
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

  double p_atm = *S->GetScalarData("atmospheric_pressure");

  // check for winter (FIXME -- evergreens!)
  if (!TranspirationPeriod_(S->time())) {
    result->PutScalar(0.);
    return;
  }
    
  // result, on the subsurface
  const AmanziMesh::Mesh& subsurf_mesh = *result->Mesh();
  Epetra_MultiVector& result_v = *result->ViewComponent("cell", false);

  for (int pft=0; pft!=result_v.NumVectors(); ++pft) {
    for (int sc=0; sc!=potential_trans.MyLength(); ++sc) {
      double column_total = 0.;
      double f_root_total = 0.;
      double f_wp_total = 0.;
      double var_dz = 0.;
      for (auto c : subsurf_mesh.cells_of_column(sc)) {
        column_total += f_wp[0][c] * f_root[pft][c] * cv[0][c];
        result_v[pft][c] = f_wp[0][c] * f_root[pft][c];
	if (f_wp[0][c] * f_root[pft][c] > 0)
	  var_dz += cv[0][c];
      }
      // if (column_total <= 0. && potential_trans[pft][sc] > 0.) {
      //   Errors::Message message("TranspirationDistributionEvaluator: Broken run, non-zero transpiration draw but no cells with some roots are above the wilting point.");
      //   Exceptions::amanzi_throw(message);
      // }
    
      if (column_total > 0.) {
        double coef = potential_trans[pft][sc] * surf_cv[0][sc] / column_total;
        if (limiter_.get()) {
          auto column_total_vector = std::vector<double>(1, column_total / surf_cv[0][sc]);
          double limiting_factor = (*limiter_)(column_total_vector);
          AMANZI_ASSERT(limiting_factor >= 0.);
          AMANZI_ASSERT(limiting_factor <= 1.);
          coef *= limiting_factor;
        }
        
        for (auto c : subsurf_mesh.cells_of_column(sc)) {
          result_v[pft][c] *= coef;
          if (limiter_local_) {
            result_v[pft][c] *= f_wp[0][c];
          }
	}
	
      } else {
        for (auto c : subsurf_mesh.cells_of_column(sc)) {
          result_v[pft][c] = 0.;
        }
      }

      //assert (result_v[0][0] <= potential_trans[0][0]);
// THIS ENFORCES no limiter
// #ifdef ENABLE_DBC
//       double new_col_total = 0.;
//       for (auto c : subsurf_mesh.cells_of_column(sc)) {
//         new_col_total += result_v[pft][c] * cv[0][c];
//       }
//       AMANZI_ASSERT(std::abs(new_col_total - trans_total[pft][sc]*surf_cv[0][sc]) < 1.e-8);
// #endif
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
TranspirationDistributionEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
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
      ->SetComponent("cell", AmanziMesh::CELL, npfts_)
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
      ->AddComponent("cell", AmanziMesh::CELL, npfts_);
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
TranspirationDistributionEvaluator::TranspirationPeriod_(double time) {
  double time_of_year = fmod(time, year_duration_);
  if (leaf_on_time_ <= time < leaf_off_time_) {
    return true;
  } else if (leaf_on_time_ < leaf_off_time_) {
    return false;
  } else if (leaf_off_time_ <= time < leaf_on_time_) {
    return false;
  } else {
    return true;
  }
}


} //namespace
} //namespace
} //namespace
