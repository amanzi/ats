/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*

*/

#include "biomass_evaluator.hh"
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {

BiomassEvaluator::BiomassEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  last_update_ = -1;
  InitializeFromPlist_();
}

BiomassEvaluator::BiomassEvaluator(const BiomassEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other),
    nspecies_(other.nspecies_),
    type_(other.type_),
    biomass_key_(other.biomass_key_),
    stem_density_key_(other.stem_density_key_),
    stem_height_key_(other.stem_height_key_),
    stem_diameter_key_(other.stem_diameter_key_),
    plant_area_key_(other.plant_area_key_),
    elev_key_(other.elev_key_),
    msl_key_(other.msl_key_)
{
  last_update_ = other.last_update_;
  update_frequency_ = other.update_frequency_;

  alpha_n.resize(nspecies_);
  alpha_h.resize(nspecies_);
  alpha_a.resize(nspecies_);
  alpha_d.resize(nspecies_);

  beta_n.resize(nspecies_);
  beta_h.resize(nspecies_);
  beta_a.resize(nspecies_);
  beta_d.resize(nspecies_);

  Bmax.resize(nspecies_);
  zmax.resize(nspecies_);
  zmin.resize(nspecies_);
  for (int i = 0; i < nspecies_; i++) {
    alpha_n[i] = other.alpha_n[i];
    alpha_h[i] = other.alpha_h[i];
    alpha_d[i] = other.alpha_d[i];
    alpha_a[i] = other.alpha_a[i];
    beta_n[i] = other.beta_n[i];
    beta_h[i] = other.beta_h[i];
    beta_d[i] = other.beta_d[i];
    beta_a[i] = other.beta_a[i];
    Bmax[i] = other.Bmax[i];
    zmax[i] = other.zmax[i];
    zmin[i] = other.zmin[i];
  }
}

Teuchos::RCP<Evaluator>
BiomassEvaluator::Clone() const
{
  return Teuchos::rcp(new BiomassEvaluator(*this));
}

void
BiomassEvaluator::InitializeFromPlist_()
{
  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  my_keys_.clear(); // clear and re-insert to ensure proper order 
  
  biomass_key_ = Keys::getKey(domain_name, "biomass");
  my_keys_.emplace_back(KeyTag{biomass_key_, tag});

  stem_density_key_ = Keys::getKey(domain_name, "stem_density");
  my_keys_.emplace_back(KeyTag{stem_density_key_, tag});

  stem_height_key_ = Keys::getKey(domain_name, "stem_height");
  my_keys_.emplace_back(KeyTag{stem_height_key_, tag});

  stem_diameter_key_ = Keys::getKey(domain_name, "stem_diameter");
  my_keys_.emplace_back(KeyTag{stem_diameter_key_, tag});

  plant_area_key_ = Keys::getKey(domain_name, "plant_area");
  my_keys_.emplace_back(KeyTag{plant_area_key_, tag});

  nspecies_ = plist_.get<int>("number of vegitation species", 1);
  //species_names_ = plist_.get<Teuchos::Array<std::string> >("species names").toVector();

  last_update_ = -1.;
  update_frequency_ = plist_.get<double>("update frequency", -1);

  type_ = plist_.get<int>("type");
  alpha_n = plist_.get<Teuchos::Array<double>>("alpha n").toVector();
  alpha_h = plist_.get<Teuchos::Array<double>>("alpha h").toVector();
  alpha_a = plist_.get<Teuchos::Array<double>>("alpha a").toVector();
  alpha_d = plist_.get<Teuchos::Array<double>>("alpha d").toVector();

  beta_n = plist_.get<Teuchos::Array<double>>("beta n").toVector();
  beta_h = plist_.get<Teuchos::Array<double>>("beta h").toVector();
  beta_a = plist_.get<Teuchos::Array<double>>("beta a").toVector();
  beta_d = plist_.get<Teuchos::Array<double>>("beta d").toVector();

  Bmax = plist_.get<Teuchos::Array<double>>("Bmax").toVector();
  zmax = plist_.get<Teuchos::Array<double>>("zmax").toVector();
  zmin = plist_.get<Teuchos::Array<double>>("zmin").toVector();

  elev_key_ = Keys::getKey(domain_name, "elevation");
  //msl_key_ = "msl";
  dependencies_.insert(KeyTag{elev_key_, tag});
}


// bool
// BiomassEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S, Key request)
// {
//   if ((update_frequency_ > 0) && (last_update_ >= 0)) {
//     double time = S->get_time();
//     if (requests_.find(request) == requests_.end()) {
//       requests_.insert(request);
//       if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
//         *vo_->os() << my_keys_[0] << " has changed, but no need to update... " << std::endl;
//       }
//       return true;
//     }
//     if (time - last_update_ < update_frequency_) return false;
//   }

//   bool chg = EvaluatorSecondaryMonotypeCV::HasFieldChanged(S, request);
//   if (chg) last_update_ = S->get_time();

//   return chg;
// }

void
BiomassEvaluator::Evaluate_(const State& S,
                            const std::vector<CompositeVector*>& results)
{
  Epetra_MultiVector& biomass = *results[0]->ViewComponent("cell");
  Epetra_MultiVector& stem_density = *results[1]->ViewComponent("cell");
  Epetra_MultiVector& stem_height = *results[2]->ViewComponent("cell");
  Epetra_MultiVector& stem_diameter = *results[3]->ViewComponent("cell");
  Epetra_MultiVector& plant_area = *results[4]->ViewComponent("cell");
  Tag tag = my_keys_.front().second;

  const Epetra_MultiVector& elev = *S.GetPtr<CompositeVector>(elev_key_, tag)->ViewComponent("cell", false);

  int ncells = biomass.MyLength();

  double MSL = S.Get<double>("msl", tag);

  for (int n = 0; n < nspecies_; n++) {
    AMANZI_ASSERT((zmax[n] - zmin[n]) > 1e-6);
    switch (type_) {
    case 1:
      for (int c = 0; c < ncells; c++) {
        double z_b = elev[0][c] - MSL;
        if ((z_b > zmin[n]) && (z_b < zmax[n])) {
          biomass[n][c] = Bmax[n] * (zmax[n] - z_b) / (zmax[n] - zmin[n]);
        } else {
          biomass[n][c] = 0.;
        }
      }
      break;
    case 2:
      for (int c = 0; c < ncells; c++) {
        double z_b = elev[0][c] - MSL;
        if (z_b >= zmax[n]) {
          biomass[n][c] = Bmax[n];
        } else if ((z_b > zmin[n]) && (z_b < zmax[n])) {
          biomass[n][c] = Bmax[n] * (z_b - zmin[n]) / (zmax[n] - zmin[n]);
        } else if (z_b <= zmin[n]) {
          biomass[n][c] = 0.;
        }
      }
      break;
    }
    for (int c = 0; c < ncells; c++) {
      stem_diameter[n][c] = alpha_d[n] * std::pow(biomass[n][c], beta_d[n]);
      stem_height[n][c] = alpha_h[n] * std::pow(biomass[n][c], beta_h[n]);
      stem_density[n][c] = alpha_n[n] * std::pow(biomass[n][c], beta_n[n]);
      plant_area[n][c] = alpha_a[n] * std::pow(biomass[n][c], beta_a[n]);
    }
  }
}


} // namespace Amanzi
