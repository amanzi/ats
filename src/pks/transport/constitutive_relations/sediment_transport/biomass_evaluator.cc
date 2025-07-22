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
  InitializeFromPlist_();
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

  nspecies_ = plist_.get<int>("number of vegetation species", 1);

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

  elev_key_ = Keys::readKey(plist_, domain_name, "elevation", "elevation");
  dependencies_.insert(KeyTag{elev_key_, tag});

  msl_key_ = Keys::readKey(plist_, domain_name, "mean sea level", "mean_sea_level");
  dependencies_.insert(KeyTag{msl_key_, tag});
}


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
  const Epetra_MultiVector& msl = *S.GetPtr<CompositeVector>(msl_key_, tag)->ViewComponent("cell", false);

  int ncells = biomass.MyLength();

  for (int n = 0; n < nspecies_; n++) {
    AMANZI_ASSERT((zmax[n] - zmin[n]) > 1e-6);
    switch (type_) {
    case 1:
      for (int c = 0; c < ncells; c++) {
        double z_b = elev[0][c] - msl[0][c];
        if ((z_b > zmin[n]) && (z_b < zmax[n])) {
          biomass[n][c] = Bmax[n] * (zmax[n] - z_b) / (zmax[n] - zmin[n]);
        } else {
          biomass[n][c] = 0.;
        }
      }
      break;
    case 2:
      for (int c = 0; c < ncells; c++) {
        double z_b = elev[0][c] - msl[0][c];
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
