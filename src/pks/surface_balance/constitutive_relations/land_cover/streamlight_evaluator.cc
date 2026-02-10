/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: 
*/

/*
  Drainage rate.

*/

#include "Teuchos_ParameterList.hpp"
#include "streamlight_evaluator.hh"
#include "streamlight_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

StreamlightEvaluator::StreamlightEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("streamlight model parameters");
  model_ = Teuchos::rcp(new StreamlightModel(sublist));
  InitializeFromPlist_();
}


Teuchos::RCP<Evaluator>
StreamlightEvaluator::Clone() const
{
  return Teuchos::rcp(new StreamlightEvaluator(*this));
}


// Initialize by setting up dependencies
void
StreamlightEvaluator::InitializeFromPlist_()
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  my_keys_.clear();

  qSWin_key_ = Keys::readKey(
    plist_, domain_, "incoming shortwave radiation", "incoming_shortwave_radiation");
  dependencies_.insert(KeyTag{ qSWin_key_, tag });

  lai_key_ = Keys::readKey(plist_, domain_, "leaf area index", "leaf_area_index");
  dependencies_.insert(KeyTag{ lai_key_, tag });

  ponded_depth_key_ = Keys::readKey(plist_, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ ponded_depth_key_, tag });

  lat_key_ = Keys::readKey(plist_, domain_, "latitude", "latitude");
  dependencies_.insert(KeyTag{ lat_key_, tag });

  lon_key_ = Keys::readKey(plist_, domain_, "longitude", "longitude");
  dependencies_.insert(KeyTag{ lon_key_, tag });

  start_date_ = plist_.get<std::string>("start date of forcings [MM-DD]", "01-01");
  days_offset_ = model_->DaysFromJan1(start_date_);

  Key swinc_key = Keys::readKey(plist_, domain_, "shortwave incoming gpp", "shortwave_incoming_gpp");
  my_keys_.emplace_back(KeyTag{ swinc_key, tag });
  Key stream_key = Keys::readKey(plist_, domain_, "stream gpp", "stream_gpp");
  my_keys_.emplace_back(KeyTag{ stream_key, tag });
  Key streambed_key = Keys::readKey(plist_, domain_, "streambed gpp", "streambed_gpp");
  my_keys_.emplace_back(KeyTag{ streambed_key, tag });  
}


void
StreamlightEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  auto mesh = S.GetMesh(domain_);
  Tag tag = my_keys_.front().second;

  const auto& lai = *S.Get<CompositeVector>(lai_key_, tag).ViewComponent("cell", false);
  const auto& qSWin = *S.Get<CompositeVector>(qSWin_key_, tag).ViewComponent("cell", false);
  const auto& ponded_depth = *S.Get<CompositeVector>(ponded_depth_key_, tag).ViewComponent("cell", false);
  const auto& lat = *S.Get<CompositeVector>(lat_key_, tag).ViewComponent("cell", false);
  const auto& lon = *S.Get<CompositeVector>(lon_key_, tag).ViewComponent("cell", false);

  Epetra_MultiVector& gpp_swinc = *results[0]->ViewComponent("cell", false);
  Epetra_MultiVector& gpp_stream = *results[1]->ViewComponent("cell", false);
  Epetra_MultiVector& gpp_streambed = *results[2]->ViewComponent("cell", false);

  double time_sec = S.get_time();
  std::pair<int,double> doy_hour = model_->DoyHourFromSeconds(time_sec, days_offset_);
  int doy = doy_hour.first;
  double hour = doy_hour.second;
  
  // Loop through each cell
  AmanziMesh::Entity_ID ncells = ponded_depth.MyLength();
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    // get_solar_angles
    StreamlightModel::SolarAngles solang = model_->CalcSolarAngles(doy, hour, lat[0][c], lon[0][c]);
    // corrected solar azimuth
    double solar_azimuth_shade2 = model_->CorrectSolarAzimuth(solang, doy, hour, lat[0][c], lon[0][c]);
    // get_radiative_transfer_estimates_cn_1998
    StreamlightModel::RadiativeTransfer radtrans = 
      model_->GetRadiativeTransferEstimatesCN1998(solang, doy, qSWin[0][c], lai[0][c]);
    // get_riparian_stream_shading
    std::pair<double,double> fraction_shade_bank_veg = 
      model_->GetRiparianStreamShading(solar_azimuth_shade2, ponded_depth[0][c], solang.solar_altitude);
    // get_energy_stream
    StreamlightModel::EnergyStream estrm = model_->GetEnergyStream(radtrans, solang, 
      fraction_shade_bank_veg.first, fraction_shade_bank_veg.second, qSWin[0][c], ponded_depth[0][c]);
    // consolidate_metrics, convert to GPP
    StreamlightModel::PARtoGPP gpp = model_->ConvertPARtoGPP(estrm, qSWin[0][c]);
    gpp_swinc[0][c] = gpp.GPP_sw_inc_gO2_m2d;
    gpp_stream[0][c] = gpp.GPP_stream_gO2_m2d;
    gpp_streambed[0][c] = gpp.GPP_streambed_gO2_m2d;
  }
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
