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

#include <iostream>
#include <cmath>
#include "Teuchos_ParameterList.hpp"
#include "streamlight_evaluator.hh"
#include "streamlight_model.hh"

const double DEG2RAD = M_PI / 180.0;
const double RAD2DEG = 180.0 / M_PI;

// WGS84 constants
const double a = 6378137.0;          // semi-major axis
const double f = 1.0/298.257223563;  // flattening
const double e = sqrt(2*f - f*f);    // eccentricity

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

struct LatLon {
    double lat;
    double lon;
};

// Inverse LCC (ellipsoidal)
LatLon lcc_to_latlon(double x, double y) {
    // Parameters from your project
    double phi1 = 25.0 * DEG2RAD;       // 1st standard parallel
    double phi2 = 60.0 * DEG2RAD;       // 2nd standard parallel
    double phi0 = 42.5 * DEG2RAD;       // latitude of origin
    double lambda0 = -100.0 * DEG2RAD;  // central meridian
    double x0 = 0.0;                     // false easting
    double y0 = 0.0;                     // false northing

    // Compute n
    double m1 = cos(phi1) / sqrt(1 - e*e * sin(phi1)*sin(phi1));
    double m2 = cos(phi2) / sqrt(1 - e*e * sin(phi2)*sin(phi2));

    auto t = [](double phi) { return tan(M_PI/4 - phi/2) / pow((1 - e*sin(phi))/(1 + e*sin(phi)), e/2); };

    double t1 = t(phi1);
    double t2 = t(phi2);
    double t0 = t(phi0);

    double n = log(m1/m2) / log(t1/t2);
    double F = m1 / (n * pow(t1, n));
    double rho0 = a * F * pow(t0, n);

    // Compute rho and theta
    double dx = x - x0;
    double dy = rho0 - (y - y0);
    double rho = sqrt(dx*dx + dy*dy);
    double theta = atan2(dx, dy);

    // Iteratively solve for latitude
    double t_ = pow(rho/(a*F), 1/n);
    double phi = M_PI/2 - 2*atan(t_);
    for(int i=0; i<5; i++) {  // usually 3-5 iterations are enough
        phi = M_PI/2 - 2*atan(t_ * pow((1 - e*sin(phi))/(1 + e*sin(phi)), e/2));
    }

    double lambda = lambda0 + theta/n;

    LatLon result;
    result.lat = phi * RAD2DEG;
    result.lon = lambda * RAD2DEG;
    return result;
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

  start_date_ = plist_.get<std::string>("start date of forcings [MM-DD]", "01-01");
  days_offset_ = model_->DaysFromJan1(start_date_);

  Key swinc_key = Keys::readKey(plist_, domain_, "shortwave incoming gpp", "shortwave_incoming_gpp");
  my_keys_.emplace_back(KeyTag{ swinc_key, tag });
  Key stream_key = Keys::readKey(plist_, domain_, "stream gpp", "stream_gpp");
  my_keys_.emplace_back(KeyTag{ stream_key, tag });
  Key streambed_key = Keys::readKey(plist_, domain_, "streambed gpp", "streambed_gpp");
  my_keys_.emplace_back(KeyTag{ streambed_key, tag });  

  Key swinc_key = Keys::readKey(plist_, domain_, "sw incoming", "sw_incoming");
  my_keys_.emplace_back(KeyTag{ swinc_key, tag });
  Key streamsurf_key = Keys::readKey(plist_, domain_, "stream surface", "stream_surface");
  my_keys_.emplace_back(KeyTag{ streamsurf_key, tag });
  Key streambed_key = Keys::readKey(plist_, domain_, "stream bed", "stream_bed");
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
    double cx = mesh->getCellCentroid(c)[0];
    double cy = mesh->getCellCentroid(c)[1];
    LatLon ll = lcc_to_latlon(cx, cy);
    double lat = ll.lat;
    double lon = ll.lon;

    // get_solar_angles
    StreamlightModel::SolarAngles solang = model_->CalcSolarAngles(doy, hour, lat, lon);
    
    // corrected solar azimuth
    double solar_azimuth_shade2 = model_->CorrectSolarAzimuth(solang, doy, hour, lat, lon);
    
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
