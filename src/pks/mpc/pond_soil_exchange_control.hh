/* Copyright 2010-202x participating institutions. BSD 3-Clause; see COPYRIGHT. */
#ifndef ATS_POND_SOIL_EXCHANGE_CONTROL_HH_
#define ATS_POND_SOIL_EXCHANGE_CONTROL_HH_
#include <algorithm>
#include <cmath>
#include <limits>
namespace Amanzi {
// Two finite stores connected by conductance G relax on tau=C0*Cs/(G*(C0+Cs)).
// The sensible lower envelope protects against leaving the latent-heat interval
// during a nonlinear trial; it changes timestep selection, never physical fluxes.
struct PondSoilExchangeDecision {
  double tau = std::numeric_limits<double>::infinity();
  double apparent_tau = std::numeric_limits<double>::infinity();
  double contrast = 0.;
  double exchange_number = 0.;
  double timestep = std::numeric_limits<double>::infinity();
};
inline PondSoilExchangeDecision PondSoilExchangeLimit(
    double depth, double pond_apparent, double pond_sensible, double soil_capacity,
    double conductance, double min_depth, double min_contrast,
    double max_exchange_number, double reference_timestep) {
  PondSoilExchangeDecision result;
  for (double value : {depth, pond_apparent, pond_sensible, soil_capacity, conductance})
    if (!std::isfinite(value)) return result;
  if (!(depth > min_depth && pond_apparent > 0. && pond_sensible > 0. &&
        soil_capacity > 0. && conductance > 0.)) return result;
  result.contrast = soil_capacity / pond_apparent;
  const double cp = std::min(pond_apparent, pond_sensible);
  result.tau = (cp / (1. + cp / soil_capacity)) / conductance;
  result.apparent_tau = (pond_apparent / (1. + pond_apparent / soil_capacity)) / conductance;
  result.exchange_number = reference_timestep / result.tau;
  if (result.contrast < min_contrast || !(result.tau > 0.) || !std::isfinite(result.tau) ||
      result.exchange_number <= max_exchange_number) return result;
  // Dyadic caps allow a prescribed hourly boundary to be reached exactly.
  double cap = reference_timestep;
  while (cap > max_exchange_number * result.tau) cap *= 0.5;
  result.timestep = cap;
  return result;
}
} // namespace Amanzi
#endif
