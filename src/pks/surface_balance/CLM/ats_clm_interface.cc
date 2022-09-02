/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * C-to-F90 wrapper for CLM
 *----------------------------------------------------------------------------*/


#include "ats_clm_interface.hh"
#include "ats_clm_interface_private.hh"

namespace ATS {
namespace CLM {


//
// Begin initialization, allocating space for driver, grid.
// ------------------------------------------------------------------
// Input:
//   ncells     | number of grid cells in the subsurface
//   ncolumns   | number of columns.  Note it is currently assumed that
//              |  ncolumns divides ncells evenly, with a fixed number of
//              |  cells per column.  (FIXME)
//   startcode  | 1 if new run, 2 if restart (not currently supported?)
//   rank       | rank (for logging purposes only)
//   verbosity  | 0 (None) - 1 (Low) - 2 (High) - 3 (Extreme)
//   
int init(int ncells, int ncolumns, int startcode, int rank, int verbosity) {
  int col_inds[ncolumns][2];
  int count = 0;
  int ncells_per = ncells / ncolumns;
  for (int i=0; i!=ncolumns; ++i) {
    col_inds[i][0] = count+1;
    count += ncells_per;
    col_inds[i][1] = count;
  }
  ats_clm_init(&ncells, &ncolumns, &col_inds[0][0], &startcode, &rank, &verbosity);
  return 0;
}


//
// Set the reference, time 0, in years.
// ------------------------------------------------------------------
// Input:
//   zero_time  | The time from which all time is measured. [year]
//              |  i.e. 2010.0 for midnight Jan 1, 2010.  Used in
//              |  determining sun/shade factors for canopy.
//   
int set_zero_time(double time) {
  // expects zero time in years
  ats_clm_zero_time(&time);
  return 0;
}


//
// Sets the initial, assumed uniform, state.
// ------------------------------------------------------------------
// Input:
//   temperature| Soil, snow, and water temperature [K]
//   snow depth | Initial snow depth [m]
//   
int set_initial_state(double temperature, double snow_depth) {
  ats_clm_initial_state(&temperature, &snow_depth);
  return 0;
}


//
// Set soil properties
// ------------------------------------------------------------------
// Input:
//   latlon     | lat/lon.  +Northern and Eastern hemispheres. [degrees]
//              |  Array of shape [ncolumns, 2]
//   sand,clay  | soil texture, fractions (must range from 0-1) [-] Size ncells.
//   color_index| index into soil color models?  See clm1d_varcon albsat and albdry
//              |  Size ncolumns.
//   f_ground   | Fraction of each land type.  Land types are set in drv_vegp.dat
//              |  Array of shape [ncolumns, NUM_LC_CLASSES]
//              |  Note only the dominant land type is used, as maxt is
//              |  hard-coded 1 and this assumption is used in a few places in drv.
//   
int set_ground_properties(double* latlon,
                          const Epetra_MultiVector& sand, const Epetra_MultiVector& clay,
                          const std::vector<int>& color_index,
                          double* fractional_ground) {
  ats_to_clm_ground_properties(latlon, sand[0], clay[0],
          color_index.data(), fractional_ground);
  return 0;
}


//
// sets the cell thicknesses in the vertical
// ------------------------------------------------------------------
// Input:
//   dz         | vector of cell dz [m]
//
int set_dz(const Epetra_MultiVector& dz) {
  ats_to_clm_dz(dz[0]);
  return 0;
}


int set_wc(const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation) {
  ats_to_clm_wc(porosity[0], saturation[0]);
  return 0;
}

int set_tksat_from_porosity(const Epetra_MultiVector& porosity) {
  ats_to_clm_tksat_from_porosity(porosity[0]);
  return 0;
}

int set_pressure(const Epetra_MultiVector& pressure, double patm) {
  ats_to_clm_pressure(pressure[0], &patm);
  return 0;
}

int set_et_controls(int beta_type, int veg_water_stress_type,
                    double wilting_point, double field_capacity,
                    double res_sat) {
  ats_to_clm_et_controls(&beta_type, &veg_water_stress_type,
                         &wilting_point, &field_capacity, &res_sat);
  return 0;
}

int set_met_data(const Epetra_MultiVector& qSW, const Epetra_MultiVector& qLW,
                 const Epetra_MultiVector& pRain, const Epetra_MultiVector& pSnow,
                 const Epetra_MultiVector& air_temp, const Epetra_MultiVector& vp_air,
                 const Epetra_MultiVector& wind_u, double patm) {
  // MOVE the unit conversions here to ats_clm.F90 to be consistent with everything else FIXME
  Epetra_MultiVector precip(pRain);
  precip.Update(1000.,pSnow,1000.); // converts m/s --> mm/s 
  Epetra_MultiVector wind_y(wind_u);
  wind_y.PutScalar(0.);
  Epetra_MultiVector patm_v(vp_air);
  patm_v.PutScalar(patm);
  
  ats_to_clm_met_data(qSW[0], qLW[0], precip[0], air_temp[0], vp_air[0], wind_u[0], wind_y[0], patm_v[0]);
  return 0;
}



int setup_begin() {
  ats_clm_setup_begin();
  return 0;
}

int setup_end() {
  ats_clm_setup_end();
  return 0;
}

int advance_time(int step, double time, double dt) {
  ats_clm_advance_time(&step, &time, &dt);
  return 0;
}




// get output
int get_total_energy_fluxes(Epetra_MultiVector& latent_heat, Epetra_MultiVector& sensible_heat,
                            Epetra_MultiVector& lw_out, Epetra_MultiVector& conducted_e) {
  clm_to_ats_total_energy_fluxes(latent_heat[0], sensible_heat[0], lw_out[0], conducted_e[0]);
  return 0;
}
int get_ground_energy_fluxes(Epetra_MultiVector& latent_heat, Epetra_MultiVector& sensible_heat,
                            Epetra_MultiVector& lw_out, Epetra_MultiVector& conducted_e) {
  clm_to_ats_ground_energy_fluxes(latent_heat[0], sensible_heat[0], lw_out[0], conducted_e[0]);
  return 0;
}

int get_mass_fluxes(Epetra_MultiVector& evap_total, Epetra_MultiVector& evap_ground,
                     Epetra_MultiVector& evap_soil, Epetra_MultiVector& evap_canopy,
                     Epetra_MultiVector& tran_veg, Epetra_MultiVector& influx,
                     Epetra_MultiVector& irrigation, Epetra_MultiVector& inst_irrigation,
                    Epetra_MultiVector& irrigation_flag, Epetra_MultiVector& tran_soil) {
  clm_to_ats_mass_fluxes(evap_total[0], evap_ground[0], evap_soil[0], evap_canopy[0],
                         tran_veg[0], influx[0],
                         irrigation[0], inst_irrigation[0], irrigation_flag[0],
                         tran_soil[0]);
  evap_total.Scale(1.e-3); // to m/s
  evap_ground.Scale(1.e-3); // to m/s
  evap_soil.Scale(1.e-3); // to m/s
  evap_canopy.Scale(1.e-3); // to m/s
  tran_veg.Scale(1.e-3); // to m/s
  influx.Scale(1.e-3); // to m/s
  irrigation.Scale(1.e-3); // to m/s
  inst_irrigation.Scale(1.e-3); // to 1/s ??
  tran_soil.Scale(1.e-3); // to 1/s
  return 0;
}

int get_diagnostics(Epetra_MultiVector& swe, Epetra_MultiVector& snow_depth,
                    Epetra_MultiVector& canopy_storage, Epetra_MultiVector& Tskin,
                    Epetra_MultiVector& Tveg, Epetra_MultiVector& Tsoil) {
  clm_to_ats_diagnostics(swe[0], snow_depth[0], canopy_storage[0],
                         Tskin[0], Tveg[0], Tsoil[0]);

  // swe, canopy storage in mm, convert to m
  swe.Scale(1.e-3);
  canopy_storage.Scale(1.e-3);

  // CLM seems to think snow_depth is already in meters, but i'm not sure I
  // believe this... FIXME      
  return 0;
}

int get_total_mass_fluxes(Epetra_MultiVector& qW_surf, Epetra_MultiVector& qW_subsurf) {
  clm_to_ats_total_mass_fluxes(qW_surf[0], qW_subsurf[0]);
  qW_surf.Scale(1.e-3); // convert to [m/s]
  qW_subsurf.Scale(1.e-3); // convert to [1/s]
  return 0;
}




} // namespace CLM
} // namespace ATS
