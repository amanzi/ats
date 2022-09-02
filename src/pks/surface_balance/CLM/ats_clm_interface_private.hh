#ifndef ATS_CLM_INTERFACE_PRIVATE_HH_
#define ATS_CLM_INTERFACE_PRIVATE_HH_


extern "C" {
  void ats_to_clm_dz(const double* dz);
  void ats_to_clm_ground_properties(const double* latlon,
          const double* sand, const double* clay, const int* color_index,
          const double* fractional_ground);

  void ats_clm_init(int* ncells, int* ncolumns, int* col_inds, int* startcode,
                    int* rank, int* verbosity);
  void ats_clm_setup_begin();
  void ats_clm_setup_end();
  void ats_clm_advance_time(int* step, double* time, double* dt);

  void clm_to_ats_total_mass_fluxes(double* qW_surf, double* qW_subsurf);
  void clm_to_ats_diagnostics(double* swe, double* snow_depth, double* canopy_storage,
          double* Tskin, double* Tveg, double* Tsoil);
  void clm_to_ats_mass_fluxes(double* evap_total, double* evap_ground,
          double* evap_soil, double* evap_canopy, double* tran_veg,
          double* influx, double* irrigation, double* inst_irrigation, double* irrigation_flag,
          double* tran_soil);
  void clm_to_ats_total_energy_fluxes(double* latent_heat, double* sensible_heat,
          double* lw_out, double* conducted_e);
  void clm_to_ats_ground_energy_fluxes(double* latent_heat, double* sensible_heat,
          double* lw_out, double* conducted_e);

  void ats_to_clm_et_controls(int* beta_type, int* veg_water_stress_type,
          double* wilting_point, double* field_capacity,
          double* res_sat);
  void ats_to_clm_wc(double* porosity, double* sat);
  void ats_to_clm_pressure(double* pressure, double* patm);
  void ats_to_clm_tksat_from_porosity(double* porosity);
  void ats_to_clm_met_data(double* qSW, double* qLW, double* precip,
                           double* air_temp, double* vp_air,
                           double* wind_x, double* wind_y, double* patm);

  void ats_clm_zero_time(double* zero_year);
  void ats_clm_initial_state(double* temperature, double* swe);

}


#endif

