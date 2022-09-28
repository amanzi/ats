
#include "elm_ats_driver.hh"
#include "elm_ats_api.h"

#ifdef __cplusplus
extern "C"{
#endif

// allocate, call constructor and cast ptr to opaque ELM_ATSDriver_ptr
ELM_ATSDriver_ptr ats_create_c(MPI_Fint *f_comm, const char *input_filename) {
  return reinterpret_cast<ELM_ATSDriver_ptr>(ATS::createELM_ATSDriver(f_comm, input_filename));
}

// reinterpret as elm_ats_driver and delete (calls destructor)
void ats_delete_c(ELM_ATSDriver_ptr ats) {
  auto ats_ptr = reinterpret_cast<ATS::ELM_ATSDriver*>(ats);
  ats_ptr->finalize();
  delete ats_ptr;
}

// call driver setup()
void ats_setup_c(ELM_ATSDriver_ptr ats) {
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->setup();
}

// call driver initialize()
void ats_initialize_c(ELM_ATSDriver_ptr ats, double *t, double *patm, double *soilp) {
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->initialize(*t, patm, soilp);
}

// call driver advance()
void ats_advance_c(ELM_ATSDriver_ptr ats, double *dt) {
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance(*dt);
}

// call driver advance_test()
void ats_advance_test_c(ELM_ATSDriver_ptr ats) {
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance_test();
}

void ats_set_soil_hydrologic_properties(ELM_ATSDriver_ptr ats,
        double *porosity,
        double *hksat,
        double *CH_bsw,
        double *CH_smpsat,
        double *CH_sr) {
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
    ->set_soil_hydrologic_properties(porosity, hksat, CH_bsw, CH_smpsat, CH_sr);
}

// call driver set_sources()
void ats_set_potential_sources_c(ELM_ATSDriver_ptr ats,
                     double const *soil_infiltration,
                     double const *soil_evaporation,
                     double const *transpiration) {
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
    ->set_potential_sources(soil_infiltration, soil_evaporation, transpiration);
}

void ats_get_actual_sources_c(ELM_ATSDriver_ptr ats,
                     double *soil_infiltration,
                     double *soil_evaporation,
                     double *transpiration) {
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
    ->get_actual_sources(soil_infiltration, soil_evaporation, transpiration);
}

// call driver get_waterstate()
void ats_get_waterstate_c(ELM_ATSDriver_ptr ats,
                          double *ponded_depth,
                          double *soil_pressure,
                          double *soil_pot,
                          double *sat_liq,
                          double *sat_ice) {
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
    ->get_waterstate(ponded_depth, soil_pressure, soil_pot, sat_liq, sat_ice);
}

// call driver get_mesh_info()
void ats_get_mesh_info_c(ELM_ATSDriver_ptr ats,
                         int *ncols_local,
                         int *ncols_global,
                         int *ncells_per_col,
                         double *lat,
                         double *lon,
                         double *elev,
                         double *surf_area,
                         double *dz,
                         double *depth) {
  reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
    ->get_mesh_info(*ncols_local, *ncols_global, *ncells_per_col,
                    lat, lon, elev, surf_area, dz, depth);
}
#ifdef __cplusplus
}
#endif
