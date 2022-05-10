
#include "elm_ats_driver.hh"
#include "elm_ats_api.h"

#ifdef __cplusplus
extern "C"{
#endif

// allocate, call constructor and cast ptr to opaque ELM_ATS_DRIVER
ELM_ATS_DRIVER ats_create() {
  return reinterpret_cast<ELM_ATS_DRIVER>(new ATS::ELM_ATSDriver());
}
// reinterpret as elm_ats_driver and delete (calls destructor)
void ats_delete(ELM_ATS_DRIVER ats) {
  auto ats_ptr = reinterpret_cast<ATS::ELM_ATSDriver*>(ats);
  ats_ptr->finalize();
  delete ats_ptr;
}
// call driver setup()
void ats_setup(ELM_ATS_DRIVER ats, MPI_Fint *f_comm, const char *input_filename) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->setup(f_comm, input_filename);
}
// call driver initialize()
void ats_initialize(ELM_ATS_DRIVER ats){
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->initialize();
}
// call driver advance()
void ats_advance(ELM_ATS_DRIVER ats, double *dt) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance(dt);
}
// call driver advance_test()
void ats_advance_test(ELM_ATS_DRIVER ats) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance_test();
}
// call driver set_sources()
void ats_set_sources(ELM_ATS_DRIVER ats, double *soil_infiltration, double *soil_evaporation,
  double *root_transpiration, int *ncols, int *ncells) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->set_sources(soil_infiltration, soil_evaporation, root_transpiration, ncols, ncells);
}
// call driver get_waterstate()
void ats_get_waterstate(ELM_ATS_DRIVER ats, double *surface_pressure, double *soil_pressure,
  double *saturation, int *ncols, int *ncells) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->get_waterstate(surface_pressure, soil_pressure, saturation, ncols, ncells);
}
// call driver get_mesh_info()
void ats_get_mesh_info(ELM_ATS_DRIVER ats, int *ncols_local, int *ncols_global, int *ncells_per_col,
  double *dz, double *depth, double *elev, double *surf_area_m2, double *lat, double *lon) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->get_mesh_info(ncols_local, ncols_global,ncells_per_col, dz, depth, elev, surf_area_m2, lat, lon);
}
#ifdef __cplusplus
}
#endif
