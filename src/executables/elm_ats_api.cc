
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
void ats_advance(ELM_ATS_DRIVER ats, double *dt, bool visout, bool chkout) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance(dt, visout, chkout);
}
// call driver advance_test()
void ats_advance_test(ELM_ATS_DRIVER ats) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance_test();
}
// call driver advance_elmstep()
void ats_advance_elmstep(ELM_ATS_DRIVER ats, double *dt_elm, bool visout, bool chkout) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance_elmstep(dt_elm, visout, chkout);
}
//
// call driver set_mesh()
void ats_set_mesh(ELM_ATS_DRIVER ats,
  double *surf_gridsX, double *surf_gridsY, double *surf_gridsZ, double *col_verticesZ,
  const int len_gridsX, const int len_gridsY, const int len_vertices) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->set_mesh(surf_gridsX, surf_gridsY, surf_gridsZ, col_verticesZ,
    len_gridsX, len_gridsY, len_vertices);
}
// call driver set_materials()
void ats_set_materials(ELM_ATS_DRIVER ats,
  double *porosity, double *hksat, double *CH_bsw, double *CH_smpsat, double *CH_sr,
  double *eff_porosity, double *zwt){
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->set_materials(porosity, hksat, CH_bsw, CH_smpsat, CH_sr, eff_porosity, zwt);
}
// call driver set_initialconditions()
void ats_set_initialconditions(ELM_ATS_DRIVER ats, double *start_t,
  double *patm, double *soilpressure, double *wtd, bool visout){
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->set_initialconditions(start_t, patm, soilpressure, wtd, visout);
}
// call driver set_boundaryconditions
void ats_set_boundaryconditions(ELM_ATS_DRIVER ats){
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->set_boundaryconditions();
}
// call driver set_sources()
void ats_set_sources(ELM_ATS_DRIVER ats, double *soil_infiltration, double *soil_evaporation,
  double *pft_transpiration, double *root_transpiration, double *soil_drainage, int *ncols, int *ncells) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->set_sources(soil_infiltration, soil_evaporation, root_transpiration, ncols, ncells);
}
//
// call driver get_waterstate()
void ats_get_waterstate(ELM_ATS_DRIVER ats, double *surface_pd, double *soil_pressure, double *soil_psi,
  double *sat_liq, double *sat_ice, int *ncols, int *ncells) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->get_waterstate(surface_pd, soil_pressure, soil_psi, sat_liq, sat_ice, ncols, ncells);
}
void ats_get_waterflux(ELM_ATS_DRIVER ats, double *soil_infiltration, double *soil_evaporation,
  double *root_transpiration, int *ncols, int *ncells) {
  return reinterpret_cast<ATS::ELM_ATSDriver*>(ats)
  ->get_waterflux(soil_infiltration, soil_evaporation, root_transpiration, ncols, ncells);
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
