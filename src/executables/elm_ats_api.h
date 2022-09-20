
#ifdef __cplusplus
extern "C" {
    class ELM_ATSDriver;
    typedef ELM_ATSDriver *ELM_ATS_DRIVER;
#else
    typedef struct ELM_ATS_DRIVER *ELM_ATS_DRIVER;
#endif

// allocate, call constructor and cast ptr to opaque ELM_ATS_DRIVER
ELM_ATS_DRIVER ats_create();
// reinterpret as elm_ats_driver and delete (calls destructor)
void ats_delete(ELM_ATS_DRIVER ats);
// call driver setup()
void ats_setup(ELM_ATS_DRIVER ats, MPI_Fint *f_comm, const char *input_filename);
// call driver initialize()
void ats_initialize(ELM_ATS_DRIVER ats);
// call driver advance(dt)
void ats_advance(ELM_ATS_DRIVER ats, double *dt, bool visout=false, bool chkout=false);
// call driver advance_test()
void ats_advance_test(ELM_ATS_DRIVER ats);
// call driver advance_elmstep()
void ats_advance_elmstep(ELM_ATS_DRIVER ats, double *dt_elm, bool visout=false, bool chkout=false);
//
//
void ats_set_mesh(ELM_ATS_DRIVER ats,
  double *surf_gridsX, double *surf_gridsY,
  double *surf_gridsZ, double *col_verticesZ,
  const int len_gridsX, const int len_gridsY, const int len_verticesZ);
void ats_set_materials(ELM_ATS_DRIVER ats,
  double *porosity, double* hksat, double *CH_bsw, double *CH_smpsat, double *CH_sr,
  double *eff_porosity, double *zwt);
void ats_set_initialconditions(ELM_ATS_DRIVER ats, double *start_t,
  double *patm, double *soilpressure, double *wtd, bool visout=false);
void ats_set_boundaryconditions(ELM_ATS_DRIVER ats);      // (TODO)
// call driver set_sources()
// soil_infiltration & soil_evaporation are 1D arrays of length ncols
// root_transpiration is a 1D array array of length (ncells)
void ats_set_sources(ELM_ATS_DRIVER ats, double *soil_infiltration, double *soil_evaporation,
  double *pft_transpiration, double *root_transpiration, double *soil_drainage, int *ncols, int *ncells);
//
// call driver get_waterstate()
// surface_pressure is a 1D array of length ncols
// soil_pressure & saturation are 1D arrays array of length (ncells)
void ats_get_waterstate(ELM_ATS_DRIVER ats, double *surface_pressure, double *soil_pressure, double *soil_psi,
  double *sat_liq, double *sat_ice, int *ncols, int *ncells);
void ats_get_waterflux(ELM_ATS_DRIVER ats, double *soil_infiltration, double *soil_evaporation, double *root_transpiration,
  int *ncols, int *ncells);
// call driver get_mesh_info()
// ncols_local, ncols_global, and ncells_per_col are scalars
// dz & depth are 1D arrays array of length (ncells) - these could likely only be ncells_per_col long
// elev, surf_area_m2, lat, lon are 1D arrays of length ncols - lat lon necessary for every cell?
void ats_get_mesh_info(ELM_ATS_DRIVER ats, int *ncols_local, int *ncols_global, int *ncells_per_col,
  double *dz, double *depth, double *elev, double *surf_area_m2, double *lat, double *lon);

#ifdef __cplusplus
}
#endif
