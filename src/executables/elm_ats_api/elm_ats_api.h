
#ifndef ELM_ATS_API_HH_
#define ELM_ATS_API_HH_

#ifdef __cplusplus
extern "C" {

// opaque pointer
// external caller only sees *ELM_ATSDriver_ptr - similar to void*, but better type safety 
// ATS resolves ELM_ATSDriver_ptr as real ELM_ATSDriver during linking
class ELM_ATSDriver;
typedef ELM_ATSDriver *ELM_ATSDriver_ptr;

#else
// calling code should not dereference the pointer to the ATS object
// pointer hidden behind typedef to discourage
typedef struct ELM_ATSDriver_ptr *ELM_ATSDriver_ptr;

#endif

// allocate, call constructor and cast ptr to opaque ELM_ATSDriver_ptr
ELM_ATSDriver_ptr ats_create_c(MPI_Fint *f_comm, const char *input_filename);

// reinterpret as elm_ats_driver and delete (calls destructor)
void ats_delete_c(ELM_ATSDriver_ptr ats);

// call driver setup()
void ats_setup_c(ELM_ATSDriver_ptr ats);

// call driver initialize()
void ats_initialize_c(ELM_ATSDriver_ptr ats, double *t, double *patm, double *soilp);

// call driver advance(dt)
void ats_advance_c(ELM_ATSDriver_ptr ats, double *dt);

// call driver advance_test()
void ats_advance_test_c(ELM_ATSDriver_ptr ats);

// set material properties
void ats_set_soil_hydrologic_properties_c(ELM_ATSDriver_ptr ats,
        double* porosity,
        double* hydraulic_conductivity,
        double* clapp_horn_b,
        double* clapp_horn_smpsat,
        double* clapp_horn_sr);

// call driver set_sources()
// soil_infiltration & soil_evaporation are 1D arrays of length ncols
// root_transpiration is a 1D array array of length (ncells)
void ats_set_potential_sources_c(ELM_ATSDriver_ptr ats,
        double const *soil_infiltration,
        double const *soil_evaporation,
        double const *root_transpiration);

void ats_get_actual_sources_c(ELM_ATSDriver_ptr ats,
        double *soil_infiltration,
        double *soil_evaporation,
        double *root_transpiration);

// call driver get_waterstate()
// surface_pressure is a 1D array of length ncols
// soil_pressure & saturation are 1D arrays array of length (ncells)
void ats_get_waterstate_c(ELM_ATSDriver_ptr ats,
                          double *surface_ponded_depth,
                          double *soil_pressure,
                          double *soil_psi,
                          double *sat_liq,
                          double *sat_ice);

// call driver get_mesh_info()
// ncols_local, ncols_global, and ncells_per_col are scalars
// dz & depth are 1D arrays array of length (ncells) - these could likely only be ncells_per_col long
// elev, surf_area_m2, lat, lon are 1D arrays of length ncols - lat lon necessary for every cell?
void ats_get_mesh_info_c(ELM_ATSDriver_ptr ats,
                         int *ncols_local,
                         int *ncols_global,
                         int *ncells_per_col,
                         double *lat,
                         double *lon,
                         double *elev,
                         double *surf_area,
                         double *dz,
                         double *depth);
#ifdef __cplusplus
}
#endif

// include guard
#endif
