/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Joe Beisman
           Fenming Yuan
           Ethan Coon
*/
//! Wrapper for driving ATS from ELM.

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

// call driver advance(dt)
void ats_advance_c(ELM_ATSDriver_ptr ats,
                   double const * const dt,
                   bool const * const checkpoint,
                   bool const * const visualize);

// call driver advance_test()
void ats_advance_test_c(ELM_ATSDriver_ptr ats);

//
// Called prior to run
// -----------------------------------------------------------------------------

// call driver get_mesh_info()
// ncols_local, ncols_global, and ncells_per_col are scalars
// dz & depth are 1D arrays array of length (ncells) - these could likely only be ncells_per_col long
// elev, surf_area_m2, lat, lon are 1D arrays of length ncols - lat lon necessary for every cell?
void ats_get_mesh_info_c(ELM_ATSDriver_ptr ats,
                         int * const ncols_local,
                         int * const ncols_global,
                         double * const lat,
                         double * const lon,
                         double * const elev,
                         double * const surf_area,
                         int * const pft,
                         int * const nlevgrnd,
                         double * const depth);

// call driver setup()
void ats_setup_c(ELM_ATSDriver_ptr ats);

// call driver initialize()
void ats_initialize_c(ELM_ATSDriver_ptr ats,
                      double const * const t,
                      double const * const patm,
                      double const * const soilp);

// set material parameters, which are constant in time
void ats_set_soil_hydrologic_parameters_c(ELM_ATSDriver_ptr ats,
        double const * const base_porosity,
        double const * const hydraulic_conductivity,
        double const * const clapp_horn_b,
        double const * const clapp_horn_smpsat,
        double const * const clapp_horn_sr);

// set veg parameters, which are constant in time
void ats_set_veg_parameters_c(ELM_ATSDriver_ptr ats,
        double const * const mafic_potential_full_turgor,
        double const * const mafic_potential_wilt_point);

//
// Called prior to timestep advance
// -----------------------------------------------------------------------------
// set hydrologic properties, non-constant in time
void ats_set_soil_hydrologic_properties_c(ELM_ATSDriver_ptr ats,
        double const * const effective_porosity);

// set veg properties, non-constant in time
void ats_set_veg_properties_c(ELM_ATSDriver_ptr ats,
        double const * const rooting_fraction);

// call driver set_sources()
// soil_infiltration & soil_evaporation are 1D arrays of length ncols
// root_transpiration is a 1D array array of length (ncells)
void ats_set_potential_sources_c(ELM_ATSDriver_ptr ats,
        double const * const surface_infiltration,
        double const * const surface_evaporation,
        double const * const subsurface_transpiration);


//
// Called after timestep advance
// -----------------------------------------------------------------------------
// Get the water solution after the step is complete
// surface_pressure is a 1D array of length ncols
// soil_pressure & saturation are 1D arrays array of length (ncells)
void ats_get_waterstate_c(ELM_ATSDriver_ptr ats,
                          double * const surface_ponded_depth,
                          double * const soil_pressure,
                          double * const soil_psi,
                          double * const sat_liq,
                          double * const sat_ice);

void ats_get_water_fluxes_c(ELM_ATSDriver_ptr ats,
                            double * const soil_infiltration,
                            double * const evaporation,
                            double * const transpiration,
                            double * net_subsurface_fluxes,
                            double * net_runon);

#ifdef __cplusplus
}
#endif

// include guard
#endif
