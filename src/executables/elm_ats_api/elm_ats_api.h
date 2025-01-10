/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Joe Beisman
           Fenming Yuan
           Ethan Coon

This file defines a C interface to the ATS library
for use with LSMs.
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
ELM_ATSDriver_ptr ats_create(MPI_Fint *f_comm, const char *input_filename);

// reinterpret as elm_ats_driver and delete (calls destructor)
void ats_delete(ELM_ATSDriver_ptr ats);

// call driver advance(dt)
void ats_advance(ELM_ATSDriver_ptr ats,
                   double const * const dt,
                   bool const * const do_vis,
                   bool const * const do_chkp);

// call driver advance_test()
void ats_advance_test(ELM_ATSDriver_ptr ats);

//
// Called prior to run
// -----------------------------------------------------------------------------

// call driver get_mesh_info()
void ats_get_mesh_info(ELM_ATSDriver_ptr ats,
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
void ats_setup(ELM_ATSDriver_ptr ats);

// call driver initialize()
void ats_initialize(ELM_ATSDriver_ptr ats,
                      double const * const start_time,
                      double const * const soil_water_content,
                      double const * const soil_pressure);

// set material parameters, which are constant in time
void ats_set_soil_hydrologic_parameters(ELM_ATSDriver_ptr ats,
        double const * const base_porosity,
        double const * const hydraulic_conductivity,
        double const * const clapp_horn_b,
        double const * const clapp_horn_smpsat,
        double const * const clapp_horn_sr);

//
// Called prior to timestep advance
// -----------------------------------------------------------------------------
// set hydrologic properties, non-constant in time
void ats_set_soil_hydrologic_properties(ELM_ATSDriver_ptr ats,
        double const * const effective_porosity);

// set veg properties, non-constant in time
void ats_set_veg_properties(ELM_ATSDriver_ptr ats,
        double const * const rooting_fraction);

// call driver set_sources()
void ats_set_sources(ELM_ATSDriver_ptr ats,
        double const * const surface_source,
        double const * const potential_evaporation,
        double const * const potential_transpiration);


//
// Called after timestep advance
// -----------------------------------------------------------------------------
// Get the water solution after the step is complete
// surface_pressure is a 1D array of length ncols
// soil_pressure & saturation are 1D arrays array of length (ncells)
void ats_get_waterstate(ELM_ATSDriver_ptr ats,
                          double * const ponded_depth,
                          double * const water_table_depth,
                          double * const mass_water_content);

void ats_get_water_fluxes(ELM_ATSDriver_ptr ats,
                            double * const infiltration,
                            double * const evaporation,
                            double * const transpiration,
                            double * const root_flux,
                            double * net_subsurface_fluxes,
                            double * net_runon);

#ifdef __cplusplus
}
#endif

// include guard
#endif
