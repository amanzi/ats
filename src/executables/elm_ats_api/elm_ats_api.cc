/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Joe Beisman
           Fenming Yuan
           Ethan Coon
*/
//! Wrapper for driving ATS from ELM.

#include "elm_ats_driver.hh"
#include "elm_ats_api.h"

#ifdef __cplusplus
extern "C"
{
#endif

  // allocate, call constructor and cast ptr to opaque ELM_ATSDriver_ptr
  ELM_ATSDriver_ptr ats_create(MPI_Fint* f_comm, const char* input_filename)
  {
    return reinterpret_cast<ELM_ATSDriver_ptr>(ATS::createELM_ATSDriver(f_comm, input_filename));
  }

  // reinterpret as elm_ats_driver and delete (calls destructor)
  void ats_delete(ELM_ATSDriver_ptr ats)
  {
    auto ats_ptr = reinterpret_cast<ATS::ELM_ATSDriver*>(ats);
    Kokkos::finalize();
    //ats_ptr->finalize();
    delete ats_ptr;
  }

  // call driver advance()
  void ats_advance(ELM_ATSDriver_ptr ats,
                   double const* const dt,
                   bool const* const checkpoint,
                   bool const* const visualize)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance(*dt, *checkpoint, *visualize);
  }

  // call driver advance_test()
  void ats_advance_test(ELM_ATSDriver_ptr ats)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->advance_test();
  }

  void ats_get_mesh_info(ELM_ATSDriver_ptr ats,
                         int* const ncols_local,
                         int* const ncols_global,
                         double* const lat,
                         double* const lon,
                         double* const elev,
                         double* const surf_area,
                         int* const pft,
                         int* const nlevgrnd,
                         double* const depth)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->get_mesh_info(
      *ncols_local, *ncols_global, lat, lon, elev, surf_area, pft, *nlevgrnd, depth);
  }


  // call driver setup()
  void ats_setup(ELM_ATSDriver_ptr ats)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->setup();
  }

  // call driver initialize()
  void ats_initialize(ELM_ATSDriver_ptr ats,
                      double const* const t,
                      double const* const patm,
                      double const* const soilp)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->initialize(*t, patm, soilp);
  }


  void ats_set_soil_hydrologic_parameters(ELM_ATSDriver_ptr ats,
                                          double const* const base_porosity,
                                          double const* const hydraulic_conductivity,
                                          double const* const clapp_horn_b,
                                          double const* const clapp_horn_smpsat,
                                          double const* const clapp_horn_sr)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->set_soil_hydrologic_parameters(
      base_porosity, hydraulic_conductivity, clapp_horn_b, clapp_horn_smpsat, clapp_horn_sr);
  }


  void ats_set_veg_parameters(ELM_ATSDriver_ptr ats,
                              double const* const mafic_potential_full_turgor,
                              double const* const mafic_potential_wilt_point)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->set_veg_parameters(mafic_potential_full_turgor,
                                                                   mafic_potential_wilt_point);
  }


  void ats_set_soil_hydrologic_properties(ELM_ATSDriver_ptr ats,
                                          double const* const effective_porosity)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->set_soil_hydrologic_properties(effective_porosity);
  }


  // set veg properties, non-constant in time
  void ats_set_veg_properties(ELM_ATSDriver_ptr ats, double const* const rooting_fraction)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->set_veg_properties(rooting_fraction);
  }


  // call driver set_sources()
  void ats_set_sources(ELM_ATSDriver_ptr ats,
                       double const* const surface_infiltration,
                       double const* const surface_evaporation,
                       double const* const subsurface_transpiration)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->set_potential_sources(
      surface_infiltration, surface_evaporation, subsurface_transpiration);
  }


  void ats_get_waterstate(ELM_ATSDriver_ptr ats,
                          double* const surface_ponded_depth,
                          double* const water_table_depth,
                          double* const soil_pressure,
                          double* const soil_psi,
                          double* const sat_liq,
                          double* const sat_ice)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->get_waterstate(
      surface_ponded_depth, water_table_depth, soil_pressure, soil_psi, sat_liq, sat_ice);
  }


  void ats_get_water_fluxes(ELM_ATSDriver_ptr ats,
                            double* const soil_infiltration,
                            double* const evaporation,
                            double* const transpiration,
                            double* net_subsurface_fluxes,
                            double* net_runon)
  {
    reinterpret_cast<ATS::ELM_ATSDriver*>(ats)->get_water_fluxes(
      soil_infiltration, evaporation, transpiration, net_subsurface_fluxes, net_runon);
  }

#ifdef __cplusplus
}
#endif
