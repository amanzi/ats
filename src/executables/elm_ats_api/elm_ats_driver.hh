/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Joe Beisman
*/
//! Simulation control from ELM.

/*!

The expected order of evaluation is:

..code::

   ELM_ATSDriver ats;
   ats.setup();
   ats.get_mesh_info();
   ats.set_soil_parameters();
   ats.set_veg_parameters();
   ats.initialize();

   for (step) {
     ats.set_soil_properties();
     ats.set_veg_properties();
     ats.advance();
     ats.get_water_state();
     ats.get_water_fluxes();
   }

*/

#pragma once

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "Key.hh"
#include "coordinator.hh"

namespace ATS {

using namespace Amanzi;

class ELM_ATSDriver : public Coordinator {

public:

  ELM_ATSDriver(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                const Comm_ptr_type& comm,
                int npfts = 17);

  void finalize();
  void advance(double dt, bool checkpoint, bool vis);
  void advance_test();

  void get_mesh_info(int& ncols_local,
                     int& ncols_global,
                     double * const lat,
                     double * const lon,
                     double * const elev,
                     double * const surf_area,
                     int * const pft,
                     int& nlevgrnd,
                     double * const depth);

  void setup();

  void initialize(double t,
                  double const * const p_atm,
                  double const * const pressure);

  void set_soil_hydrologic_parameters(double const * const base_porosity,
          double const * const hydraulic_conductivity,
          double const * const clapp_horn_b,
          double const * const clapp_horn_smpsat,
          double const * const clapp_horn_sr);
  void set_veg_parameters(double const * const mafic_potential_full_turgor,
          double const * const mafic_potential_wilt_point);
  void set_soil_hydrologic_properties(double const * const effective_porosity);
  void set_veg_properties(double const * const rooting_fraction);
  void set_potential_sources(double const * const surface_infiltration,
                             double const * const surface_evaporation,
                             double const * const subsurface_transpiration);

  void get_waterstate(double * const surface_ponded_depth,
                      double * const water_table_depth,
                      double * const soil_pressure,
                      double * const soil_psi,
                      double * const sat_liq,
                      double * const sat_ice);

  void get_water_fluxes(double * const soil_infiltration,
                        double * const evaporation,
                        double * const transpiration,
                        double * const net_subsurface_fluxes,
                        double * const net_runon);

 private:
  void copyToSurf_(double const * const in, const Key& key, Key owner="");
  void copyToSub_(double const * const in, const Key& key, Key owner="");
  void copyFromSurf_(double * const out, const Key& key) const;
  void copyFromSub_(double * const out, const Key& key) const;

  void initZero_(const Key& key);

 private:
  Teuchos::RCP<Teuchos::ParameterList> elm_list_;

  int ncolumns_;
  int ncells_per_col_;
  int npfts_;

  Key domain_subsurf_;
  Key domain_surf_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_subsurf_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_surf_;

  Key lat_key_;
  Key lon_key_;
  Key elev_key_;

  Key base_poro_key_;
  Key perm_key_;
  Key ch_b_key_;
  Key ch_smpsat_key_;
  Key ch_sr_key_;

  // Key poro_key_;
  Key root_frac_key_;

  Key pot_evap_key_;
  Key pot_trans_key_;
  Key pot_infilt_key_;

  Key pd_key_;
  Key wtd_key_;
  Key pres_key_;
  Key pc_key_;
  Key sat_liq_key_;
  Key sat_ice_key_;

  Key infilt_key_;
  Key trans_key_;
  Key evap_key_;

  Key surf_mol_dens_key_;
  Key surf_mass_dens_key_;
  Key subsurf_mol_dens_key_;
  Key subsurf_mass_dens_key_;

  Key surf_cv_key_;
  Key cv_key_;
};


//
// Nonmember constructor/factory reads file, converts comm to the right type.
//
ELM_ATSDriver*
createELM_ATSDriver(MPI_Fint *f_comm, const char *infile, int npfts=17);

} // namespace ATS




