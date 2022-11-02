/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Joe Beisman
*/
//! Simulation control from ELM.

/*
ELM-ATS driver:
Provides ATS functionality to an external caller like ELM
contains:
setup
initialize
advance
set_sources
get_waterstate
get_mesh_info
*/

#pragma once

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "Key.hh"
#include "coordinator.hh"

namespace ATS {

class ELM_ATSDriver : public Coordinator {

public:

  ELM_ATSDriver(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                const Amanzi::Comm_ptr_type& comm);

  void setup();
  void initialize(double t, double *p_atm, double *pressure);
  void finalize();

  // Note advance is an overload here
  void advance(double dt);
  void advance_test();

  void get_mesh_info(int& ncols_local,
                     int& ncols_global,
                     int& ncells_per_col,
                     double *lat,
                     double *lon,
                     double *elev,
                     double *surf_area,
                     double *dz,
                     double *depth);

  // set material properties
  void set_soil_hydrologic_properties(double* porosity,
          double* hydraulic_conductivity,
          double* clapp_horn_b,
          double* clapp_horn_smpsat,
          double* clapp_horn_sr);

  void set_potential_sources(double const *soil_infiltration,
                             double const *soil_evaporation,
                             double const *root_transpiration);

  void get_actual_sources(double *soil_infiltration,
                          double *soil_evaporation,
                          double *root_transpiration);

  void get_waterstate(double *ponded_depth,
                      double *soil_pressure,
                      double *soil_pot,
                      double *sat_liq,
                      double *sat_ice);

 private:
  void copyToSurf_(Epetra_MultiVector& target, double const * source);
  void copyToSub_(Epetra_MultiVector& target, double const * source);
  void copyFromSurf_(double* target, const Epetra_MultiVector& source);
  void copyFromSub_(double* target, const Epetra_MultiVector& source);

 private:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_subsurf_;
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_surf_;
  Amanzi::Key domain_subsurf_;
  Amanzi::Key domain_surf_;

  Amanzi::Key infilt_key_;
  Amanzi::Key pot_evap_key_;
  Amanzi::Key pot_trans_key_;
  Amanzi::Key evap_key_;
  Amanzi::Key trans_key_;

  Amanzi::Key pd_key_;
  Amanzi::Key pres_key_;
  Amanzi::Key satl_key_;

  Amanzi::Key poro_key_;
  Amanzi::Key perm_key_;
  Amanzi::Key elev_key_;
  Amanzi::Key surf_cv_key_;

  Amanzi::Key surf_mol_dens_key_;
  Amanzi::Key surf_mass_dens_key_;
  Amanzi::Key subsurf_mol_dens_key_;
  Amanzi::Key subsurf_mass_dens_key_;

  int ncolumns_{0};
  int ncells_per_col_{0};
};


//
// Nonmember constructor/factory reads file, converts comm to the right type.
//
ELM_ATSDriver*
createELM_ATSDriver(MPI_Fint *f_comm, const char *infile);

} // namespace ATS




