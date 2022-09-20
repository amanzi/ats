/*
ELM-ATS Driver:
Provides an interface to ATS functionality for ELM 
*/
#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "VerboseObject.hh"

#include "elm_ats_coordinator.hh"
#include "State.hh"

#include "elm_ats_data.hh"
#include "elm_ats_plist.hh"

namespace ATS {

class ELM_ATSDriver {

public:
  // default constructor and destructor
  ELM_ATSDriver() {};
  ~ELM_ATSDriver() = default;

  // methods
  void setup(MPI_Fint *f_comm, const char *input_filename);
  void initialize();
  void advance(double *dt, bool visout=false, bool chkout=false);
  void advance_test();
  void advance_elmstep(double *dt_elm, bool visout=false, bool chkout=false);
  void finalize();

  void set_mesh(double *surf_gridsX, double *surf_gridsY, double *surf_gridsZ, double *col_verticesZ,
    const int len_gridsX, const int len_gridsY, const int len_vertices);
  void set_materials(double *porosity, double *hksat, double *CH_bsw, double *CH_smpsat, double *CH_sr,
    double *eff_porosity, double *zwt);
  void set_initialconditions(double *start_t, double *patm, double *soilpressure, double *wtd, bool visout=false);
  void set_boundaryconditions();
  void set_sources(double *soil_infiltration, double *soil_evaporation,
    double *root_transpiration, int *ncols, int *ncells);

  void get_waterstate(double *surface_pressure, double *soil_pressure, double *soil_psi,
    double *saturation, double *saturation_ice, int *ncols, int *ncells);
  void get_waterflux(double *soil_infiltration, double *soil_evaporation, double *root_transpiration,
	int *ncols, int *ncells);
  void get_mesh_info(int *ncols_local, int *ncols_global, int *ncells_per_col, double *dz,
    double *depth, double *elev, double *surf_area_m2, double *lat, double *lon);

private:

  void col_depth(double *dz, double *depth);
  double HfunctionSmooth(double x, double x_1, double x_0, bool derivative=false);

  std::unique_ptr<ELM_ATSCoordinator> elm_coordinator_;
  Teuchos::RCP<Amanzi::State> S_;
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_subsurf_;
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_surf_;

  Amanzi::Key domain_sub_;
  Amanzi::Key domain_srf_;
  Amanzi::Key sub_src_key_;
  Amanzi::Key srf_src_key_;
  Amanzi::Key pres_key_;
  Amanzi::Key pc_key_;
  Amanzi::Key pd_key_;
  Amanzi::Key satl_key_;
  Amanzi::Key por_key_;
  Amanzi::Key elev_key_;
  Amanzi::Key watl_key_;

  Amanzi::Key srf_mol_dens_key_;
  Amanzi::Key srf_mass_dens_key_;
  Amanzi::Key sub_mol_dens_key_;
  Amanzi::Key sub_mass_dens_key_;

  Teuchos::RCP<Teuchos::ParameterList> pks_plist_;
  Teuchos::RCP<Teuchos::ParameterList> state_plist_;
  Amanzi::Key subpk_key_;
  Amanzi::Key srfpk_key_;
  Amanzi::Key sub_pv_key_;
  Amanzi::Key srf_pv_key_;

  Amanzi::Key srf_rain_key_, srf_pev_key_, srf_evp_key_;

  bool plist_visout_;

  int ncolumns_;
  int ncol_cells_;

  elm_data elmdata_;

};


// include here temporarily during development
// maybe place into AmanziComm.hh
// or leave as local function
#include "AmanziTypes.hh"
#ifdef TRILINOS_TPETRA_STACK

#ifdef HAVE_MPI
#include "Teuchos_MpiComm.hpp"
#else
#include "Teuchos_SerialComm.hpp"
#endif

#else // Epetra stack

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#endif // trilinos stack

inline
Amanzi::Comm_ptr_type setComm(MPI_Comm comm) {
#ifdef TRILINOS_TPETRA_STACK
#ifdef HAVE_MPI
  return Teuchos::rcp(new Teuchos::MpiComm<int>(comm));
#else
  return Teuchos::rcp(new Teuchos::SerialComm<int>());
#endif
#else
#ifdef HAVE_MPI
  return Teuchos::rcp(new Epetra_MpiComm(comm));
#else
  return Teuchos::rcp(new Epetra_SerialComm());
#endif
#endif
}

} // namespace

