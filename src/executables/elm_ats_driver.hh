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

#pragma once

#include "Key.hh"
#include "coordinator.hh"


namespace Amanzi {
class State;
};

namespace ATS {

class ELM_ATSDriver : public Coordinator {

public:
  
  ELM_ATSDriver();
  ~ELM_ATSDriver() = default;

  virtual void setup(MPI_Fint *f_comm, const char *input_filename);
  virtual void initialize();
  virtual void advance(double *dt);
  void advance_test();
  virtual void finalize();
  void set_sources(double *soil_infiltration, double *soil_evaporation, double *root_transpiration,
    int *ncols, int *ncells);
  void get_waterstate(double *surface_pressure, double *soil_pressure, double *saturation,
    int *ncols, int *ncells);
  void get_mesh_info(int *ncols_local, int *ncols_global, int *ncells_per_col, double *dz,
    double *depth, double *elev, double *surf_area_m2, double *lat, double *lon);
private:

  void col_depth(double *dz, double *depth);

  std::unique_ptr<ELM_ATSCoordinator> elm_coordinator_;
  Teuchos::RCP<Amanzi::State> S_;
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_subsurf_;
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_surf_;

  Amanzi::Key domain_sub_;
  Amanzi::Key domain_srf_;
  Amanzi::Key sub_src_key_;
  Amanzi::Key srf_src_key_;
  Amanzi::Key pres_key_;
  Amanzi::Key pd_key_;
  Amanzi::Key satl_key_;
  Amanzi::Key por_key_;
  Amanzi::Key elev_key_;

  Amanzi::Key srf_mol_dens_key_;
  Amanzi::Key srf_mass_dens_key_;
  Amanzi::Key sub_mol_dens_key_;
  Amanzi::Key sub_mass_dens_key_;

  int ncolumns_{0};
  int ncol_cells_{0};
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




