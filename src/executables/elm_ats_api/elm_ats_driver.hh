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

namespace Impl {


enum class Indexer_kind {
  SCALAR = 0,
  ELM,
  ATS
};
enum class Domain_kind {
  SURF,
  SUBSURF
};


template<Indexer_kind ikind, Domain_kind dkind, typename T> struct Indexer;

//
// Note, only const version valid
//
template<Domain_kind dkind, typename T>
struct Indexer<Indexer_kind::SCALAR, dkind, T> {
  Indexer(T val) : val_(val) {}

  template<typename D = double>
  typename std::enable_if<dkind == Domain_kind::SUBSURF, const D&>::type
  operator()(const int col, const int cic) const { return val_; }

  template<typename D = double>
  typename std::enable_if<dkind == Domain_kind::SURF, const D&>::type
  operator()(const int i) const { return val_; }

private:
  double val_;
};


template<Domain_kind dkind, typename T>
struct Indexer<Indexer_kind::ELM, dkind, T> {

  Indexer(T* vals, int ncells_per_col=-1)
    : vals_(vals),
      ncells_per_col_(ncells_per_col) {}

  // for use by subsurface
  template<typename D = T>
  typename std::enable_if<dkind == Domain_kind::SUBSURF, const D&>::type
  operator()(const int col, const int cic) const {
    return vals_[col * ncells_per_col_ + cic];
  }

  template<typename D = T>
  typename std::enable_if<(dkind == Domain_kind::SUBSURF) & !std::is_const<T>::value, D&>::type
  operator()(const int col, const int cic) {
    return vals_[col * ncells_per_col_ + cic];
  }

  // for use by surface
  template<typename D = T>
  typename std::enable_if<dkind == Domain_kind::SURF, const D&>::type
  operator()(const int i) const {
    return vals_[i];
  }

  template<typename D = T>
  typename std::enable_if<(dkind == Domain_kind::SURF) & !std::is_const<T>::value, D&>::type
  operator()(const int i) {
    return vals_[i];
  }

 private:
  T* vals_;
  int ncells_per_col_;
};


template<Domain_kind dkind, typename T>
struct Indexer<Indexer_kind::ATS, dkind, T> {

  Indexer(T* vals, const Amanzi::AmanziMesh::Mesh& mesh)
    : vals_(vals),
      mesh_(mesh) {}

  // for use by subsurface
  template<typename D = const double>
  typename std::enable_if<dkind == Domain_kind::SUBSURF, D&>::type
  operator()(const int col, const int cic) const {
    return (*vals_)[0][mesh_.cells_of_column(col)[cic]];
  }

  template<typename D = double>
  typename std::enable_if<(dkind == Domain_kind::SUBSURF) & !std::is_const<T>::value, D&>::type
  operator()(const int col, const int cic) {
    return (*vals_)[0][mesh_.cells_of_column(col)[cic]];
  }

  // for use by surface
  template<typename D = const double>
  typename std::enable_if<dkind == Domain_kind::SURF, D&>::type
  operator()(const int i) const {
    return (*vals_)[0][i];
  }

  template<typename D = double>
  typename std::enable_if<(dkind == Domain_kind::SURF) & !std::is_const<T>::value, D&>::type
  operator()(const int i) {
    return (*vals_)[0][i];
  }

 private:
  T* vals_;
  const Amanzi::AmanziMesh::Mesh& mesh_;
};


template<Domain_kind dkind>
Indexer<Indexer_kind::SCALAR, dkind, double>
create_Indexer(double val) {
  return Indexer<Indexer_kind::SCALAR, dkind, double>(val);
}

template<Domain_kind dkind>
Indexer<Indexer_kind::ELM, dkind, double>
create_Indexer(double* val) {
  return Indexer<Indexer_kind::ELM, dkind, double>(val);
}

template<Domain_kind dkind>
Indexer<Indexer_kind::ELM, dkind, const double>
create_Indexer(const double* val) {
  return Indexer<Indexer_kind::ELM, dkind, const double>(val);
}

template<Domain_kind dkind>
Indexer<Indexer_kind::ATS, dkind, Epetra_MultiVector>
create_Indexer(Epetra_MultiVector* val) {
  return Indexer<Indexer_kind::ELM, dkind, Epetra_MultiVector>(val);
}

template<Domain_kind dkind>
Indexer<Indexer_kind::ELM, dkind, const Epetra_MultiVector>
create_Indexer(const Epetra_MultiVector* val) {
  return Indexer<Indexer_kind::ELM, dkind, const Epetra_MultiVector>(val);
}


template<typename Indexer, typename cIndexer>
void copySurf(Indexer&& out, cIndexer&& in, int count) {
  for (int i=0; i!=count; ++i) {
    out(i) = in(i);
  }
}


template<typename Indexer, typename cIndexer>
void copySubsurf(Indexer&& out, cIndexer&& in, int ncols, int ncells_per_col) {
  for (int col=0; col!=ncols; ++col) {
    for (int i=0; i!=ncells_per_col; ++i) {
      out(col, i) = in(col, i);
    }
  }
}

} // namespace Impl


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
  template<typename Out, typename In>
  void copySurf_(Out& out, const In& in) {
    Impl::copySurf(create_SurfIndexer(out),
                   create_SurfIndexer(in),
                   ncolumns_);
  }

  template<typename Out, typename In>
  void copySubsurf_(Out& out, const In& in) {
    Impl::copySubsurf(create_SubsurfIndexer(out),
                      create_SubsurfIndexer(in),
                      ncolumns_,
                      ncells_per_col_);
  }


  // oh how I miss C++17 now... no partial specialization, no constexpr, no
  // return auto makes this hard
  Impl::Indexer<Impl::Indexer_kind::SCALAR, Impl::Domain_kind::SURF, double>
  create_SurfIndexer(double val) {
    return Impl::Indexer<Impl::Indexer_kind::SCALAR, Impl::Domain_kind::SURF, double>(val);
  }
  Impl::Indexer<Impl::Indexer_kind::SCALAR, Impl::Domain_kind::SUBSURF, double>
  create_SubsurfIndexer(double val) {
    return Impl::Indexer<Impl::Indexer_kind::SCALAR, Impl::Domain_kind::SUBSURF, double>(val);
  }

  Impl::Indexer<Impl::Indexer_kind::ELM, Impl::Domain_kind::SURF, double>
  create_SurfIndexer(double* val) {
    return Impl::Indexer<Impl::Indexer_kind::ELM, Impl::Domain_kind::SURF, double>(val);
  }
  Impl::Indexer<Impl::Indexer_kind::ELM, Impl::Domain_kind::SURF, const double>
  create_SurfIndexer(double const* val) {
    return Impl::Indexer<Impl::Indexer_kind::ELM, Impl::Domain_kind::SURF, const double>(val);
  }
  Impl::Indexer<Impl::Indexer_kind::ELM, Impl::Domain_kind::SUBSURF, double>
  create_SubsurfIndexer(double* val) {
    return Impl::Indexer<Impl::Indexer_kind::ELM, Impl::Domain_kind::SUBSURF, double>(val, ncells_per_col_);
  }
  Impl::Indexer<Impl::Indexer_kind::ELM, Impl::Domain_kind::SUBSURF, const double>
  create_SubsurfIndexer(double const* val) {
    return Impl::Indexer<Impl::Indexer_kind::ELM, Impl::Domain_kind::SUBSURF, const double>(val, ncells_per_col_);
  }

  Impl::Indexer<Impl::Indexer_kind::ATS, Impl::Domain_kind::SURF, Epetra_MultiVector>
  create_SurfIndexer(Epetra_MultiVector& val) {
    return Impl::Indexer<Impl::Indexer_kind::ATS, Impl::Domain_kind::SURF, Epetra_MultiVector>(&val, *mesh_surf_);
  }
  Impl::Indexer<Impl::Indexer_kind::ATS, Impl::Domain_kind::SURF, const Epetra_MultiVector>
  create_SurfIndexer(Epetra_MultiVector const& val) {
    return Impl::Indexer<Impl::Indexer_kind::ATS, Impl::Domain_kind::SURF, const Epetra_MultiVector>(&val, *mesh_surf_);
  }
  Impl::Indexer<Impl::Indexer_kind::ATS, Impl::Domain_kind::SUBSURF, Epetra_MultiVector>
  create_SubsurfIndexer(Epetra_MultiVector& val) {
    return Impl::Indexer<Impl::Indexer_kind::ATS, Impl::Domain_kind::SUBSURF, Epetra_MultiVector>(&val, *mesh_subsurf_);
  }
  Impl::Indexer<Impl::Indexer_kind::ATS, Impl::Domain_kind::SUBSURF, const Epetra_MultiVector>
  create_SubsurfIndexer(Epetra_MultiVector const& val) {
    return Impl::Indexer<Impl::Indexer_kind::ATS, Impl::Domain_kind::SUBSURF, const Epetra_MultiVector>(&val, *mesh_subsurf_);
  }


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




