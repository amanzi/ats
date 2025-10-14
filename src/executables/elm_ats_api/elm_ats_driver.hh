/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Joe Beisman
*/

#pragma once

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "Key.hh"
#include "coordinator.hh"

#include "ats_variables.hh"

namespace ATS {

using namespace Amanzi;

class ELM_ATSDriver : public Coordinator {

 public:

  struct MeshInfo {
    int ncols_local;
    int ncols_global;
    int nlevgrnd;
  };

  ELM_ATSDriver(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                const Teuchos::RCP<Teuchos::Time>& wallclock_timer,
                const Teuchos::RCP<const Teuchos::Comm<int>>& teuchos_comm,
                const Amanzi::Comm_ptr_type& comm,
                int npfts = 17);

  MeshInfo getMeshInfo();
  void setup();
  void initialize();
  void advance(double dt, bool checkpoint, bool vis);
  void advanceTest();
  void finalize();

  void setScalar(const ELM::VarID& scalar_id, double in);
  double getScalar(const ELM::VarID& scalar_id);

  void setField(const ELM::VarID& var_id, double const * const in);
  void getField(const ELM::VarID& var_id, double * const out);
  double const * getFieldPtr(const ELM::VarID& var_id);
  double * getFieldPtrW(const ELM::VarID& var_id);

 private:
  template<ELM::VarID var_id>
  void setField_(double const * const in);


  template<ELM::VarID var_id>
  void getField_(double * const out);

  void copyToSurf_(double const * const in, const Key& key);
  void copyToSub_(double const * const in, const Key& key);
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

  Key gross_water_source_key_;
  Key pot_evap_key_;
  Key pot_trans_key_;

  Key pres_key_;
  Key wc_key_;
  Key surf_wc_key_;

  Key evap_key_;
  Key col_trans_key_;
  Key col_baseflow_key_;
  Key col_runoff_key_;

  Key surf_mol_dens_key_;
  Key surf_mass_dens_key_;
  Key subsurf_mol_dens_key_;
  Key subsurf_mass_dens_key_;

  Key surf_cv_key_;
  Key cv_key_;

  std::map<ELM::VarID, Key> key_map_;

};


//
// Nonmember constructor/factory reads file, converts comm to the right type.
//
ELM_ATSDriver*
createELM_ATSDriver(MPI_Fint *f_comm, const char *infile, int npfts=17);

} // namespace ATS




