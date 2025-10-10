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

  void setScalar(const ScalarID& scalar_id, double in);
  double getScalar(const ScalarID& scalar_id);

  void setField(const VarID& var_id, double const * const in);
  void getField(const VarID& var_id, double * const out);
  double const * getFieldPtr(const VarID& var_id);
  double * getFieldPtrW(const VarID& var_id);

 private:
  template<VarID var_id>
  void setField_(double const * const in);


  template<VarID var_id>
  void getField_(double * const out);

  void copyToSurf_(double const * const in, const Key& key);
  void copyToSub_(double const * const in, const Key& key);
  void copyFromSurf_(double * const out, const Key& key) const;
  void copyFromSub_(double * const out, const Key& key) const;

  void initPressureFromWC_(double const * const elm_water_content);
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
  Key poro_key_;
  Key root_frac_key_;
  Key pot_evap_key_;
  Key pot_trans_key_;
  Key pot_infilt_key_;
  Key pd_key_;
  Key wtd_key_;
  Key pres_key_;
  Key wc_key_;
  Key pc_key_;
  Key sat_key_;
  //Key sat_gas_key_;
  //Key sat_ice_key_;
  Key infilt_key_;
  Key trans_key_;
  Key evap_key_;
  Key total_trans_key_;
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




