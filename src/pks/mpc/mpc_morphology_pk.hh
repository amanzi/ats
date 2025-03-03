/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  PK for coupling of surface and subsurface transport PKs
*/

#ifndef ATS_AMANZI_MORPHOLOGY_PK_HH_
#define ATS_AMANZI_MORPHOLOGY_PK_HH_

#include "Teuchos_RCP.hpp"

#include "mpc_subcycled.hh"
#include "pk_physical_bdf_default.hh"
#include "PK.hh"
#include "Debugger.hh"

namespace Amanzi {

class Morphology_PK : public MPCSubcycled {
 public:
  Morphology_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& soln);
  ~Morphology_PK() {}

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  //virtual void set_dt(double dt);
  virtual void Setup();
  virtual void Initialize();

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  virtual void CommitStep(double t_old, double t_new, const Tag& tag);

  std::string name() { return name_; }

 protected:
  void Initialize_MeshVertices_(const Teuchos::Ptr<State>& S,
                                Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                Key vert_field_key, const Amanzi::Tag field_tag);

  void Update_MeshVertices_(const Teuchos::Ptr<State>& S, const Amanzi::Tag field_tag);

  void FlowAnalyticalSolution_(const Teuchos::Ptr<State>& S, double time);

  Key domain_, domain_3d_, domain_ss_;
  Key vertex_coord_key_, vertex_coord_key_3d_, vertex_coord_key_ss_;
  Key elevation_increase_key_;

  Teuchos::RCP<Epetra_MultiVector> dz_accumul_;

  Teuchos::RCP<PK_BDF_Default> flow_pk_;
  Teuchos::RCP<PK> sed_transport_pk_;

  double master_dt_, slave_dt_;
  double dt_MPC_, dt_sample_;
  double MSF_; // morphology scaling factor

  Teuchos::RCP<AmanziMesh::Mesh> mesh_, mesh_3d_, mesh_ss_;
  Teuchos::RCP<EvaluatorPrimaryCV> deform_eval_;
  Key erosion_rate_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> flow_db_;
  Teuchos::RCP<Debugger> trans_db_;

  // factory registration
  static RegisteredPKFactory<Morphology_PK> reg_;
};

} // namespace Amanzi
#endif
