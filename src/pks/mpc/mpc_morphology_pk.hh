/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This PK couples surface flow and sediment transport and provides capability to model
  geomorphological changes of the surface elevation and slopes.
  The elevation change is defined as follows:

  .. math::

  \Delta Z = \frac{\delta t}{(1-\phi_s)\rho_s} (Q_t + Q_s - Q_e + Q_db)


  where :math:'\phi_s' is soil porosity and  :math:'\rho_s' - sediment density

  The surface and subsurface meshes has to be defined as 'deformable mesh'

  <Parameter name="deformable mesh" type="bool" value="true" />


*/

#ifndef ATS_AMANZI_MORPHOLOGY_PK_HH_
#define ATS_AMANZI_MORPHOLOGY_PK_HH_

#include "Teuchos_RCP.hpp"

#include "mpc_flow_transport.hh"
#include "pk_physical_bdf_default.hh"
#include "PK.hh"
#include "Debugger.hh"

namespace Amanzi {

class Morphology_PK : public MPCFlowTransport {
 public:
  Morphology_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  void parseParameterList() override;
  void Setup() override;
  void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

 protected:
  void Update_MeshVertices_(const Teuchos::Ptr<State>& S, const Tag& tag);

  Key domain_, domain_3d_, domain_ss_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_, mesh_3d_, mesh_ss_;

  Key vertex_coord_key_;
  Key elevation_increase_key_;
  Key dens_key_, pd_key_; // used in transport's subcycled quantites, need interpolations

  double MSF_; // morphology scaling factor

  // factory registration
  static RegisteredPKFactory<Morphology_PK> reg_;
};

} // namespace Amanzi
#endif
