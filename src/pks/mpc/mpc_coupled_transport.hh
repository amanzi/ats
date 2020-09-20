/*
  This is the mpc_pk component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  PK for coupling of surface and subsurface transport PKs
*/

#ifndef ATS_AMANZI_COUPLEDTRANSPORT_PK_HH_
#define ATS_AMANZI_COUPLEDTRANSPORT_PK_HH_

#include "Teuchos_RCP.hpp"

//#include "pk_mpcsubcycled_ats.hh"
#include "weak_mpc.hh"
#include "PK.hh"
#include "transport_ats.hh"

namespace Amanzi {

  class CoupledTransport_PK: public WeakMPC{

  public:
    CoupledTransport_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln);
    ~CoupledTransport_PK(){}

    // PK methods
    // -- dt is the minimum of the sub pks
    virtual double get_dt();
    //virtual void set_dt(double dt);
    virtual void Setup(const Teuchos::Ptr<State>& S);
    virtual void Initialize(const Teuchos::Ptr<State>& S);

    // -- advance each sub pk from t_old to t_new.
    virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);
    //virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

    std::string name() { return name_;}
    int num_aqueous_component();


  private:

    void InterpolateCellVector(const Epetra_MultiVector& v0, const Epetra_MultiVector& v1,
                               double dt_int, double dt, Epetra_MultiVector& v_int) ;

    void ComputeVolumeDarcyFlux(const Teuchos::Ptr<State>& S);

    Teuchos::RCP<const AmanziMesh::Mesh> mesh_, surf_mesh_;
    std::string passwd_;

    Teuchos::RCP<Teuchos::ParameterList> surface_transport_list_;
    Teuchos::RCP<Teuchos::ParameterList> subsurface_transport_list_;
    int subsurf_id_, surf_id_;

    Key subsurface_flux_key_, surface_flux_key_;
    Key surface_name_, subsurface_name_;
    Key mass_darcy_key, surf_mass_darcy_key;
    Key vol_darcy_key, surf_vol_darcy_key;
    Key mol_density_key, surf_mol_density_key;


    Teuchos::RCP<Transport::Transport_ATS> subsurf_pk_, surf_pk_;

    // factory registration
    static RegisteredPKFactory<CoupledTransport_PK> reg_;
};

}  // namespace Amanzi
#endif
