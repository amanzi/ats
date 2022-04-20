/*
  This is the mpc_pk component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

*/

#include "mpc_coupled_transport.hh"


namespace Amanzi {

CoupledTransport_PK::CoupledTransport_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) :
  PK(pk_tree_or_fe_list, global_list, S, soln),
  WeakMPC(pk_tree_or_fe_list, global_list, S, soln)
{
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = global_list->sublist("verbose object");
  vo_ =  Teuchos::rcp(new VerboseObject("Coupled TransportPK", vlist));
  name_ = Keys::cleanPListName(plist_->name());

  Teuchos::Array<std::string> pk_order = plist_->get<Teuchos::Array<std::string> >("PKs order");

  for (int i=0; i<2; i++){
    Teuchos::RCP<Teuchos::ParameterList> list = Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), pk_order[i]);
    Key dom_name = list->get<std::string>("domain name", "domain");
    if (dom_name == "domain"){
      subsurface_transport_list_ = list;
      subsurface_name_ = dom_name;
      mesh_ = S->GetMesh(subsurface_name_);
      subsurf_id_ = i;
      subsurface_tcc_key_ = Keys::readKey(*list, dom_name, "concentration", "total_component_concentration");
    }else{
      surface_transport_list_ = list;
      surface_name_ = dom_name;
      surf_mesh_ = S->GetMesh(surface_name_);
      surf_id_ = i;
      surface_tcc_key_ = Keys::readKey(*list, dom_name, "concentration", "total_component_concentration");
    }
  }

  subsurface_flux_key_ =  plist_->get<std::string>("flux_key",
          Keys::getKey(subsurface_name_, "water_flux"));

  surface_flux_key_ =  plist_->get<std::string>("flux_key",
          Keys::getKey(surface_name_, "water_flux"));

}

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double CoupledTransport_PK::get_dt()
{
  double surf_dt = sub_pks_[surf_id_]->get_dt();
  double subsurf_dt = sub_pks_[subsurf_id_]->get_dt();

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH){
    *vo_->os() << surface_name_<<"  transport dt = " << surf_dt << std::endl
               << subsurface_name_<< " transport dt = " << subsurf_dt << std::endl;
  }
  double dt = std::min(surf_dt, subsurf_dt);
  set_dt(dt);
  return dt;
}


void CoupledTransport_PK::Setup(const Teuchos::Ptr<State>& S)
{
  WeakMPC::Setup(S);

  subsurf_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(sub_pks_[subsurf_id_]);
  AMANZI_ASSERT(subsurf_pk_ != Teuchos::null);
  surf_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(sub_pks_[surf_id_]);
  AMANZI_ASSERT(surf_pk_ != Teuchos::null);

  //SetupCouplingConditions();
}

int CoupledTransport_PK::num_aqueous_component()
{
  int num_aq_comp = subsurf_pk_->num_aqueous_component();
  if (num_aq_comp != surf_pk_->num_aqueous_component()){
    Errors::Message message("CoupledTransport_PK:: numbers aqueous component does not match.");
    Exceptions::amanzi_throw(message);
  }
  return num_aq_comp;
}

// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool CoupledTransport_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  double dt_MPC = S_inter_->final_time() - S_inter_->initial_time();


  AMANZI_ASSERT(subsurf_pk_ != Teuchos::null);
  AMANZI_ASSERT(surf_pk_ != Teuchos::null);

  // In the case of rain sources these rain source are mixed with solutes on the surface
  // to provide BC for subsurface domain
  surf_pk_->AdvanceStep(t_old, t_new, reinit);
  subsurf_pk_->AdvanceStep(t_old, t_new, reinit);

  const Epetra_MultiVector& surf_tcc = *S_inter_->GetFieldCopyData(surface_tcc_key_, "subcycling")->ViewComponent("cell",false);
  const Epetra_MultiVector& tcc = *S_inter_->GetFieldCopyData(subsurface_tcc_key_, "subcycling")->ViewComponent("cell",false);

  const std::vector<std::string>&  component_names_sub = subsurf_pk_->component_names();
  int num_components =  subsurf_pk_->num_aqueous_component();
  std::vector<double> mass_subsurface(num_components, 0.), mass_surface(num_components, 0.);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM){

    for (int i=0; i<num_components; i++){
      *vo_->os() << component_names_sub[i]<< ":";
      mass_subsurface[i] = Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(sub_pks_[subsurf_id_])
        ->ComputeSolute(tcc, i);
      mass_surface[i] = Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(sub_pks_[surf_id_])
        ->ComputeSolute(surf_tcc, i);
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() <<subsurface_name_<<" = " << mass_subsurface[i] << " mol";
      *vo_->os() <<", "<<surface_name_<<" = " << mass_surface[i]<< " mol";
      *vo_->os() <<", ToTaL =" << mass_surface[i]+mass_subsurface[i]<< " mol" <<std::endl;
    }
  }
  return fail;
}


/* *******************************************************************
 * Interpolate linearly in time between two values v0 and v1. The time
 * is measuared relative to value v0; so that v1 is at time dt. The
 * interpolated data are at time dt_int.
 ******************************************************************* */
void CoupledTransport_PK::InterpolateCellVector(
  const Epetra_MultiVector& v0, const Epetra_MultiVector& v1,
  double dt_int, double dt, Epetra_MultiVector& v_int)
{
  double a = dt_int / dt;
  double b = 1.0 - a;
  v_int.Update(b, v0, a, v1, 0.);
}

void CoupledTransport_PK::SetupCouplingConditions()
{
  Teuchos::ParameterList& bc_list = subsurface_transport_list_->sublist("boundary conditions").sublist("concentration");
  Teuchos::ParameterList& src_list = surface_transport_list_->sublist("source terms").sublist("component mass source");

  if (!bc_list.isSublist("BC coupling")){
    Teuchos::ParameterList& bc_coupling = bc_list.sublist("BC coupling");
    bc_coupling.set<std::string>("spatial distribution method", "domain coupling");
    bc_coupling.set<std::string>("submodel", "conserved quantity");
    std::vector<std::string> regs(1);
    regs[0] = "surface";
    bc_coupling.set<Teuchos::Array<std::string> >("regions", regs);
    Teuchos::ParameterList& tmp = bc_coupling.sublist("fields");
    tmp.set<std::string>("conserved_quantity_key", "surface-total_component_quantity");
    tmp.set<std::string>("field_out_key", surface_tcc_key_);
  }

  if (!src_list.isSublist("surface coupling")){
    Teuchos::ParameterList& src_coupling = src_list.sublist("surface coupling");
    src_coupling.set<std::string>("spatial distribution method", "domain coupling");
    src_coupling.set<std::string>("submodel", "rate");
    std::vector<std::string> regs = {"surface domain"};
    //regs[0] = surface_name_;
    src_coupling.set<Teuchos::Array<std::string> >("regions", regs);
    Teuchos::ParameterList& tmp = src_coupling.sublist("fields");
    tmp.set<std::string>("flux_key", subsurface_flux_key_);
    tmp.set<std::string>("copy_flux_key", "next_timestep");
    tmp.set<std::string>("field_in_key", surface_tcc_key_);
    tmp.set<std::string>("field_out_key", subsurface_tcc_key_);
    tmp.set<std::string>("copy_field_out_key", "default");
    tmp.set<std::string>("copy_field_in_key", "default");
  }
}

}  // namespace Amanzi
