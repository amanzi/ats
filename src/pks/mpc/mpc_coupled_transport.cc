/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  This is the mpc_pk component of the Amanzi code.

*/

#include "mpc_coupled_transport.hh"

namespace Amanzi {

MPCCoupledTransport::MPCCoupledTransport(Teuchos::ParameterList& pk_tree,
                                         const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                         const Teuchos::RCP<State>& S,
                                         const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_list, S, soln), WeakMPC(pk_tree, global_list, S, soln)
{}


void
MPCCoupledTransport::Setup()
{
  name_ss_ = sub_pks_[0]->name();
  name_surf_ = sub_pks_[1]->name();

  // see bug amanzi/ats#125 this is probably backwards
  pk_ss_ = Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(sub_pks_[0]);
  pk_surf_ = Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(sub_pks_[1]);

  if (pk_ss_ == Teuchos::null || pk_surf_ == Teuchos::null) {
    Errors::Message msg("MPCCoupledTransport expects to only couple PKs of type \"transport ATS\"");
    Exceptions::amanzi_throw(msg);
  }

  SetupCouplingConditions_();
  WeakMPC::Setup();
}

// bug, see amanzi/ats#125
bool
MPCCoupledTransport::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = pk_surf_->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;
  fail = pk_ss_->AdvanceStep(t_old, t_new, reinit);
  return fail;
}


int
MPCCoupledTransport::get_num_aqueous_component()
{
  int num_aq_comp = pk_ss_->get_num_aqueous_component();
  if (num_aq_comp != pk_surf_->get_num_aqueous_component()) {
    Errors::Message msg("MPCCoupledTransport:: numbers aqueous component does not match.");
    Exceptions::amanzi_throw(msg);
  }
  return num_aq_comp;
}


void
MPCCoupledTransport::SetupCouplingConditions_()
{
  Key domain_ss = pks_list_->sublist(name_ss_).get<std::string>("domain name", "domain");
  Key domain_surf = pks_list_->sublist(name_surf_).get<std::string>("domain name", "surface");

  auto& bc_list =
    pks_list_->sublist(name_ss_).sublist("boundary conditions").sublist("concentration");
  auto& src_list =
    pks_list_->sublist(name_surf_).sublist("source terms").sublist("component mass source");

  Key ss_flux_key =
    Keys::readKey(pks_list_->sublist(name_ss_), domain_ss, "water flux", "water_flux");
  Key surf_flux_key =
    Keys::readKey(pks_list_->sublist(name_surf_), domain_surf, "water flux", "water_flux");

  Key ss_tcc_key = Keys::readKey(pks_list_->sublist(name_ss_), domain_ss, "concentration");
  Key surf_tcc_key = Keys::readKey(pks_list_->sublist(name_surf_), domain_surf, "concentration");
  Key surf_tcq_key = Keys::readKey(
    pks_list_->sublist(name_surf_), domain_surf, "conserved quantity", "total_component_quantity");

  if (!bc_list.isSublist("BC coupling")) {
    Teuchos::ParameterList& bc_coupling = bc_list.sublist("BC coupling");
    bc_coupling.set<std::string>("spatial distribution method", "domain coupling");
    bc_coupling.set<std::string>("submodel", "conserved quantity");
    std::vector<std::string> regs(1);
    regs[0] = "surface";
    bc_coupling.set<Teuchos::Array<std::string>>("regions", regs);
    Teuchos::ParameterList& tmp = bc_coupling.sublist("fields");
    tmp.set<std::string>("conserved quantity key", surf_tcq_key);
    tmp.set<std::string>("conserved quantity copy key", tag_next_.get());
    tmp.set<std::string>("external field key", surf_tcc_key);
    tmp.set<std::string>("external field copy key", tag_next_.get());
  }

  if (!src_list.isSublist("surface coupling")) {
    Teuchos::ParameterList& src_coupling = src_list.sublist("surface coupling");
    src_coupling.set<std::string>("spatial distribution method", "domain coupling");
    src_coupling.set<std::string>("submodel", "rate");
    std::vector<std::string> regs = { "surface domain" };
    //regs[0] = surface_name_;
    src_coupling.set<Teuchos::Array<std::string>>("regions", regs);
    Teuchos::ParameterList& tmp = src_coupling.sublist("fields");
    tmp.set<std::string>("flux key", ss_flux_key);
    // NOTE: FIXME --ETC amanzi/amanzi#646
    tmp.set<std::string>("flux copy key", Tags::NEXT.get());
    tmp.set<std::string>("conserved quantity key", surf_tcc_key);
    tmp.set<std::string>("conserved quantity copy key", tag_next_.get());
    tmp.set<std::string>("external field key", ss_tcc_key);
    tmp.set<std::string>("external field copy key", tag_next_.get());
  }
}

} // namespace Amanzi
