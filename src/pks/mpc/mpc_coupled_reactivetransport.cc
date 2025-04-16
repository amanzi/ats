/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Process kernel for coupling of Transport PK and Chemistry PK.
*/

#include "mpc_coupled_reactivetransport.hh"
#include "chem_pk_helpers.hh"
#include "PK_Helpers.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCCoupledReactiveTransport::MPCCoupledReactiveTransport(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_list, S, soln),
    WeakMPC(pk_tree, global_list, S, soln)
{}

void
MPCCoupledReactiveTransport::parseParameterList()
{
  // tweak the sub-PK parameter lists
  auto coupled_pk_names = plist_->get<Teuchos::Array<std::string>>("PKs order");
  auto chem_names = getSubPKPlist_(0)->get<Teuchos::Array<std::string>>("PKs order");
  auto transport_names = getSubPKPlist_(1)->get<Teuchos::Array<std::string>>("PKs order");

  domain_ = pks_list_->sublist(transport_names[0]).get<std::string>("domain name", "domain");
  domain_surf_ = pks_list_->sublist(transport_names[1]).get<std::string>("domain name", "surface");

  mol_frac_key_ = Keys::readKey(pks_list_->sublist(transport_names[0]), domain_,
                           "primary variable", "molar_mixing_ratio");
  mol_frac_surf_key_ = Keys::readKey(pks_list_->sublist(transport_names[1]), domain_surf_,
                           "primary variable", "molar_mixing_ratio");

  tcc_key_ = Keys::readKey(pks_list_->sublist(chem_names[0]), domain_,
                           "primary variable", "total_component_concentration");
  tcc_surf_key_ = Keys::readKey(pks_list_->sublist(chem_names[1]), domain_surf_,
                           "primary variable", "total_component_concentration");

  mol_dens_key_ = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");
  mol_dens_surf_key_ =
    Keys::readKey(*plist_, domain_surf_, "surface molar density liquid", "molar_density_liquid");

  if (tcc_key_ == mol_frac_key_) {
    Errors::Message msg;
    msg << "Chemistry and Transport may not be given the same primary variable name (\"" << tcc_key_
        << "\") -- rename one or the other.";
    Exceptions::amanzi_throw(msg);
  }
  if (tcc_surf_key_ == mol_frac_surf_key_) {
    Errors::Message msg;
    msg << "Chemistry and Transport may not be given the same primary variable name (\"" << tcc_surf_key_
        << "\") -- rename one or the other.";
    Exceptions::amanzi_throw(msg);
  }

  // this MPC accesses chemistry and transport primary variables
  pks_list_->sublist(chem_names[0]).set<std::string>("primary variable password", name_);
  pks_list_->sublist(chem_names[1]).set<std::string>("primary variable password", name_);
  pks_list_->sublist(transport_names[0]).set<std::string>("primary variable password", name_);
  pks_list_->sublist(transport_names[1]).set<std::string>("primary variable password", name_);

  // Only reaction PKs set IC, but all need the list to be present in all PK_Physical_Default.
  pks_list_->sublist(transport_names[0]).sublist("initial conditions");
  pks_list_->sublist(transport_names[1]).sublist("initial conditions");

  cast_sub_pks_();

  coupled_chemistry_pk_->parseParameterList();
#ifdef ALQUIMIA_ENABLED
  transport_pk_->setChemEngine(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_));
  transport_pk_surf_->setChemEngine(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_surf_));
#endif
  coupled_transport_pk_->parseParameterList();
}


void
MPCCoupledReactiveTransport::Setup()
{
  // must Setup transport first to get alias for saturation, etc set up correctly
  coupled_transport_pk_->Setup();
  coupled_chemistry_pk_->Setup();

  requireEvaluatorAtNext(mol_dens_key_, tag_next_, *S_)
    .SetMesh(S_->GetMesh(domain_))
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  requireEvaluatorAtNext(mol_dens_surf_key_, tag_next_, *S_)
    .SetMesh(S_->GetMesh(domain_surf_))
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}


void
MPCCoupledReactiveTransport::cast_sub_pks_()
{
  coupled_transport_pk_ = Teuchos::rcp_dynamic_cast<MPCCoupledTransport>(sub_pks_[1]);
  AMANZI_ASSERT(coupled_transport_pk_ != Teuchos::null);

  coupled_chemistry_pk_ = Teuchos::rcp_dynamic_cast<WeakMPC>(sub_pks_[0]);
  AMANZI_ASSERT(coupled_chemistry_pk_ != Teuchos::null);

  transport_pk_ =
    Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(coupled_transport_pk_->get_subpk(0));
  AMANZI_ASSERT(transport_pk_ != Teuchos::null);
  transport_pk_surf_ =
    Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(coupled_transport_pk_->get_subpk(1));
  AMANZI_ASSERT(transport_pk_surf_ != Teuchos::null);

  chemistry_pk_ =
    Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(coupled_chemistry_pk_->get_subpk(0));
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);
  chemistry_pk_surf_ =
    Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(coupled_chemistry_pk_->get_subpk(1));
  AMANZI_ASSERT(chemistry_pk_surf_ != Teuchos::null);

  AMANZI_ASSERT(transport_pk_->domain() == chemistry_pk_->domain());
  AMANZI_ASSERT(transport_pk_surf_->domain() == chemistry_pk_surf_->domain());
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void
MPCCoupledReactiveTransport::Initialize()
{
  // initialize chemistry, including geochemical ICs
  coupled_chemistry_pk_->Initialize();
  AMANZI_ASSERT(S_->GetRecord(tcc_key_, tag_next_).initialized());
  AMANZI_ASSERT(S_->GetRecord(tcc_surf_key_, tag_next_).initialized());

  // Compute mol frac from concentration
  //
  // NOTE: this requires that Reactive-Transport is done last, or at least
  // after the density of water can be evaluated.  This could be problematic
  // for, e.g., salinity intrusion problems where water density is a function
  // of concentration itself, but should work for all other problems?
  convertConcentrationToMolFrac(*S_, {tcc_key_, tag_next_},
          {mol_frac_key_, tag_next_}, {mol_dens_key_, tag_next_}, name());
  S_->GetRecordW(mol_frac_key_, tag_next_, name()).set_initialized();

  convertConcentrationToMolFrac(*S_, {tcc_surf_key_, tag_next_},
          {mol_frac_surf_key_, tag_next_}, {mol_dens_surf_key_, tag_next_}, name());
  S_->GetRecordW(mol_frac_surf_key_, tag_next_, name()).set_initialized();

  coupled_transport_pk_->Initialize();
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
MPCCoupledReactiveTransport::AdvanceStep(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // First we do a transport step.
  bool fail = coupled_transport_pk_->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // move from mol_frac@next to tcc@current
  convertMolFracToConcentration(*S_, {mol_frac_key_, tag_next_},
          {tcc_key_, tag_current_}, {mol_dens_key_, tag_next_}, name());
  convertMolFracToConcentration(*S_, {mol_frac_surf_key_, tag_next_},
          {tcc_surf_key_, tag_current_}, {mol_dens_surf_key_, tag_next_}, name());

  // Next to chemistry step
  fail = coupled_chemistry_pk_->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // move from tcc@next to mol_frac@next
  convertConcentrationToMolFrac(*S_, {tcc_key_, tag_next_},
          {mol_frac_key_, tag_next_}, {mol_dens_key_, tag_next_}, name());
  convertConcentrationToMolFrac(*S_, {tcc_surf_key_, tag_next_},
          {mol_frac_surf_key_, tag_next_}, {mol_dens_surf_key_, tag_next_}, name());
  return fail;
};


} // namespace Amanzi
