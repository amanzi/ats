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

#include "mpc_reactivetransport.hh"
#include "PK_Helpers.hh"
#include "chem_pk_helpers.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCReactiveTransport::MPCReactiveTransport(Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                           const Teuchos::RCP<State>& S,
                                           const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_list, S, soln),
    WeakMPC(pk_tree, global_list, S, soln) {}


void
MPCReactiveTransport::parseParameterList()
{
  // 0 is chem, 1 is transport
  cast_sub_pks_();

  domain_ = getSubPKPlist_(1)->get<std::string>("domain name", "domain");
  tcc_key_ = Keys::readKey(*getSubPKPlist_(0), domain_, "primary variable",
                           "total_component_concentration");
  mol_frac_key_ = Keys::readKey(*getSubPKPlist_(1), domain_, "primary variable",
                           "molar_ratio");
  if (tcc_key_ == mol_frac_key_) {
    Errors::Message msg;
    msg << "Chemistry and Transport may not be given the same primary variable name (\"" << tcc_key_
        << "\") -- rename one or the other.";
    Exceptions::amanzi_throw(msg);
  }

  mol_dens_key_ = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");

  // Only chemistry PK needs to set the IC -- we use chemistry to process geochemical
  // initial conditions -- but the "initial conditions" list must be present
  // for all PK_Physical PKs, so just touch it here to make sure it exists.
  getSubPKPlist_(1)->sublist("initial conditions");

  // this MPC need access to both primary variables
  getSubPKPlist_(0)->set<std::string>("primary variable password", name_);
  getSubPKPlist_(1)->set<std::string>("primary variable password", name_);

  WeakMPC::parseParameterList();
}


void
MPCReactiveTransport::Setup()
{
  // must Setup transport first to get alias for saturation, etc set up correctly
  transport_pk_->Setup();
  chemistry_pk_->Setup();

  requireEvaluatorAtNext(mol_dens_key_, tag_next_, *S_)
    .SetMesh(S_->GetMesh(domain_))
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}

void
MPCReactiveTransport::cast_sub_pks_()
{
  // cast and call parse on chemistry
  chemistry_pk_ = Teuchos::rcp_dynamic_cast<AmanziChemistry::Alquimia_PK>(sub_pks_[0]);
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);

  // now chem engine is set and we can hand it to transport
  transport_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(sub_pks_[1]);
  AMANZI_ASSERT(transport_pk_ != Teuchos::null);
  transport_pk_->setChemEngine(chemistry_pk_);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void
MPCReactiveTransport::Initialize()
{
  // initialize chemistry, including geochemical ICs
  chemistry_pk_->Initialize();
  AMANZI_ASSERT(S_->GetRecord(tcc_key_, tag_next_).initialized());

  // Compute mol frac from concentration
  //
  // NOTE: this requires that Reactive-Transport is done last, or at least
  // after the density of water can be evaluated.  This could be problematic
  // for, e.g., salinity intrusion problems where water density is a function
  // of concentration itself, but should work for all other problems?
  convertConcentrationToMolFrac(*S_, {tcc_key_, tag_next_},
          {mol_frac_key_, tag_next_}, {mol_dens_key_, tag_next_}, name());
  S_->GetRecordW(mol_frac_key_, tag_next_, name()).set_initialized();

  // initialize transport
  transport_pk_->Initialize();
}



// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
MPCReactiveTransport::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // First we do a transport step.
  bool fail = transport_pk_->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // move from mol_frac@next to tcc@current
  convertMolFracToConcentration(*S_, {mol_frac_key_, tag_next_},
          {tcc_key_, tag_current_}, {mol_dens_key_, tag_next_}, name());

  // Next to chemistry step
  fail = chemistry_pk_->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // move from tcc@next to mol_frac@next
  convertConcentrationToMolFrac(*S_, {tcc_key_, tag_next_},
          {mol_frac_key_, tag_next_}, {mol_dens_key_, tag_next_}, name());
  return fail;
};


} // namespace Amanzi
