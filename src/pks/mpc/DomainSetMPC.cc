/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A DomainSet coupler, couples a bunch of domains of the same structure.

------------------------------------------------------------------------- */

#include "DomainSetMPC.hh"

namespace Amanzi {

DomainSetMPC::DomainSetMPC(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution)
    : MPC<PK>(pk_tree, global_list, S, solution),
      PK(pk_tree, global_list, S, solution)
{
  // grab the list of subpks
  auto subpks = this->plist_->template get<Teuchos::Array<std::string> >("PKs order");
  std::string dsname = subpks.back();
  subpks.pop_back();

  KeyTriple triple;
  bool is_ds = Keys::splitDomainSet(dsname, triple);
  if (!is_ds || !S->HasDomainSet(std::get<0>(triple))) {
    Errors::Message msg;
    msg << "DomainSetMPC: \"" << dsname << "\" should be a domain-set of the form DOMAIN_*-PK_NAME";
    Exceptions::amanzi_throw(msg);
  }

  // add for the various sub-pks based on IDs
  auto ds = S->GetDomainSet(std::get<0>(triple));
  for (auto& name_id : *ds) {
    subpks.push_back(Keys::getKey(name_id.first, std::get<2>(triple)));
  }
  this->plist_->template set("PKs order", subpks);

  // construct the sub-PKs on COMM_SELF
  MPC<PK>::init_(S, getCommSelf());
}


// must communicate dts since columns are serial
double DomainSetMPC::get_dt() {
  double dt = 1.0e99;
  for (const auto& pk : sub_pks_) {
    dt = std::min<double>(dt,pk->get_dt());
  }
  
  double dt_local = dt;
  solution_->Comm()->MinAll(&dt_local, &dt, 1);
  
  return dt;
}

// -----------------------------------------------------------------------------
// Set timestep for sub PKs 
// -----------------------------------------------------------------------------
void DomainSetMPC::set_dt( double dt) {
  for (const auto& pk : sub_pks_) {
    pk->set_dt(dt);
  }

};


//-------------------------------------------------------------------------------------
// Semi coupled thermal hydrology
bool 
DomainSetMPC::AdvanceStep(double t_old, double t_new, bool reinit) {
  int nfailed = 0;
  for (const auto& pk : sub_pks_) {
    bool fail = pk->AdvanceStep(t_old, t_new, reinit);
    if (fail) {
      nfailed++;
      break;
    }
  }

  int nfailed_global(0);
  solution_->Comm()->SumAll(&nfailed, &nfailed_global, 1);
  if (nfailed_global) return true;
  return false;
}
  

} // namespace Amanzi
