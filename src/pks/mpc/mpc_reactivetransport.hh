/*
  This is the mpc_pk component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling of Transport_PK and Chemistry_PK.
*/


#ifndef AMANZI_REACTIVETRANSPORT_PK_ATS_HH_
#define AMANZI_REACTIVETRANSPORT_PK_ATS_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "PK.hh"
#include "transport_ats.hh"
#include "Chemistry_PK.hh"
#include "PK_MPCAdditive.hh"

namespace Amanzi {

class MPCReactiveTransport : public WeakMPC {
 public:
  MPCReactiveTransport(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt() override;
  virtual void set_dt(double dt) override;

  // -- standard PK functions
  virtual void Setup() override;
  virtual void Initialize() override;
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  virtual void CommitStep(double t_old, double t_new, ) override;

  void ConvertConcentrationToAmanzi(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                                    const Epetra_MultiVector& mol_den,
                                    const Epetra_MultiVector& tcc_ats,
                                    Epetra_MultiVector& tcc_amanzi);

  void ConvertConcentrationToATS(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                                 const Epetra_MultiVector& mol_den,
                                 const Epetra_MultiVector& tcc_amanzi,
                                 Epetra_MultiVector& tcc_ats);

  bool AdvanceChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                        const Epetra_MultiVector& mol_den,
                        Teuchos::RCP<Epetra_MultiVector> tcc_copy,
                        double t_old, double t_new, bool reinit = false);


  std::string name() { return "reactive transport"; }



protected:
  virtual void cast_sub_pks_();

protected:
  int transport_pk_index_, chemistry_pk_index_;
  Teuchos::RCP<Teuchos::ParameterList> rt_pk_list_;
  bool chem_step_succeeded_;
  bool transport_subcycling_;
  double dTtran_, dTchem_;

  Key domain_;
  Key tcc_key_;
  Key mol_den_key_;

  Teuchos::RCP<Teuchos::Time> chem_timer_;
  Teuchos::RCP<Teuchos::Time> alquimia_timer_;

private:

  // storage for the component concentration intermediate values
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration_stor_;
  Teuchos::RCP<Transport::Transport_ATS> tranport_pk_;
  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chemistry_pk_;

  // factory registration
  static RegisteredPKFactory<MPCReactiveTransport> reg_;
};

}  // namespace Amanzi
#endif
