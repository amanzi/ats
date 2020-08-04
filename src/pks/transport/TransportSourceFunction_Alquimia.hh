/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_SOURCE_FUNCTION_ALQUIMIA_HH_
#define AMANZI_TRANSPORT_SOURCE_FUNCTION_ALQUIMIA_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "FunctionTabularString.hh"
#include "TransportDomainFunction.hh"

#ifdef ALQUIMIA_ENABLED
#include "Alquimia_PK.hh"
#include "ChemistryEngine.hh"

namespace Amanzi {
namespace Transport {

class TransportSourceFunction_Alquimia : public TransportDomainFunction {
 public:
  TransportSourceFunction_Alquimia(const Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                   Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
                                   Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine);
  ~TransportSourceFunction_Alquimia();
  
  void Compute(double t_old, double t_new);

  void set_mol_dens_data_( Teuchos::Ptr<const Epetra_MultiVector> data) {mol_dens_data_ = data;}
  void set_liquid_src_data_( Teuchos::Ptr<const Epetra_MultiVector> data) {src_liq_data_ = data;}
  // require by the case class
  virtual std::string name() const { return "volume"; } 

 private:
  void Init_(const std::vector<std::string> &regions);


 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::Ptr<const Epetra_MultiVector> mol_dens_data_;
  Teuchos::Ptr<const Epetra_MultiVector> src_liq_data_;

  // string function of geochemical conditions
  Teuchos::RCP<FunctionTabularString> f_;

  // Chemistry state and engine.
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;

  // Containers for interacting with the chemistry engine.
  AlquimiaState alq_state_;
  AlquimiaProperties alq_mat_props_;
  AlquimiaAuxiliaryData alq_aux_data_;
  AlquimiaAuxiliaryOutputData alq_aux_output_;
};

}  // namespace Transport
}  // namespace Amanzi

#endif


#endif
