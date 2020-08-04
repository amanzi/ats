/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "TransportSourceFunction_Alquimia.hh"

#ifdef ALQUIMIA_ENABLED

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Constructor of BCs for Alquimia.
****************************************************************** */
TransportSourceFunction_Alquimia::TransportSourceFunction_Alquimia(
    const Teuchos::ParameterList& plist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
    Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine)
    : mesh_(mesh),
      chem_pk_(chem_pk),
      chem_engine_(chem_engine)
{
  // Check arguments.
  if (chem_engine_ != Teuchos::null) {
    chem_engine_->InitState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
    chem_engine_->GetPrimarySpeciesNames(tcc_names_);
  } else {
    Errors::Message msg;
    msg << "Geochemistry is off, but a geochemical condition was requested.";
    Exceptions::amanzi_throw(msg); 
  }

  // Get the regions assigned to this geochemical condition. We do not
  // check for region overlaps here, since a better way is to derive from 
  // the generic mesh function.
  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();
  std::vector<double> times = plist.get<Teuchos::Array<double> >("times").toVector();
  std::vector<std::string> conditions = plist.get<Teuchos::Array<std::string> >("geochemical conditions").toVector();

  // Function of geochemical conditions and the associates regions.
  f_ = Teuchos::rcp(new FunctionTabularString(times, conditions));
  Init_(regions);
}


/* ******************************************************************
* Delegating destructor.
****************************************************************** */
TransportSourceFunction_Alquimia::~TransportSourceFunction_Alquimia()
{
  chem_engine_->FreeState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
}


/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void TransportSourceFunction_Alquimia::Init_(const std::vector<std::string>& regions)
{
  for (int i = 0; i < regions.size(); ++i) {
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(regions[i], AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL, &block);
    int nblock = block.size();

    // Now get the cells that are attached to these faces.
    for (int n = 0; n < nblock; ++n) {
      int c = block[n];
      value_[c].resize(chem_engine_->NumPrimarySpecies());
    }
  }
}


/* ******************************************************************
* Evaluate values at time.
****************************************************************** */
void TransportSourceFunction_Alquimia::Compute(double t_old, double t_new) 
{
  std::string cond_name = (*f_)(t_new);

  //if (src_liq_data_==Teuchos::null)

  // Loop over sides and evaluate values.
  for (auto it = begin(); it != end(); ++it) {
    int cell = it->first; 

    // Dump the contents of the chemistry state into our Alquimia containers.
    chem_pk_->CopyToAlquimia(cell, alq_mat_props_, alq_state_, alq_aux_data_);

    double liquid_vol_L = alq_mat_props_.volume * (*src_liq_data_)[0][cell] * 1000 ;
     
    // Enforce the condition.
    chem_engine_->EnforceCondition(cond_name, t_new, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

    // Move the concentrations into place.
    std::vector<double>& values = it->second;
    

    for (int i = 0; i < values.size(); i++) {

      double moles = alq_state_.total_mobile.data[i] * liquid_vol_L; // conversion from mol/L to mol
      values[i] = moles/(alq_mat_props_.volume);
    }    
  }

}

}  // namespace Transport
}  // namespace Amanzi

#endif

