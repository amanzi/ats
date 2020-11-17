/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Source term evaluator for enthalpy of mass source.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "advected_energy_source_evaluator.hh"

namespace Amanzi {
namespace Energy {

// constructor format for all derived classes
AdvectedEnergySourceEvaluator::AdvectedEnergySourceEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  InitializeFromPlist_();
}

AdvectedEnergySourceEvaluator::AdvectedEnergySourceEvaluator(const AdvectedEnergySourceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    internal_enthalpy_key_(other.internal_enthalpy_key_),
    external_enthalpy_key_(other.external_enthalpy_key_),
    mass_source_key_(other.mass_source_key_),
    internal_density_key_(other.internal_density_key_),
    external_density_key_(other.external_density_key_),
    cell_vol_key_(other.cell_vol_key_),
    conducted_source_key_(other.conducted_source_key_),
    include_conduction_(other.include_conduction_),
    source_units_(other.source_units_)
{}

Teuchos::RCP<FieldEvaluator>
AdvectedEnergySourceEvaluator::Clone() const {
  return Teuchos::rcp(new AdvectedEnergySourceEvaluator(*this));
}

// Required methods from SecondaryVariableFieldEvaluator
void
AdvectedEnergySourceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  const Epetra_MultiVector& int_enth = *S->GetFieldData(internal_enthalpy_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& ext_enth = *S->GetFieldData(external_enthalpy_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& mass_source = *S->GetFieldData(mass_source_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S->GetFieldData(cell_vol_key_)
      ->ViewComponent("cell",false);

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  
  if (source_units_ == SOURCE_UNITS_METERS_PER_SECOND) {
    const Epetra_MultiVector& int_dens = *S->GetFieldData(internal_density_key_)
                                         ->ViewComponent("cell",false);
    const Epetra_MultiVector& ext_dens = *S->GetFieldData(external_density_key_)
                                         ->ViewComponent("cell",false);

    unsigned int ncells = res.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      if (mass_source[0][c] > 0.) { // positive indicates increase of water in surface
        // upwind, take external values
        res[0][c] = mass_source[0][c] * ext_dens[0][c] * ext_enth[0][c];
      } else {
        // upwind, take internal values
        res[0][c] = mass_source[0][c] * int_dens[0][c] * int_enth[0][c];
      }
    }

  } else {
    unsigned int ncells = res.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      if (mass_source[0][c] > 0.) { // positive indicates increase of water in surface
        // upwind, take external values
        res[0][c] = mass_source[0][c] * ext_enth[0][c];
      } else {
        // upwind, take internal values
        res[0][c] = mass_source[0][c] * int_enth[0][c];
        //std::cout << "Advected E-source air-surf = " << mass_source[0][c] << " * " << int_enth[0][c] << " = " << res[0][c] << std::endl;
      }
    }
  }

  if (source_units_ == SOURCE_UNITS_MOLS_PER_SECOND) {
    unsigned int ncells = res.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      res[0][c] /= cv[0][c];
    }
  }
  

  if (include_conduction_) {
    const Epetra_MultiVector& cond = *S->GetFieldData(conducted_source_key_)
        ->ViewComponent("cell",false);
    unsigned int ncells = res.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      res[0][c] += cond[0][c];
    }
  }
}

void
AdvectedEnergySourceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  if (include_conduction_ && wrt_key == conducted_source_key_) {
    *result->ViewComponent("cell",false) = *S->GetFieldData(cell_vol_key_)
        ->ViewComponent("cell",false);
  } else {
    result->PutScalar(0.);
  }
}

void
AdvectedEnergySourceEvaluator::InitializeFromPlist_() {
  std::string domain = Keys::getDomain(my_key_);

  internal_enthalpy_key_ = Keys::readKey(plist_, domain, "internal enthalpy", "enthalpy");
  external_enthalpy_key_ = Keys::readKey(plist_, domain, "external enthalpy", "mass_source_enthalpy");
  mass_source_key_ = Keys::readKey(plist_, domain, "mass source", "mass_source");

  dependencies_.insert(internal_enthalpy_key_);
  dependencies_.insert(external_enthalpy_key_);
  dependencies_.insert(mass_source_key_);

  // this handles both surface fluxes (in m/s) and subsurface fluxes (in mol/s)
  std::string source_units = plist_.get<std::string>("mass source units");
  if (source_units == "mol s^-1") {
    source_units_ = SOURCE_UNITS_MOLS_PER_SECOND;
  } else if (source_units == "m s^-1") {
    source_units_ = SOURCE_UNITS_METERS_PER_SECOND;
  } else if (source_units == "mol m^-2 s^-1" ||
             source_units == "mol m^-3 s^-1") {
    source_units_ = SOURCE_UNITS_MOLS_PER_SECOND_PER_METERSD;
  } else {
    Errors::Message message;
    message << "AdvectedEnergySourceEvaluator: "
            << my_key_
            << ": invalid units \""
            << source_units
            << "\" for \"mass source units\", valid are \"mol s^-1\", \"m s^-1\", \"mol m^-2 s^-1\","
            << " and \"mol m^-3 s^-1\".";
    Exceptions::amanzi_throw(message);
  }

  if (source_units_ == SOURCE_UNITS_METERS_PER_SECOND) {
    internal_density_key_ = Keys::readKey(plist_, domain, "internal molar density", "molar_density_liquid");
    external_density_key_ = Keys::readKey(plist_, domain, "external molar density", "source_molar_density");
    dependencies_.insert(internal_density_key_);
    dependencies_.insert(external_density_key_);
  }

  // this enables the addition of provided diffusive fluxes as well
  include_conduction_ = plist_.get<bool>("include conduction");
  if (include_conduction_) {
    conducted_source_key_ = plist_.get<std::string>("conducted energy source key",
            Keys::getKey(domain, "conducted_energy_source"));
    dependencies_.insert(conducted_source_key_);
  }

  cell_vol_key_ = plist_.get<std::string>("cell volume key",
          Keys::getKey(domain, "cell_volume"));

}


} //namespace
} //namespace

