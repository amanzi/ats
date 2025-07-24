/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Source term evaluator for enthalpy of mass source.

*/

#include "advected_energy_source_evaluator.hh"

namespace Amanzi {
namespace Energy {

// constructor format for all derived classes
AdvectedEnergySourceEvaluator::AdvectedEnergySourceEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  InitializeFromPlist_();
}

Teuchos::RCP<Evaluator>
AdvectedEnergySourceEvaluator::Clone() const
{
  return Teuchos::rcp(new AdvectedEnergySourceEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
AdvectedEnergySourceEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& int_enth =
    *S.GetPtr<CompositeVector>(internal_enthalpy_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& ext_enth =
    *S.GetPtr<CompositeVector>(external_enthalpy_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& water_source =
    *S.GetPtr<CompositeVector>(water_source_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& cv =
    *S.GetPtr<CompositeVector>(cell_vol_key_, tag)->ViewComponent("cell", false);

  Epetra_MultiVector& res = *result[0]->ViewComponent("cell", false);

  if (source_units_ == SOURCE_UNITS_METERS_PER_SECOND) {
    const Epetra_MultiVector& int_dens =
      *S.GetPtr<CompositeVector>(internal_density_key_, tag)->ViewComponent("cell", false);
    const Epetra_MultiVector& ext_dens =
      *S.GetPtr<CompositeVector>(external_density_key_, tag)->ViewComponent("cell", false);

    unsigned int ncells = res.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) {
      if (water_source[0][c] > 0.) { // positive indicates increase of water in surface
        // upwind, take external values
        res[0][c] = water_source[0][c] * ext_dens[0][c] * ext_enth[0][c];
      } else {
        // upwind, take internal values
        res[0][c] = water_source[0][c] * int_dens[0][c] * int_enth[0][c];
      }
    }

  } else {
    unsigned int ncells = res.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) {
      if (water_source[0][c] > 0.) { // positive indicates increase of water in surface
        // upwind, take external values
        res[0][c] = water_source[0][c] * ext_enth[0][c];
      } else {
        // upwind, take internal values
        res[0][c] = water_source[0][c] * int_enth[0][c];
        //std::cout << "Advected E-source air-surf = " << water_source[0][c] << " * " << int_enth[0][c] << " = " << res[0][c] << std::endl;
      }
    }
  }

  if (source_units_ == SOURCE_UNITS_MOLS_PER_SECOND) {
    unsigned int ncells = res.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) {
      res[0][c] /= cv[0][c];
    }
  }

  if (include_conduction_) {
    const Epetra_MultiVector& cond =
      *S.GetPtr<CompositeVector>(conducted_source_key_, tag)->ViewComponent("cell", false);
    unsigned int ncells = res.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) {
      res[0][c] += cond[0][c];
    }
  }
}

void
AdvectedEnergySourceEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  if (include_conduction_ && wrt_key == conducted_source_key_) {
    *result[0]->ViewComponent("cell", false) =
      *S.GetPtr<CompositeVector>(cell_vol_key_, tag)->ViewComponent("cell", false);
  } else {
    result[0]->PutScalar(0.);
  }
}

void
AdvectedEnergySourceEvaluator::InitializeFromPlist_()
{
  std::string domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  internal_enthalpy_key_ = Keys::readKey(plist_, domain, "internal enthalpy", "enthalpy");
  dependencies_.insert(KeyTag{ internal_enthalpy_key_, tag });

  external_enthalpy_key_ =
    Keys::readKey(plist_, domain, "external enthalpy", "water_source_enthalpy");
  dependencies_.insert(KeyTag{ external_enthalpy_key_, tag });

  water_source_key_ = Keys::readKey(plist_, domain, "water source", "water_source");
  dependencies_.insert(KeyTag{ water_source_key_, tag });

  // this handles both surface fluxes (in m/s) and subsurface fluxes (in mol/s)
  std::string source_units = plist_.get<std::string>("water source units");
  if (source_units == "mol s^-1") {
    source_units_ = SOURCE_UNITS_MOLS_PER_SECOND;
  } else if (source_units == "m s^-1") {
    source_units_ = SOURCE_UNITS_METERS_PER_SECOND;
  } else if (source_units == "mol m^-2 s^-1" || source_units == "mol m^-3 s^-1") {
    source_units_ = SOURCE_UNITS_MOLS_PER_SECOND_PER_METERSD;
  } else {
    Errors::Message message;
    message
      << "AdvectedEnergySourceEvaluator: " << my_keys_.front().first << ": invalid units \""
      << source_units
      << "\" for \"mass source units\", valid are \"mol s^-1\", \"m s^-1\", \"mol m^-2 s^-1\","
      << " and \"mol m^-3 s^-1\".";
    Exceptions::amanzi_throw(message);
  }

  if (source_units_ == SOURCE_UNITS_METERS_PER_SECOND) {
    internal_density_key_ =
      Keys::readKey(plist_, domain, "internal molar density", "molar_density_liquid");
    dependencies_.insert(KeyTag{ internal_density_key_, tag });

    external_density_key_ =
      Keys::readKey(plist_, domain, "external molar density", "source_molar_density");
    dependencies_.insert(KeyTag{ external_density_key_, tag });
  }

  // this enables the addition of provided diffusive fluxes as well
  include_conduction_ = plist_.get<bool>("include conduction");
  if (include_conduction_) {
    conducted_source_key_ = plist_.get<std::string>(
      "conducted energy source key", Keys::getKey(domain, "conducted_energy_source"));
    dependencies_.insert(KeyTag{ conducted_source_key_, tag });
  }

  cell_vol_key_ = Keys::readKey(plist_, domain, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cell_vol_key_, tag });
}


} // namespace Energy
} // namespace Amanzi
