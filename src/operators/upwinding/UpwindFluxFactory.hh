/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Factory for taking coefficients for div-grad operators from cells to faces.
#ifndef AMANZI_UPWINDING_FACTORY_HH_
#define AMANZI_UPWINDING_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "upwinding.hh"

namespace Amanzi {
namespace Operators {
namespace UpwindFactory {

Teuchos::RCP<Upwinding> Create(Teuchos::ParameterList& oplist,
                               State& S,
                               const std::string& pkname,
                               const Tag& tag,
                               const Key& flux_key);

} // namespace UpwindFactory
} // namespace Operators
} // namespace Amanzi

#endif
