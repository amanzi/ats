/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#include <set>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "MultiFunction.hh"

#include "Mesh.hh"
#include "transport_ats.hh"

namespace Amanzi {
namespace Transport {

/* ****************************************************************
* Find place of the given component in a multivector.
**************************************************************** */
int
Transport_ATS::FindComponentNumber(const std::string component_name)
{
  int ncomponents = component_names_.size();
  for (int i = 0; i < ncomponents; i++) {
    if (component_names_[i] == component_name) return i;
  }
  Errors::Message msg("TransportExplicit_PK: component \"");
  msg << component_name << "\" was requested, but this is not a known component for this PK.";
  Exceptions::amanzi_throw(msg);
  return -1;
}

} // namespace Transport
} // namespace Amanzi
