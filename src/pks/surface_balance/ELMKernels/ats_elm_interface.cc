/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

*/

#include <cmath>

#include "pk_helpers.hh"
#include "elm_kokkos_interface.hh"

#define NUM_LC_CLASSES 18

namespace Amanzi {
namespace LandPhysics {

  Kernels::Kernels(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& solution): {}

} // namespace LandPhysics
} // namespace ATS
