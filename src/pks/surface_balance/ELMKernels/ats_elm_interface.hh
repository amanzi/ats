
#ifndef ATS_ELM_KERNELS_INTERFACE_HH_
#define ATS_ELM_KERNELS_INTERFACE_HH_

#include <cstdint>
#include <vector>

#include "Epetra_MultiVector.h"

#define NUM_LC_CLASSES 18

namespace ATS {
namespace LandPhysics {

  class Kernels {
  public:
    Kernels();
    Kernels(Teuchos::ParameterList& pk_tree,
           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
           const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

protected:

  Teuchos::RCP<ELM::ELMInterface> elm_kernels_;

  }; // Kernels

} // namespace LandPhysics
} // namespace ATS
