/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  Transport PK

*/

#ifndef AMANZI_ATS_SEDIMENTTRANSPORT_PK_HH_
#define AMANZI_ATS_SEDIMENTTRANSPORT_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Explicit_TI_FnBase.hh"
//#include "MaterialProperties.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "ReconstructionCellLinear.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "PK_PhysicalExplicit.hh"
#include "DenseVector.hh"

#include <string>

// Transport
//#include "TransportDomainFunction.hh"
//#include "SedimentTransportDefs.hh"
#include "transport_ats.hh"
#include "pk_physical_default.hh"

/* ******************************************************************
The transport PK receives a reduced (optional) copy of a physical
state at time n and returns a different state at time n+1.

Unmodified physical quantaties in the returned state are the smart
pointers to the original variables.
****************************************************************** */

namespace Amanzi {
namespace Transport {

class SedimentTransport_PK : public Transport_ATS {
 public:
  SedimentTransport_PK(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& glist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln);

  ~SedimentTransport_PK() = default;

  void parseParameterList() override;
  void Initialize() override;

 protected:

  void SetupTransport_() override;
  void AddSourceTerms_(double t0, double t1,
        Epetra_MultiVector& conserve_qty,
        int n0,
        int n1)  override;;

  
  Key sd_trapping_key_, sd_settling_key_, sd_erosion_key_, horiz_mixing_key_,  sd_organic_key_;
  Key elevation_increase_key_;
  Key biomass_key_, porosity_key_;
  Key plant_area_key_, stem_diameter_key_, stem_height_key_, stem_density_key_;

 private:
  double sediment_density_;


 private:
  // factory registration
  static RegisteredPKFactory<SedimentTransport_PK> reg_;
  
};

} // namespace Transport
} // namespace Amanzi

#endif
