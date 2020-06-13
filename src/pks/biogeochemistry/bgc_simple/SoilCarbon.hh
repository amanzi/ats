/*

Soil carbon data structures class

Author: Ethan Coon (ecoon@lanl.gov)
        Chonggang Xu (cxu@lanl.gov)

Licencse: BSD
*/

#ifndef ATS_BGC_SOIL_CARBON_HH_
#define ATS_BGC_SOIL_CARBON_HH_

#include "Teuchos_RCP.hpp"

#include "SoilCarbonParameters.hh"

namespace Amanzi {
namespace BGC {

class SoilCarbon {

 public:

  SoilCarbon(const Teuchos::RCP<const SoilCarbonParameters> params_) :
      params(params_),
      SOM(params_->nPools),
      nPools(params_->nPools) {}

  SoilCarbon(const Teuchos::RCP<const SoilCarbonParameters> params_,
             double * som) :
      params(params_),
      SOM(View, som, params_->nPools),
      nPools(params_->nPools) {}

 public:
  int nPools;
  Epetra_SerialDenseVector SOM;
  Teuchos::RCP<const SoilCarbonParameters> params;

 private:
  SoilCarbon(const SoilCarbon& other); // not implemented
  SoilCarbon& operator=(const SoilCarbon& other); // not implemented
};


} // namespace
} // namespace

#endif
