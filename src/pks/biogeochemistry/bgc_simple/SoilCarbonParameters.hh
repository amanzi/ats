/*

Soil carbon parameters data structures class

Author: Ethan Coon (ecoon@lanl.gov)
        Chonggang Xu (cxu@lanl.gov)

Licencse: BSD
*/

#ifndef ATS_BGC_SOIL_CARBON_PARAMS_HH_
#define ATS_BGC_SOIL_CARBON_PARAMS_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace BGC {

class SoilCarbonParameters {

 public:

  SoilCarbonParameters(int npools, double percent_sand);
  SoilCarbonParameters(int npools, Teuchos::ParameterList& plist);

 public:
  Epetra_SerialDenseVector RespF;
  Epetra_SerialDenseVector TurnoverRates;
  Epetra_SerialDenseMatrix Tij;
  int nPools;

 protected:
  void InitCentury_(double percent_sand);
  void Init_(Teuchos::ParameterList& plist);

 private:
  SoilCarbonParameters(const SoilCarbonParameters& other); // not implemented
  SoilCarbonParameters& operator=(const SoilCarbonParameters& other); // not implemented
};


} // namespace
} // namespace

#endif
