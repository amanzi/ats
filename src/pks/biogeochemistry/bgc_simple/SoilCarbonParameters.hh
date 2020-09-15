/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Chonggang Xu (cxu@lanl.gov)
*/
/*!

Parameters for the Century model for biogeochemistry.

.. _soil-carbon-spec:
.. admonition:: soil-carbon-spec

    * `"percent sand`" ``[double]`` Soil makeup.  Ranges from [0,100]

*/

#ifndef ATS_BGC_SOIL_CARBON_PARAMS_HH_
#define ATS_BGC_SOIL_CARBON_PARAMS_HH_

#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

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
