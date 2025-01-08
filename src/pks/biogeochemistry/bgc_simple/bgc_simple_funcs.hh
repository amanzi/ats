/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Chonggang Xu (cxu@lanl.gov)
*/

/*
Main functions for biogeochemistry on a column.

*/


#ifndef ATS_BGC_SIMPLE_FUNCS_HH_
#define ATS_BGC_SIMPLE_FUNCS_HH_

#include <vector>

#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"

#include "utils.hh"
#include "PFT.hh"
#include "SoilCarbon.hh"

namespace Amanzi {
namespace BGC {

void
BGCAdvance(double t,
           double dt,
           double gridarea,
           double cryoturbation_coef,
           const MetData& met,
           const Epetra_SerialDenseVector& SoilTArr,
           const Epetra_SerialDenseVector& SoilArr,
           const Epetra_SerialDenseVector& SoilDArr,
           const Epetra_SerialDenseVector& SoilThicknessArr,
           std::vector<Teuchos::RCP<PFT>>& pftarr,
           std::vector<Teuchos::RCP<SoilCarbon>>& soilcarr,
           Epetra_SerialDenseVector& SoilCO2Arr,
           Epetra_SerialDenseVector& TransArr,
           double& sw_shaded);

void
Cryoturbate(double dt,
            const Epetra_SerialDenseVector& SoilTArr,
            const Epetra_SerialDenseVector& SoilDArr,
            const Epetra_SerialDenseVector& SoilThicknessArr,
            std::vector<Teuchos::RCP<SoilCarbon>>& soilcarr,
            double diffusion_coef);

void
Cryoturbate(double dt,
            const Epetra_SerialDenseVector& SoilTArr,
            const Epetra_SerialDenseVector& SoilDArr,
            const Epetra_SerialDenseVector& SoilThicknessArr,
            std::vector<Teuchos::RCP<SoilCarbon>>& soilcarr,
            std::vector<double>& diffusion_coefs);


} // namespace BGC
} // namespace Amanzi

#endif
