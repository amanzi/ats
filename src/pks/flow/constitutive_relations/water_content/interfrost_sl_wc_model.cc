/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The interfrost sl water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Interfrost water content portion sl.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "interfrost_sl_wc_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
InterfrostSlWcModel::InterfrostSlWcModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
InterfrostSlWcModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{}


// main method
double
InterfrostSlWcModel::WaterContent(double phi, double sl, double nl, double ni, double cv) const
{
  return cv * phi * sl * (-ni + nl);
}

double
InterfrostSlWcModel::DWaterContentDPorosity(double phi, double sl, double nl, double ni, double cv)
  const
{
  return cv * sl * (-ni + nl);
}

double
InterfrostSlWcModel::DWaterContentDSaturationLiquid(double phi,
                                                    double sl,
                                                    double nl,
                                                    double ni,
                                                    double cv) const
{
  return cv * phi * (-ni + nl);
}

double
InterfrostSlWcModel::DWaterContentDMolarDensityLiquid(double phi,
                                                      double sl,
                                                      double nl,
                                                      double ni,
                                                      double cv) const
{
  return cv * phi * sl;
}

double
InterfrostSlWcModel::DWaterContentDMolarDensityIce(double phi,
                                                   double sl,
                                                   double nl,
                                                   double ni,
                                                   double cv) const
{
  return -cv * phi * sl;
}

double
InterfrostSlWcModel::DWaterContentDCellVolume(double phi,
                                              double sl,
                                              double nl,
                                              double ni,
                                              double cv) const
{
  return phi * sl * (-ni + nl);
}

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
