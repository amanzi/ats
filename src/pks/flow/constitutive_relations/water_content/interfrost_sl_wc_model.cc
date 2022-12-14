/*
  The interfrost sl water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Interfrost water content portion sl.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "interfrost_sl_wc_model.hh"

namespace Amanzi {
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
} // namespace Amanzi
