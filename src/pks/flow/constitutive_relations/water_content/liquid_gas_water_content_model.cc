/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The liquid + gas water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Water content for a two-phase, liquid+water vapor evaluator.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "liquid_gas_water_content_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
LiquidGasWaterContentModel::LiquidGasWaterContentModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
LiquidGasWaterContentModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{}


// main method
double
LiquidGasWaterContentModel::WaterContent(double phi,
                                         double sl,
                                         double nl,
                                         double sg,
                                         double ng,
                                         double omega,
                                         double cv) const
{
  return cv * phi * (ng * omega * sg + nl * sl);
}

double
LiquidGasWaterContentModel::DWaterContentDPorosity(double phi,
                                                   double sl,
                                                   double nl,
                                                   double sg,
                                                   double ng,
                                                   double omega,
                                                   double cv) const
{
  return cv * (ng * omega * sg + nl * sl);
}

double
LiquidGasWaterContentModel::DWaterContentDSaturationLiquid(double phi,
                                                           double sl,
                                                           double nl,
                                                           double sg,
                                                           double ng,
                                                           double omega,
                                                           double cv) const
{
  return cv * nl * phi;
}

double
LiquidGasWaterContentModel::DWaterContentDMolarDensityLiquid(double phi,
                                                             double sl,
                                                             double nl,
                                                             double sg,
                                                             double ng,
                                                             double omega,
                                                             double cv) const
{
  return cv * phi * sl;
}

double
LiquidGasWaterContentModel::DWaterContentDSaturationGas(double phi,
                                                        double sl,
                                                        double nl,
                                                        double sg,
                                                        double ng,
                                                        double omega,
                                                        double cv) const
{
  return cv * ng * omega * phi;
}

double
LiquidGasWaterContentModel::DWaterContentDMolarDensityGas(double phi,
                                                          double sl,
                                                          double nl,
                                                          double sg,
                                                          double ng,
                                                          double omega,
                                                          double cv) const
{
  return cv * omega * phi * sg;
}

double
LiquidGasWaterContentModel::DWaterContentDMolFracGas(double phi,
                                                     double sl,
                                                     double nl,
                                                     double sg,
                                                     double ng,
                                                     double omega,
                                                     double cv) const
{
  return cv * ng * phi * sg;
}

double
LiquidGasWaterContentModel::DWaterContentDCellVolume(double phi,
                                                     double sl,
                                                     double nl,
                                                     double sg,
                                                     double ng,
                                                     double omega,
                                                     double cv) const
{
  return phi * (ng * omega * sg + nl * sl);
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
