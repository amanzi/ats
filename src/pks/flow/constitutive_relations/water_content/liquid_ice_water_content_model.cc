/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The liquid + ice water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Water content for a two-phase, liquid+ice evaluator.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "liquid_ice_water_content_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
LiquidIceWaterContentModel::LiquidIceWaterContentModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
LiquidIceWaterContentModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{}


// main method
double
LiquidIceWaterContentModel::WaterContent(double phi,
                                         double sl,
                                         double nl,
                                         double si,
                                         double ni,
                                         double cv) const
{
  return cv * phi * (ni * si + nl * sl);
}

double
LiquidIceWaterContentModel::DWaterContentDPorosity(double phi,
                                                   double sl,
                                                   double nl,
                                                   double si,
                                                   double ni,
                                                   double cv) const
{
  return cv * (ni * si + nl * sl);
}

double
LiquidIceWaterContentModel::DWaterContentDSaturationLiquid(double phi,
                                                           double sl,
                                                           double nl,
                                                           double si,
                                                           double ni,
                                                           double cv) const
{
  return cv * nl * phi;
}

double
LiquidIceWaterContentModel::DWaterContentDMolarDensityLiquid(double phi,
                                                             double sl,
                                                             double nl,
                                                             double si,
                                                             double ni,
                                                             double cv) const
{
  return cv * phi * sl;
}

double
LiquidIceWaterContentModel::DWaterContentDSaturationIce(double phi,
                                                        double sl,
                                                        double nl,
                                                        double si,
                                                        double ni,
                                                        double cv) const
{
  return cv * ni * phi;
}

double
LiquidIceWaterContentModel::DWaterContentDMolarDensityIce(double phi,
                                                          double sl,
                                                          double nl,
                                                          double si,
                                                          double ni,
                                                          double cv) const
{
  return cv * phi * si;
}

double
LiquidIceWaterContentModel::DWaterContentDCellVolume(double phi,
                                                     double sl,
                                                     double nl,
                                                     double si,
                                                     double ni,
                                                     double cv) const
{
  return phi * (ni * si + nl * sl);
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
