/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The richards water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "richards_water_content_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
RichardsWaterContentModel::RichardsWaterContentModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
RichardsWaterContentModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{}


// main method
double
RichardsWaterContentModel::WaterContent(double phi, double sl, double nl, double cv) const
{
  return cv * nl * phi * sl;
}

double
RichardsWaterContentModel::DWaterContentDPorosity(double phi, double sl, double nl, double cv) const
{
  return cv * nl * sl;
}

double
RichardsWaterContentModel::DWaterContentDSaturationLiquid(double phi,
                                                          double sl,
                                                          double nl,
                                                          double cv) const
{
  return cv * nl * phi;
}

double
RichardsWaterContentModel::DWaterContentDMolarDensityLiquid(double phi,
                                                            double sl,
                                                            double nl,
                                                            double cv) const
{
  return cv * phi * sl;
}

double
RichardsWaterContentModel::DWaterContentDCellVolume(double phi,
                                                    double sl,
                                                    double nl,
                                                    double cv) const
{
  return nl * phi * sl;
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
