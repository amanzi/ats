/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The three phase water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Water content for a three-phase, gas+liquid+ice evaluator.

*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "three_phase_water_content_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
ThreePhaseWaterContentModel::ThreePhaseWaterContentModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
ThreePhaseWaterContentModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{}


// main method
double
ThreePhaseWaterContentModel::WaterContent(double phi,
                                          double sl,
                                          double nl,
                                          double si,
                                          double ni,
                                          double sg,
                                          double ng,
                                          double omega,
                                          double cv) const
{
  return cv * phi * (ng * omega * sg + ni * si + nl * sl);
}

double
ThreePhaseWaterContentModel::DWaterContentDPorosity(double phi,
                                                    double sl,
                                                    double nl,
                                                    double si,
                                                    double ni,
                                                    double sg,
                                                    double ng,
                                                    double omega,
                                                    double cv) const
{
  return cv * (ng * omega * sg + ni * si + nl * sl);
}

double
ThreePhaseWaterContentModel::DWaterContentDSaturationLiquid(double phi,
                                                            double sl,
                                                            double nl,
                                                            double si,
                                                            double ni,
                                                            double sg,
                                                            double ng,
                                                            double omega,
                                                            double cv) const
{
  return cv * nl * phi;
}

double
ThreePhaseWaterContentModel::DWaterContentDMolarDensityLiquid(double phi,
                                                              double sl,
                                                              double nl,
                                                              double si,
                                                              double ni,
                                                              double sg,
                                                              double ng,
                                                              double omega,
                                                              double cv) const
{
  return cv * phi * sl;
}

double
ThreePhaseWaterContentModel::DWaterContentDSaturationIce(double phi,
                                                         double sl,
                                                         double nl,
                                                         double si,
                                                         double ni,
                                                         double sg,
                                                         double ng,
                                                         double omega,
                                                         double cv) const
{
  return cv * ni * phi;
}

double
ThreePhaseWaterContentModel::DWaterContentDMolarDensityIce(double phi,
                                                           double sl,
                                                           double nl,
                                                           double si,
                                                           double ni,
                                                           double sg,
                                                           double ng,
                                                           double omega,
                                                           double cv) const
{
  return cv * phi * si;
}

double
ThreePhaseWaterContentModel::DWaterContentDSaturationGas(double phi,
                                                         double sl,
                                                         double nl,
                                                         double si,
                                                         double ni,
                                                         double sg,
                                                         double ng,
                                                         double omega,
                                                         double cv) const
{
  return cv * ng * omega * phi;
}

double
ThreePhaseWaterContentModel::DWaterContentDMolarDensityGas(double phi,
                                                           double sl,
                                                           double nl,
                                                           double si,
                                                           double ni,
                                                           double sg,
                                                           double ng,
                                                           double omega,
                                                           double cv) const
{
  return cv * omega * phi * sg;
}

double
ThreePhaseWaterContentModel::DWaterContentDMolFracGas(double phi,
                                                      double sl,
                                                      double nl,
                                                      double si,
                                                      double ni,
                                                      double sg,
                                                      double ng,
                                                      double omega,
                                                      double cv) const
{
  return cv * ng * phi * sg;
}

double
ThreePhaseWaterContentModel::DWaterContentDCellVolume(double phi,
                                                      double sl,
                                                      double nl,
                                                      double si,
                                                      double ni,
                                                      double sg,
                                                      double ng,
                                                      double omega,
                                                      double cv) const
{
  return phi * (ng * omega * sg + ni * si + nl * sl);
}

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
