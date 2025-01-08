/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Chonggang Xu (cxu@lanl.gov)
*/

/*

Soil carbon parameters data structures class

*/

#include "dbc.hh"
#include "SoilCarbonParameters.hh"

namespace Amanzi {
namespace BGC {

SoilCarbonParameters::SoilCarbonParameters(int nPools_, double percent_sand)
  : nPools(nPools_), RespF(nPools_), TurnoverRates(nPools_), Tij(nPools_, nPools_)
{
  InitCentury_(percent_sand);
}


SoilCarbonParameters::SoilCarbonParameters(int nPools_, Teuchos::ParameterList& plist)
  : nPools(nPools_), RespF(nPools_), TurnoverRates(nPools_), Tij(nPools_, nPools_)
{
  Init_(plist);
}


void
SoilCarbonParameters::InitCentury_(double percent_sand)
{
  double tt = 0.85 - 0.68 * 0.01 * (100 - percent_sand);

  // initialize the respiration fraction
  RespF[0] = 0.0;
  RespF[1] = 0.55;
  RespF[2] = 0.5;
  RespF[3] = 0.5;
  RespF[4] = tt;
  RespF[5] = 0.55;
  RespF[6] = 0.55;

  // initialize baseline turnover rates, years
  TurnoverRates[0] = 4.1;
  TurnoverRates[1] = 0.06;
  TurnoverRates[2] = 0.25;
  TurnoverRates[3] = 0.25;
  TurnoverRates[4] = 0.17;
  TurnoverRates[5] = 6.1;
  TurnoverRates[6] = 270.0;

  // initialize conversion factors
  Tij[0][2] = 0.76;
  Tij[0][3] = 0.24;
  Tij[1][4] = 1.0;
  Tij[2][4] = 1.0;
  Tij[3][5] = 1.0;
  Tij[4][5] = 1.0 - 0.004 / (1.0 - tt);
  Tij[4][6] = 1.0 - Tij[4][5];
  Tij[5][4] = 0.93;
  Tij[5][6] = 0.07;
  Tij[6][4] = 1.0;
}

void
SoilCarbonParameters::Init_(Teuchos::ParameterList& plist)
{
  std::string model = plist.get<std::string>("soil carbon model", "century");
  if (model == "century") {
    double p_sand = plist.get<double>("percent sand");
    InitCentury_(p_sand);
  } else {
    AMANZI_ASSERT(0);
  }

  // special case init from plist here...
}

} // namespace BGC
} // namespace Amanzi
