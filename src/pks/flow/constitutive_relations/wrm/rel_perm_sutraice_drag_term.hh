/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#pragma once
#include <cmath>

namespace Amanzi {
namespace ATS_Physics {            
namespace Flow {
namespace SutraIceTerm {

//inline double
//dragcoef(double sl, double si, double sr, double omega)
//{
//  return std::pow(10, -omega * si / (sl + si - sr + 1e-6));
//}
//
//inline double
//d_dragcoef_dsl(double sl, double si, double sr, double omega)
//{
//  double coef = std::pow(10, -omega * si / (sl + si - sr + 1e-6));
//  return coef * std::log(10) * omega *si * std::pow(sl + si -sr + 1e-6, -2);
//}
//
//inline double
//d_dragcoef_dsi(double sl, double si, double sr, double omega)
//{
//  double coef = std::pow(10, -omega * si / (sl + si - sr + 1e-6));
//  return coef * std::log(10) * omega * (sr - sl) * std::pow(sl + si - sr + 1e-6, -2);
//}


inline double
dragcoef(double sl, double sg, double sr, double omega)
{
  return std::pow(10, -omega * (1 - sl - sg) / (1 - sg - sr + 1e-6));
}

inline double
d_dragcoef_dsl(double sl, double sg, double sr, double omega)
{
  double coef = std::pow(10, -omega * (1 - sl - sg) / (1 - sg - sr + 1e-6));
  return coef * std::log(10) * omega / (1 - sg - sr + 1e-6);
}

inline double
d_dragcoef_dsg(double sl, double sg, double sr, double omega)
{
  double coef = std::pow(10, -omega * (1 - sl - sg) / (1 - sg - sr + 1e-6));
  return coef * std::log(10) * omega * (sl - sr) * std::pow(1 - sg - sr + 1e-6, -2);
}


} // namespace SutraIceTerm
} // namespace Flow
}
} // namespace Amanzi
