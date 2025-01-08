/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Chonggang Xu (cxu@lanl.gov)
*/

/*
Utility functions for Vegetation.

*/


#ifndef ATS_BGC_VEG_HH_
#define ATS_BGC_VEG_HH_

namespace Amanzi {
namespace BGC {

// Functor for calculating QSat.
class QSat {
 public:
  QSat();

  void
  operator()(double tleafk, double pressure, double* es, double* esdT, double* qs, double* qsdT);

 private:
  double a0, a1, a2, a3, a4, a5, a6, a7, a8;
  double b0, b1, b2, b3, b4, b5, b6, b7, b8;
  double c0, c1, c2, c3, c4, c5, c6, c7, c8;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8;
};

// Determine time, in minutes of a given day of the year at a given latitude.
double
DayLength(double lat, int doy);

// Limit the highest temp?
double
HighTLim(double tleaf);

//solve the quadratic equation
void
Quadratic(double a, double b, double c, double* r1, double* r2);

// This function calculate the net photosynthetic rate based on Farquhar
// model, with updated leaf temperature based on energy balances by seperately solve light and RUBISCO-limited carboxylations
void
Photosynthesis(double PARi,
               double LUE,
               double LER,
               double pressure,
               double windv,
               double tair,
               double vp_air,
               double CO2a,
               double mp,
               double Vcmax25,
               double* A,
               double* tleaf,
               double* Resp,
               double* ET);
//// This function calculate the net photosynthetic rate based on Farquhar
// model, with updated leaf temperature based on energy balances by jointly solving light and RUBISCO-limited carboxylations
void
Photosynthesis0(double PARi,
                double LUE,
                double LER,
                double pressure,
                double windv,
                double tair,
                double vp_air,
                double CO2a,
                double mp,
                double Vcmax25,
                double* A,
                double* tleaf,
                double* Resp,
                double* ET);


} // namespace BGC
} // namespace Amanzi

#endif
