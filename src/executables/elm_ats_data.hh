/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* Saved ELM-ATS interface data

 Data shared/referenced between ELM and ATS.


 Authors: F.-M. Yuan (yuanf@ornl.gov), Joe Beisman(jjbeisman@ornl.gov) Ethan Coon (coonet@ornl.gov)
 Year: 2022
*/

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

#ifndef ELM_ATS_DATA_H
#define ELM_ATS_DATA_H

// saved data from ELM
typedef struct {
    int NX_=0;
    int NY_=0;
    int NZ_=0;
    int NCOLS_=0;
    int NCELLS_=0;

    std::vector<double> surf_gridsX_;     // elm surface grids X coord in meters converted from lon, size of at least 2 (for 1 grid)
    std::vector<double> surf_gridsY_;     // elm surface grids Y coord in meters converted from lat, size of at least 2 (for 1 grid)
    std::vector<double> surf_gridsZ_;     // elm surface grids center elevation in meters, size of at least 1 (for 1 grid)
    std::vector<double> cols_verticesZ_;  // elm soil column nodes in meters, size of 16 (for 15 layers) from top to bottom with upward positive;

    std::vector<double> porosity_;        // unit: [-]
    std::vector<double> hksat_;           // unit: [mm s-1]
    std::vector<double> CH_bsw_;          // unit: [-]
    std::vector<double> CH_smpsat_;       // unit: [-]
    std::vector<double> CH_sr_;           // unit: [-]
    std::vector<double> eff_porosity_;    // unit: [-]
    std::vector<double> zwt_;             // unit: [m]

} elm_data;

#endif
