/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

// TPLs
#include "boost/algorithm/string.hpp"

#include "transport_ats.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Re-partition components between liquid and gas phases.
******************************************************************* */
void
Transport_ATS::PrepareAirWaterPartitioning_()
{
  henry_law_ = true;
  for (int i = 0; i < num_gaseous; i++) {
    int ig = num_aqueous + i;
    std::string name_l = boost::replace_all_copy(component_names_[ig], "(g)", "(l)");

    int il = FindComponentNumber(name_l);
    air_water_map_.push_back(il);

    if (il < 0 || il >= num_aqueous) {
      henry_law_ = false;
      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "Gas component \"" << component_names_[ig]
                   << "\" has no matching liquid component \"" << name_l << "\"\n";
      }
      break;
    }
  }

  if (henry_law_) {
    Teuchos::Array<double> empty;
    kH_ = plist_->sublist("molecular diffusion")
            .get<Teuchos::Array<double>>("air-water partitioning coefficient", empty)
            .toVector();
  } else {
    air_water_map_.clear();
  }
}


/* *******************************************************************
* Re-partition components between liquid and gas phases.
******************************************************************* */
void
Transport_ATS::MakeAirWaterPartitioning_()
{
  Epetra_MultiVector& tcc_c = *tcc_tmp->ViewComponent("cell", false);
  const Epetra_MultiVector& sat_l = *ws_;

  for (int i = 0; i < num_gaseous; ++i) {
    int ig = num_aqueous + i;
    int il = air_water_map_[i];

    for (int c = 0; c < ncells_owned; c++) {
      double sl = sat_l[0][c];
      double total = tcc_c[il][c] * sl + tcc_c[ig][c] * (1.0 - sl);
      tcc_c[ig][c] = total / (1.0 + (kH_[i] - 1.0) * sl);
      tcc_c[il][c] = tcc_c[ig][c] * kH_[i];
    }
  }
}

} // namespace Transport
} // namespace Amanzi
