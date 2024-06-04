/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!  Downregulates evaporation through a dessicated zone via soil resistance,
using the Sakagucki-Zeng model for soil resistance.

*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "Key.hh"
#include "StateDefs.hh"
#include "State.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

template <class cView_type, class View_type>
class EvaporationDownregulationSoilModel {
 public:
  static const int n_dependencies = 2;

  // this will mostly be differentiated with respect to pressure for flow
  // Jacobians, but none of the terms that _really_ depend on p are actually
  // implemented.  That would require differentiating RSoil with respect to
  // s_l, s_g, etc.  But only derivatives wrt potential evaporation are
  // implemented.  That will rarely if ever be p-dependent.  Therefore, this
  // is just turned off to avoid lengthy calculations with 0.
  static const bool provides_derivatives = false;
  static const std::string eval_type;

  explicit EvaporationDownregulationSoilModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Tag tag(plist->get<std::string>("tag"));
    my_key_ = { Keys::cleanPListName(*plist), tag };
    Key domain = Keys::getDomain(my_key_.first);

    // dependencies
    pet_key_ =
      Keys::readKeyTag(*plist, domain, "potential evaporation", "potential_evaporation", tag);
    rsoil_key_ = Keys::readKeyTag(*plist, domain, "soil resistance", "soil_resistance", tag);
  }


  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    res_ = res[0];

    pet_ = deps[0];
    rsoil_ = deps[1];
  }

  void freeViews()
  {
    res_ = View_type();
    pet_ = cView_type();
    rsoil_ = cView_type();
  }

  KeyTagVector getMyKeys() const
  {
    return {
      my_key_,
    };
  }
  KeyTagVector getDependencies() const { return { pet_key_, rsoil_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    res_(i, 0) = pet_(i, 0) / (1. + rsoil_(i, 0));
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const { assert(false); }

 private:
  View_type res_;
  cView_type pet_, rsoil_;

  KeyTag my_key_;
  KeyTag pet_key_, rsoil_key_;
};


template <class cView_type, class View_type>
const std::string EvaporationDownregulationSoilModel<cView_type, View_type>::eval_type =
  "evaporation downregulation, soil resistance";

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
