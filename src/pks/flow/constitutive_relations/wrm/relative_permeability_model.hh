/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Uses WRM class instances to compute the scalar portion of the hydraulic conductance.
/*!

This computes the product,

.. math::
   k = \frac{n}{\mu} k_r

which is the scalar portion of the hydraulic conductance (excludes the absolute
permeability).  Relative permeability, k_r, is computed using a water retention
model, which provides generic methods for this quantity based on internal
parameters (e.g. van Genuchten/Mualem, Brooks & Corey, etc.).

While we call this "relative permeability", it is actually the product above,
and might better be called the scalar hydraulic conductance?

type : `"relative permeability`"

.. _relative-permeability-model-spec
.. admonition:: relative-permeability-model-spec

   * `"model parameters`" ``[wrm-spec-list]``

   KEYS
   - `"saturation`" **DOMAIN-saturation_liquid**
   - `"density`" **DOMAIN-molar_density_liquid**
   - `"viscosity`" **DOMAIN-viscosity**  
                
*/

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "StateDefs.hh"
#include "wrm_van_genuchten.hh"

namespace Amanzi {

class State;
class RelativePermeabilityEvaluator;

namespace Flow {
namespace Relations {

template <class cView_type, class View_type, class WRM_type>
class RelativePermeabilityModel {
 public:
  static const int n_dependencies = 1;
  static const bool provides_derivatives = true;
  static const std::string eval_type;

  RelativePermeabilityModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : model_(plist->sublist("model parameters"))
  {
    my_key_ = { Keys::cleanPListName(*plist), Tag{ plist->get<std::string>("tag") } };
    auto domain = Keys::getDomain(my_key_.first);
    s_key_ = Keys::readKeyTag(*plist, domain, "saturation", "saturation_liquid", my_key_.second);
  }

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    res_ = res[0];
    s_ = deps[0];
  }

  void freeViews()
  {
    res_ = View_type();
    s_ = cView_type();
  }

  KeyTagVector getMyKeys() const
  {
    return {
      my_key_,
    };
  }
  KeyTagVector getDependencies() const
  {
    return {
      s_key_,
    };
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    res_(i, 0) = model_.k_relative(s_(i, 0));
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    res_(i, 0) = model_.d_k_relative(s_(i, 0));
  }

 private:
  View_type res_;
  cView_type s_;
  KeyTag my_key_, s_key_;

  WRM_type model_;
};


template <class cView_type, class View_type, class WRM_type>
const std::string RelativePermeabilityModel<cView_type, View_type, WRM_type>::eval_type =
  "relative permeability " + WRM_type::eval_type;


template <class cView_type, class View_type>
using RelativePermeabilityVanGenuchtenModel =
  RelativePermeabilityModel<cView_type, View_type, WRMVanGenuchten>;


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
