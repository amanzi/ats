/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates incoming longwave radiation from vapor pressure and air temperature.
/*!

This computes incoming longwave, from the atmosphere, by using the Boltzmann
equation and an emissivity empirically related to vapor pressure.

.. _longwave_evaluator-spec:
.. admonition:: longwave_evaluator-spec

    DEPENDENCIES:

    * `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
    * `"vapor pressure air key`" ``[string]`` **DOMAIN-vapor_pressure_air**

*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "Key.hh"
#include "StateDefs.hh"
#include "State.hh"
#include "seb_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

template <class cView_type, class View_type>
class IncomingLongwaveRadiationModel {
 public:
  static const int n_dependencies = 2;
  static const bool provides_derivatives = false;
  static const std::string eval_type;

  explicit IncomingLongwaveRadiationModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Tag tag(plist->get<std::string>("tag"));
    my_key_ = { Keys::cleanPListName(*plist), tag };
    auto domain = Keys::getDomain(my_key_.first);

    // dependencies
    air_temp_key_ = Keys::readKeyTag(*plist, domain, "air temperature", "air_temperature", tag);
    vp_air_key_ = Keys::readKeyTag(*plist, domain, "vapor pressure air", "vapor_pressure_air", tag);

    Teuchos::ParameterList& model_list = plist->sublist("model parameters");
    scaling_ = model_list.get<double>("scaling [-]", 1);
  }


  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    res_ = res[0];
    air_temp_ = deps[0];
    vp_air_ = deps[1];
  }

  void freeViews()
  {
    res_ = View_type();
    air_temp_ = cView_type();
    vp_air_ = cView_type();
  }

  KeyTagVector getMyKeys() const { return { my_key_ }; }
  KeyTagVector getDependencies() const { return { air_temp_key_, vp_air_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    res_(i, 0) = scaling_ * Functions::atmosphereLongwaveRadiation(air_temp_(i, 0), vp_air_(i, 0));
  }

  // derivatives not currently provided
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const { assert(false); }

 private:
  KeyTag my_key_;
  KeyTag air_temp_key_, vp_air_key_;

  View_type res_;
  cView_type air_temp_, vp_air_;
  double scaling_;
};


template <class cView_type, class View_type>
const std::string IncomingLongwaveRadiationModel<cView_type, View_type>::eval_type =
  "incoming longwave radiation";

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
