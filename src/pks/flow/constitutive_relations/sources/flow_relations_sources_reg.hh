/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#pragma once

#include "Factory.hh"
#include "registration_macro.hh"

#include "soil_resistance_sakagucki_zeng_evaluator.hh"

namespace Amanzi {

template <>
REGISTER(Flow::Relations::SoilResistanceSakaguckiZengEvaluator);
template <>
REGISTER(Flow::Relations::SoilResistanceSakaguckiZengEvaluatorByMaterial);

} // namespace Amanzi
