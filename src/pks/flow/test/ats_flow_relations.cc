/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

#include "UnitTest++.h"
#include "AmanziTypes.hh"
#include "ats_flow_relations_registration.hh"

using namespace Amanzi;
using View_type = Kokkos::View<double**, DefaultHostMemorySpace>;
using cView_type = Kokkos::View<const double**, DefaultHostMemorySpace>;


SUITE(FLOW_RELATIONS)
{
  TEST(RICHARDS_WATER_CONTENT)
  {
    auto plist = Teuchos::rcp(new Teuchos::ParameterList("saturation"));
    Flow::Relations::RichardsWaterContentModel<cView_type, View_type> model(plist);
    View_type sat("saturation", 4, 1);
    sat(0, 0) = 0.;
    sat(1, 0) = .25;
    sat(2, 0) = .5;
    sat(3, 0) = 1.;

    View_type poro("porosity", 4, 1);
    Kokkos::deep_copy(poro, 0.5);
    View_type cv("cell volume", 4, 1);
    Kokkos::deep_copy(cv, 1);
    View_type dens("density", 4, 1);
    Kokkos::deep_copy(dens, 1);

    View_type wc("water_content", 4, 1);
    State s;
    model.setViews({ dens, sat, poro, cv }, { wc }, s);
    for (int i = 0; i != 4; ++i) model(i);
    CHECK_CLOSE(0.0, wc(0, 0), 1.0e-10);
    CHECK_CLOSE(0.125, wc(1, 0), 1.0e-10);
    CHECK_CLOSE(0.25, wc(2, 0), 1.0e-10);
    CHECK_CLOSE(0.5, wc(3, 0), 1.0e-10);
  }
}
