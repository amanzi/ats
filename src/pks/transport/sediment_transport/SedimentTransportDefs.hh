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

#ifndef AMANZI_TRANSPORT_CONSTANTS_HH_
#define AMANZI_TRANSPORT_CONSTANTS_HH_

namespace Amanzi {
namespace Transport {

// const int TRANSPORT_PHASE_LIQUID = 0;  // phases from 0 to 1
// const int TRANSPORT_PHASE_GAS = 1;
// const int TRANSPORT_NUMBER_PHASES = 2;

const double TRANSPORT_LARGE_TIME_STEP = 1e+99;
const double TRANSPORT_SMALL_TIME_STEP = 1e-12;

const int TRANSPORT_BC_CONSTANT_TCC = 1;
const int TRANSPORT_BC_DISPERSION_FLUX = 2;
const int TRANSPORT_BC_NULL = 3;

const int TRANSPORT_FLOW_STEADYSTATE = 1;
const int TRANSPORT_FLOW_TRANSIENT = 2;

const double TRANSPORT_CONCENTRATION_OVERSHOOT = 1e-6;
const double TRANSPORT_CONCENTRATION_INFINITY = 1e+99;

const int TRANSPORT_HEX_FACES = 6; // Hexahedron is the common element
const int TRANSPORT_HEX_NODES = 8;
const int TRANSPORT_HEX_EDGES = 12;

const int TRANSPORT_QUAD_FACES = 4; // Quadrilateral is the common element
const int TRANSPORT_QUAD_NODES = 4;
const int TRANSPORT_QUAD_EDGES = 4;

const int TRANSPORT_MAX_FACES = 14; // Kelvin's tetrakaidecahedron
const int TRANSPORT_MAX_NODES = 47; // These polyhedron parameters must
const int TRANSPORT_MAX_EDGES = 60; // be calculated in Init().

const int TRANSPORT_DISPERSION_METHOD_TPFA = 1;
const int TRANSPORT_DISPERSION_METHOD_NLFV = 2;

const int TRANSPORT_LIMITER_BARTH_JESPERSEN = 1;
const int TRANSPORT_LIMITER_TENSORIAL = 2;
const int TRANSPORT_LIMITER_KUZMIN = 3;
const double TRANSPORT_LIMITER_TOLERANCE = 1e-14;

const int TRANSPORT_INTERNAL_ERROR = 911; // contact (lipnikov@lanl.gov)

} // namespace Transport
} // namespace Amanzi

#endif
