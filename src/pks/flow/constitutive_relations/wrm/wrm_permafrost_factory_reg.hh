/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM implementations.
   ------------------------------------------------------------------------- */

#include "wrm_permafrost_factory.hh"

// explicity instantitate the static data of Factory<WRM>
template <>
Amanzi::Utils::Factory<Amanzi::Flow::WRMPermafrostModel>::map_type*
  Amanzi::Utils::Factory<Amanzi::Flow::WRMPermafrostModel>::map_;
