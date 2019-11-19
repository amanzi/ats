/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Anisotropic mechanical dispersion model.
*/

#ifndef AMANZI_MDM_BURNETT_FRIND_HH_
#define AMANZI_MDM_BURNETT_FRIND_HH_

// TPLs
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "factory.hh"
#include "Point.hh"
#include "Tensor.hh"

// Transport
#include "MDM.hh"

namespace Amanzi {
namespace Transport {

class MDM_BurnettFrind : public MDM {
 public:
  explicit MDM_BurnettFrind(Teuchos::ParameterList& plist);
  ~MDM_BurnettFrind() {};
  
  // Required methods from the base class
  // -- dispersion tensor of rank 2
  WhetStone::Tensor mech_dispersion(
      const AmanziGeometry::Point& u, int axi_symmetry, double s, double phi) const;

  // -- the model is valid if at least one parameter is not zero.
  bool is_valid() const { return (alphaL_ + alphaTH_ + alphaTV_ != 0.0); }

  // -- check model applicability
  void set_dim(int dim) { AMANZI_ASSERT(dim == 3); dim_ = dim; }

 private:
  double alphaL_, alphaTH_, alphaTV_;

  static Utils::RegisteredFactory<MDM, MDM_BurnettFrind> factory_;
};

}  // namespace Transport
}  // namespace Amanzi
 
#endif
