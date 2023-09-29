/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_

#include "Mesh_Algorithms.hh"
#include "CompositeVector.hh"

namespace Amanzi {

void
CopySurfaceToSubsurface(const CompositeVector& surf, CompositeVector& sub);

void
CopySubsurfaceToSurface(const CompositeVector& sub, CompositeVector& surf);

void
MergeSubsurfaceAndSurfacePressure(const CompositeVector& kr_surf,
                                  CompositeVector& sub_p,
                                  CompositeVector& surf_p);

//
// The next two helpers are getters and setters that can be used to get/set a
// component vector that is either FACE or BOUNDARY_FACE, selected at COMPILE
// time.  This is much better than the alternative of putting the check for
// which component one has in the innermost loop, then calling ViewComponent in
// the innermost loop.  Instead, any function can be templated with the type,
// chosen at compile-time, then used to construct the getter/setter and call
// the right templated method.  Then a check at runtime can be made at a much
// higher level (before the loop).  See for instance the implementations of
// CopySurfaceToSubsurface, or mpc_coupled_water, for example usages.
//
struct DomainFaceGetter {
  DomainFaceGetter(const AmanziMesh::Mesh& domain_mesh, const Epetra_MultiVector& vec)
    : domain_mesh_(domain_mesh), vec_(vec)
  {}

  template <AmanziMesh::Entity_ID FaceEntity>
  double get(const AmanziMesh::Entity_ID& f)
  {
    return vec_[0][f];
  }

 private:
  const AmanziMesh::Mesh& domain_mesh_;
  const Epetra_MultiVector& vec_;
};


template <>
inline double
DomainFaceGetter::get<AmanziMesh::BOUNDARY_FACE>(const AmanziMesh::Entity_ID& f)
{
  AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(domain_mesh_, f);
  return vec_[0][bf];
}


struct DomainFaceSetter {
  DomainFaceSetter(const AmanziMesh::Mesh& domain_mesh, Epetra_MultiVector& vec)
    : domain_mesh_(domain_mesh), vec_(vec)
  {}

  template <AmanziMesh::Entity_ID FaceEntity>
  double get(const AmanziMesh::Entity_ID& f)
  {
    return vec_[0][f];
  }

  template <AmanziMesh::Entity_ID FaceEntity>
  void set(const AmanziMesh::Entity_ID& f, const double& val)
  {
    vec_[0][f] = val;
  }

 private:
  const AmanziMesh::Mesh& domain_mesh_;
  Epetra_MultiVector& vec_;
};

template <>
inline double
DomainFaceSetter::get<AmanziMesh::BOUNDARY_FACE>(const AmanziMesh::Entity_ID& f)
{
  AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(domain_mesh_, f);
  return vec_[0][bf];
}

template <>
inline void
DomainFaceSetter::set<AmanziMesh::BOUNDARY_FACE>(const AmanziMesh::Entity_ID& f, const double& val)
{
  AmanziMesh::Entity_ID bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(domain_mesh_, f);
  vec_[0][bf] = val;
}


} // namespace Amanzi


#endif
