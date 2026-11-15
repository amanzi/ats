#ifndef ATS_ECOSIM_SURFACE_RESISTANCE_HH_
#define ATS_ECOSIM_SURFACE_RESISTANCE_HH_
#include "State.hh"
#include "Mesh.hh"
#include <cmath>
namespace Amanzi { namespace Energy {
// Apply the frozen EcoSIM layer-zero series resistance to the existing soil
// half-cell resistance. The same transformation is used by residual and PC.
inline void EcoSurfaceResistance(const State& S,const std::string& domain,const Tag& tag,
 const CompositeVector& conductivity,CompositeVector& output,bool derivative=false){
 auto sd=(domain.empty() || domain=="domain")?std::string("surface"):std::string("surface_")+domain;
 auto key=Keys::getKey(sd,"ecosim_layer0_thermal_resistance");
 if(!S.HasRecord(key,tag) || !output.HasComponent("face"))return;
 const auto surface=S.GetMesh(sd);const auto mesh=conductivity.Mesh();
 const auto& resistance=*S.Get<CompositeVector>(key,tag).ViewComponent("cell",false);
 auto fraction_key=Keys::getKey(sd,"ecosim_layer0_fraction");
 const Epetra_MultiVector* fraction=S.HasRecord(fraction_key,tag)?
   S.Get<CompositeVector>(fraction_key,tag).ViewComponent("cell",false).get():nullptr;
 const auto& k=*conductivity.ViewComponent("face",false);auto& result=*output.ViewComponent("face",false);
 for(int sc=0;sc<resistance.MyLength();++sc){
  int f=surface->getEntityParent(AmanziMesh::CELL,sc);auto cells=mesh->getFaceCells(f);
  double distance=AmanziGeometry::norm(mesh->getFaceCentroid(f)-mesh->getCellCentroid(cells[0]));
  double r=resistance[0][sc];
  double covered=fraction?(*fraction)[0][sc]:1.;
  // With no pond/litter store the surface temperature is an algebraic trace:
  // retain its dry Tsurface=Tsoil closure, which carries zero net heat.
  double scale=covered>0.?covered:1.;
  if(distance>0.){
   if(derivative){double factor=std::max(0.,1.-k[0][f]*r/(scale*distance));result[0][f]*=scale*factor*factor;}
   else result[0][f]=scale*k[0][f]/(1.+k[0][f]*r/distance);
  }
 }
}
} }
#endif
