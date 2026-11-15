#ifndef ATS_ECOSIM_SOURCE_ROUTING_HH_
#define ATS_ECOSIM_SOURCE_ROUTING_HH_
#include "State.hh"
#include "Teuchos_Array.hpp"
namespace Amanzi {
// EcoSIM declares this metadata during parse; physical PKs consume it during
// Setup, after all PKs have parsed, independent of their construction order.
inline Key AttachEcoSIMCellSource(State& S,const Key& primary,const Tag& tag,
 const std::string& metadata,Key source,bool& enabled){
 if(!S.HasEvaluatorList(primary))return source;
 const auto& spec=S.GetEvaluatorList(primary);
 if(!spec.isParameter(metadata))return source;
 auto extra=spec.get<std::string>(metadata);
 if(!enabled){enabled=true;return extra;}
 if(source==extra)return source;
 S.Require<CompositeVector,CompositeVectorSpace>(source,tag)
   .SetMesh(S.GetMesh(Keys::getDomain(primary)))->AddComponent("cell",AmanziMesh::CELL,1);
 S.RequireEvaluator(source,tag);
 if(S.GetEvaluator(source,tag).IsDependency(S,extra,tag))return source;
 // A user sum may already contain the raw EcoSIM source: add only its
 // temperature correction in that case, never the same reference flux twice.
 const auto& added=S.GetEvaluatorList(extra);
 if(added.isParameter("raw coupled source key")){
   auto raw=added.get<std::string>("raw coupled source key");
   if(source==raw || S.GetEvaluator(source,tag).IsDependency(S,raw,tag))
     extra=added.get<std::string>("coupled temperature correction key");
 }
 auto combined=source+"_with_ecosim_ground";
 auto& sum=S.GetEvaluatorList(combined);
 sum.set("evaluator type","additive evaluator");
 sum.set("dependencies",Teuchos::Array<std::string>{source,extra});
 return combined;
}
}
#endif
