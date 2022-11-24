/*
  The {evalNameString} model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
{docDict}
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_{ namespaceCaps } _{ evalNameCaps } _MODEL_HH_
#define AMANZI_ { namespaceCaps } _{ evalNameCaps } _MODEL_HH_

namespace Amanzi {
{
  namespace {namespace}
  {
    {
      namespace Relations {
      {
        class {
          evalClassName
        } Model{ {

          public : explicit { evalClassName } Model(Teuchos::ParameterList & plist);

        { modelMethodDeclaration }

        {
          modelDerivDeclarationList
        }

       protected:
        void InitializeFromPlist_(Teuchos::ParameterList & plist);

       protected:
        {
          paramDeclarationList
        }

      }  // namespace Relations
      }; // namespace Relations
    }
  } // namespace }
} // namespace Amanzi
} // namespace Amanzi
}
} //namespace

#endif
