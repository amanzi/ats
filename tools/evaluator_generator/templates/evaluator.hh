/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The {evalNameString} evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
{docDict}

*/

#ifndef AMANZI_{ namespaceCaps } _{ evalNameCaps } _EVALUATOR_HH_
#define AMANZI_ { namespaceCaps } _{ evalNameCaps } _EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
{
  namespace {namespace}
  {
    {
      namespace Relations {
      {
        class {
          evalClassName
        } Model;

        class {
          evalClassName
        } Evaluator : public EvaluatorSecondaryMonotypeCV{ {

          public : explicit { evalClassName } Evaluator(Teuchos::ParameterList & plist);
        { evalClassName } Evaluator(const { evalClassName } Evaluator & other);

        virtual Teuchos::RCP<Evaluator> Clone() const;

        // Required methods from EvaluatorSecondaryMonotypeCV
        virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                                    const Teuchos::Ptr<CompositeVector>& result);
        virtual void EvaluateFieldPartialDerivative_(
          const Teuchos::Ptr<State>& S, Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

        Teuchos::RCP<{ evalClassName } Model> get_model()
        {
          {
            return model_;
          }
        }

       protected:
        void InitializeFromPlist_();

        { keyDeclarationList }

        Teuchos::RCP<{ evalClassName } Model>
          model_;

       private:
        static Utils::RegisteredFactory<Evaluator, { evalClassName } Evaluator> reg_;

      } // namespace Relations
      }; // namespace Relations
    }
  } // namespace }
} // namespace Amanzi
} // namespace Amanzi
}
} //namespace

#endif
