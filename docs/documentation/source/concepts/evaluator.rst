Evaluators and the Dependency Graph
===================================

Previously, all concepts were used to describe data -- now we will
describe the concepts used to describe how we compute data.

**Evaluators** are the finest-grained unit of computation done in
Amanzi-ATS. Typically, they are used to define a term in an equation,
calculating one or more data values.  Evaluators are used for nearly
all persistent data; those data and are split into four categories:

   * **primary variables** are used to refer to the variables
     solved for by solving a PDE.  For instance, the primary variable of
     an energy equation is typically a temperature; an energy PDE would
     create a primary variable evaluator to tell the code that it intends
     to compute this value.

   * **independent variables** are functions of space and time, but
     nothing else.  Independent variable evaluators compute functions, or
     load user-provided data from files.
     
   * **secondary variables** are functions of other variables.
     Secondary variable evaluators implement models to compute
     secondary variables, e.g. water density as a function of
     temperature, or saturation as a function of pressure.

   * **aliased variables** are simply a way of creating multiple
     logical references to the same data.  This is shearly for
     performance concerns, and is only ever used across multiple tags
     of the same variable.  For instance, if a transport equation is
     being subcycled, the `"NEXT`" tag and the
     `"transport_subcycled_next`" tag may not actually need
     independent copies of the same data.  Once subcycling is complete
     and the outer step is done, they will refer to the same values,
     and can share the same data.  Therefore we create an aliased
     evaluator that indicates that `"concentration@NEXT`" points to
     the same data stored as
     `"concentration@transport_subcycled_next`".

Evaluators are combined into a single **dependency graph**, whereby a
given evaluator points to all of its dependencies.  This graph is a
directed, acyclic graph.  Leaf nodes of the graph are either primary
or independent variable evaluators; internal nodes are secondary or
aliased variable evaluators.  This is a very powerful concept, and
forms one of the two core capabilities that makes Arcos useful.

For example, a simple dependency graph might be used to describe how
the code computes water content.  At runtime, the user might define a
set of evaluators that, together, form a dependecy graph for water
content.

.. image:: static/images/dag_wc.png
  :width: 400
  :alt: Simple dependency graph for water content.

Here and throughout, primary variables are shown in red, independent
variables in blue, and secondary or aliased variables in green.
Arrows in the graph point from an evaluator to its dependencies.

