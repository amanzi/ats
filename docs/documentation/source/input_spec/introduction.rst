About the Specification
#######################

ATS, and Amanzi's "native" specificiation, is an xml file following
Trilinos's Teuchos ParameterList schema.  There are only two types of
tags used -- `"Parameter`" and `"ParameterList`".  `"Parameter`"
elements consist of `"name`", `"type`", and `"value`" attributes.
`"ParameterList`" elements use the `"name`" attribute and include
subelements that are other `"ParameterList`" and `"Parameter`"
elements.

The top-most, `"main`" list is read by the code and used to provide
all information needed to run the simulation.  This input spec is
designed for the code, not necessarily for the user.  In general,
avoid writing input files from scratch, and prefer to modify existing
demos or examples.

Additionally, the `ATS Input Spec <https://github.com/ecoon/ats_input_spec>`_
tool can be used to generate ATS input files through a python script.
Despite this tool, it is important to be able to parse and understand
the resulting xml file.

Here we document the input spec by defining what each possible element
used by the code needs to be well posed.

Specs
=====

In many cases, the input specifies data for a particular parameterized
model, and ATS supports a number of parameterizations.  For example,
initial data might be uniform (the value is required), or linear in y
(the value and its gradient are required).  Where ATS supports a
number of parameterized models for quantity Z, the available models
will be listed by name, and then will be described in the subsequent
section.  For example, the specification for an `"apple`" list might look
like:

.. _apple-spec:
.. admonition:: apple-spec

   * `"apple kind`" ``[string]`` A string defining the kind of apple;
     valid choices include `"cosmic crisp`", `"granny smith`", or
     `"golden delicious`"
    
   * `"color`" ``[string]`` **red** The color of the apple, defaults to
     `"red`".

An `"apple-spec`" might be used when defining a pie recipe:

.. _apple-pie-spec:
.. admonition:: apple-pie-spec

   * `"apple`" ``[apple-spec]``

   * `"apple volume [cups]`" ``[double]`` **1**

Note that the `"apple`" entry would be a `ParameterList`, named
`"apple`", that includes all of the parameters required to define an
`"apple-spec`".  Also, note that, wherever possible, units are defined
in the name of the parameter.    

An example of using such a specification:

.. code-block:: xml

    <ParameterList name="my apple pie">
      <Parameter name="apple volume [cups]" type="double" value="1.5"/>
      <ParameterList name="apple">
        <Parameter name="apple kind" type="string" value="granny smith"/>
        <Parameter name="color" type="string" value="green"/>
      </ParameterList>   
    </ParameterList>   

or equivalently,

.. code-block:: yaml

   my apple pie:
     apple volume [cups]: 1.5
     apple:
       apple kind: granny smith
       color: green
       

Syntax
======

* Reserved keywords and labels are `"quoted and italicized`" -- these
  labels or values of parameters in user-generated input files must
  match (using XML matching rules) the specified or allowable values.

* User-defined labels are indicated with ALL-CAPS, and are meant to
  represent a typical or default name given by a user - these can be
  names or numbers or whatever serves best the organization of the
  user input data.  Things liked PRESSURE or SURFACE-PONDED_DEPTH can
  be renamed from their defaults if it makes sense to the problem.

* Bold values are default values, and are used if the Parameter
  is not provided.

  
Units
=====

Wherever possible, units for a parameter are included in the parameter
name, in brackets, and in SI units.  When there are common alternative
units (e.g. g/mL) we sometimes provide a second parameter option with
the same name, but different units.

All units are written in the numerator.  `^` is used to indicate
powers, including negative powers.  For example, a parameter
prescribing the density of water might be: `"density of water [kg
m^-3]`".


Variable Naming
===============

Unlike parameters (which can only be spelled one way), variables are
often named with user-defined strings which have a default value.  For
instance, a Richards equation might need a variable defining the
saturation of the liquid phase.  By default, we would call this
`"saturation_liquid`", but the user may choose to override this name.
This makes variables hard to define with certainty, because the user
may change these things.  This is a feature, but to a new user, it may
seem confusing.  Prefer to use default names wherever possible, unless
there is a good reason not to do so.

Variables are named according to a very strong convention.  While
variable names may be overridden by the user, users should choose to
follow these conventions or things like visualization scripts may not
behave as expected.

A **variable name** looks like one of:

- ROOT_NAME
- DOMAIN-ROOT_NAME
- DOMAIN_SET:ID-ROOT_NAME

where:

- ROOT_NAME is the root of the variable, e.g. `"temperature`".
- When DOMAIN is not supplied, it is implied to be the "default" mesh,
  called `"domain`" in the mesh list.  Otherwise the domain name is
  the same as the mesh name (e.g. `"surface`").  This distinguishes
  between `"surface-temperature`" (the temperature on the surface
  mesh) and `"temperature`" (the temperature on the subsurface/domain
  mesh).
- DOMAIN_SET:ID is itself a DOMAIN, where the set defines the
  collection as a whole (from the mesh list) and the ID is defined by
  an index across the collection, e.g. `"column:4`".  So
  `"column:4-temperature`" is the temperature on the 4th domain of the
  domain set `"column`".

Additionally, a variable name may include a **component name** and a
**subfield name**, e.g. ROOT_NAME.COMPONENT_NAME.SUBFIELD_NAME.  The
component name refers to a component of the field, and typically is
named according to the entity on which it is defined, e.g. `"cell`" or
`"face`".  The SUBFIELD_NAME refers to a name for each degree of
freedom in the vector.  Typically this is used for n-dimensional
vector-valued fields such as velocity, where the names are e.g. `"x`",
`"y`", and `"z`", or chemical specie names.  If multiple degrees of
freedom are needed but subfield names are not provided, they are
assigned integer values indexing the degree of freedom.
  
**Tags** indicate the use of a variable at a specific time in the
discretized time interval.  Default tags include `"current`" and
`"next`" indicating the variable at the beginning and end of the
interval, respectively.  Often subcycling and other schemes will
designate special-purpose tags which may be used internally by a
subset of the equations begin solved.  Tags are combined with
variables to indicate a specific data structure,
e.g. `"surface-pressure@next`".  Note that the default tag `"`" is
equivalent to `"next`".

Lastly, **derivatives** are named using the `"d`" and the `"|`"
character, e.g. `"dsurface-water_content|dsurface-pressure`" is the
derivative of the `"water_content`" variable on the `"surface`" domain
with respect to the `"pressure`" on the same domain.

As a result of these conventions, none of the above individual
strings, (root names, domains, domain sets, or IDs) can contain any of
the following reserved characters: `:`, `-`, `|`, `@`.  The lone
exception to this rule is that subfield names may include trailing `-`
characters to indicate charge,
e.g. `total_component_concentration.NO3-`.

Avoid including spaces in variable names.  This doesn't break the
code, but it does break some visualization tools.


Default Names, Symbols, and Units
=================================

.. include:: symbol_table.rst

.. include:: symbol_table_notes.rst             
