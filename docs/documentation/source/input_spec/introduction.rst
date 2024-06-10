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
section.  For example, the specification for an `"X`" list might look
like:

.. _X-spec:
.. admonition:: X-spec

  * `"Y`" ``[string]`` **default_value** Documentation desribing Y.
  * `"Z`" ``[Z-spec]`` Model for Z, One of `"z1`" or `"z2`" (see below) 

Here, an `"X`" is defined by a `"Y`" and a `"Z`".  The `"Y`" is a
string parameter but the `"Z`" is given by a model (which will require
its own set of parameters).  The options for `"Z`" will then be
described seperately as a `"Z-spec`"


An example of using such a specification:

.. code-block:: xml

    <ParameterList name="X">
      <Parameter name="Y" type="string" value="hello"/>
      <ParameterList name="z2">
        <Parameter name="z2a" type="double" value="0.7"/>
        <Parameter name="z2b" type="int" value="3"/>
      </ParameterList>   
    </ParameterList>   
 

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

Naming
======

Variables are named according to a very strong convention.  While some
variables may be overridden by the user, users should choose to follow
these conventions or things like visualization scripts may not behave
as expected.

A variable name looks like one of:

- ROOT_NAME
- DOMAIN-ROOT_NAME
- DOMAIN_SET:ID-ROOT_NAME

where:

- When DOMAIN is not supplied, it is implied to be the "default" mesh,
  called `"domain`" in the mesh list.  Otherwise the domain name is
  the same as the mesh name (e.g. `"surface`").
- DOMAIN_SET:ID is itself a DOMAIN, where the set defines the
  collection as a whole (from the mesh list) and the ID is defined by
  an index across the collection, e.g. `"column:4`"

Tags indicate the use of a variable at a specific time in the
discretized time interval.  Default tags include `"current`" and
`"next`" indicating the variable at the beginning and end of the
interval, respectively.  Often subcycling and other schemes will
designate special-purpose tags which may be used internally by a
subset of the equations begin solved.  Tags are combined with
variables to indicate a specific data structure,
e.g. `"surface-pressure@NEXT`".

Lastly, derivatives are named using the `"d`" and the `"|`" character,
e.g. `"dsurface-water_content|dsurface-pressure`" is the derivative of
the `"water_content`" variable on the `"surface`" domain with respect
to the `"pressure`" on the same domain.

As a result of these conventions, none of the above individual strings,
(root names, domains, domain sets, or IDs) can contain any of the
following reserved characters: `:`, `-`, `|`, `@`.
  

Name and Symbol Index
=====================

.. include:: symbol_table.rst

.. include:: symbol_table_notes.rst             
