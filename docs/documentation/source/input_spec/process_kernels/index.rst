Process Kernels
###############
.. _`PK`:

Process Kernels, or PKs, are the fundamental unit of a model, and
represent a single or system of Partial Differential Equations (PDEs)
or Differential Algebraic Equations (DAEs).  PKs are broadly split
into individual equations (Physical PKs) and systems of equations,
called Multi-Process Coordinators (MPCs).

The PK tree forms the fundamental definition of the entire system of
equations to be solved by the simulator, and is represented by a
single PK or a single MPC which couples other MPCs and/or Physical
PKs.

.. toctree::
   :maxdepth: 1

   base
   physical/index
   mpcs/index
   
