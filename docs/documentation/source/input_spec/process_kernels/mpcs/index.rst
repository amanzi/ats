MPC
===

Multi-process-couplers or MPCs couple other PKs.  They also are PKs
themselves, in that they implement the PK interface.  So MPCs can also
couple other MPCs.  There are a few common "base" MPCs which do the
simplest form of coupling -- sequential and globally implicit (with a
diagonal preconditioner).  Then there are specific couplers which know
more about their coupled sub-PKs, and can do more complicated things
(for instance, adding off-diagonal block entries to the
preconditioner).

MPCs are also used to couple across domains -- for instance integrated
hydrology is a surface+subsurface flow coupler.  They also can do
fancier things like drape a bunch of subgrid columns off of a mesh, or
other things.  Think of these as the custom couplers.

.. toctree::
   :maxdepth: 1

   base
   tightly_coupled
   multiple_domain
   globalization
