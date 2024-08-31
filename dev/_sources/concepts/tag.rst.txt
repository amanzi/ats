Tags
====

Like domains which can be used to make the same variable distinct
across different parts of space (e.g. `"surface-evaporation`" is
evaporation from the `"surface`" domain, while `"snow-evaporation`" is
evaporation from the `"snow`" domain), it is useful to make a similar
distinction for the same variable but at different times.  Time
**tags** are used to indicate multiple copies of **variables** which
evolve in time.  Tags are used to label different **data**,
or copies of variables, with respect to where they sit on the
reference timestep interval.  Two tags are of particular importance --
the so-called `"global`" tags `"CURRENT`" and `"NEXT`".  These two
tags define the endpoints of the reference timestep.  Individual PDEs
may define other tags their own use; e.g. `"flow_midpoint`" (for a
midpoint time integration scheme) or `"transport_subcycled_next`"
(used to subcycle the transport PK).

Throughout, it can be assumed that the same variable at different tags
always has the same structure, e.g. a field on all cells of a given
mesh, or a single scalar value, etc.  When a step is complete, we can
always assign all variables from their global `"NEXT`" tag's data to
their global `"CURRENT`" tag's data to advance the step, e.g.:

    snow-evaporation\@CURRENT = snow-evaporation\@NEXT

Note the use of the `"@`" delimiter to separate variables from tags.

Tags become particularly important when one realizes that Amanzi-ATS
integrates in time with an adaptive timestep; any PDE is allowed to
fail at advancing to the next timestep.  In a set of coupled
PDEs, if only one fails, we must be able to back up and recover the
initial state, then try again with a smaller timestep.  Therefore, all
variables which cannot be recomputed must have at least two copies --
one for the `"CURRENT`" and one for the `"NEXT`" tags; if another PDE
fails after this PDE was successful, we must not have thrown away the
`"CURRENT`" values.
