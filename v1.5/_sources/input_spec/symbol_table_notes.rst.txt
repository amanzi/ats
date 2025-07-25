.. note::

   1. This is incorrectly named in ATS, as it is not what is
      traditionally called the relative_permeability.  As ATS works in
      a pressure basis (as opposed to a head basis), it uses Darcy
      equations that use permeability, not conductivity. The
      "diffusion" coefficient in Darcy's equation in pressure form is
      :math:`\frac{n}{\mu}k_r K`, where the density, viscosity, and
      relative permeability are scalars that are typically upwinded or
      averaged to faces, while the absolute permeability may be a
      tensor. We therefore store the first three terms together, but
      incorrectly call this "relative_permeability." Furthermore,
      because K is order :math:`10^{\{-10 \dash -15\}}`, but
      :math:`\frac{{n}}{{\mu}}` is order :math:`10^7`, we are
      multiplying a very small number by a very large number, a
      classic problem in numerics.  Therefore, we typically rescale
      both, moving 7 orders of magnitude off of the scalar and putting
      them on the tensor.  As a result, the units of this variable are
      something like: :math:`[mol m^-3 Pa^-1 s^-1 10^7]`, and typical
      values range from 0 to ~6.  Note that the rescaling factor is
      NOT stored on the absolute permeability, so permeability is in
      the typical units :math:`[m^2]`.
   2. The total energy source includes both direct sources of energy
      (e.g. warming/cooling, radiation, etc), but also sources/sinks
      of internal energy due to sources/sinks of mass.  It does not
      include fluxes of energy, e.g. diffusion or advection of energy.
   3. We use the word "extensive" to mean a quantity that is measuring
      the quantity, and is not per unit grid cell volume or area.
   4. Here the number indicates the coordinate dimension, e.g. x,y,z.
   5. Conserved quantity for a PK.
   6. We use the word "specific" to mean a quantity that is per unit
      extent, e.g. specific enthalpy is per unit mol of water, or
      specific leaf area is per unit dry mass.

