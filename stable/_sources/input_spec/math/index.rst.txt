Time integrators, solvers, and other mathematical specs
#######################################################

  Common specs for all solvers and time integrators, used in PKs.

There are three commonly used broad classes of time integration
strategies.

"Update" methods are the simplest -- they use no formal
mathematical definition of a differential equation, but instead
implicitly use a process by which variables at the new time are
directly calculated.  Typically there *is* an implied ODE or PDE here,
but it is not stated as such and time integration routines are not
used.  Examples of these are common in biogeochemistry and vegetation
models.

"Explicit" time methods are the next simplest.  These include a
variety of options from forward Euler to higher order Runge-Kutta
schemes.  These only require evaluating forward models where we have
existing of the dependencies.  If they work, these are great thanks to
their deterministic nature and lack of expensive, memory-bandwith
limited solvers.  But they only work on some types of problems.
Examples of of these include transport, where we use high order time
integration schemes to preserve fronts.

"Implicit" and semi-implicit methods instead require the evaluation of
a residual equation -- the solution is guessed at, and the residual is
calculated, which measures how far the equation is from being
satisfied.  This measure is then inverted, finding a correction to the
guess which hopefully reduces the residual.  As the residual goes to
zero, the error, a measure of the difference between the guess and the
true solution, also goes to zero.  To do this inversion, we lean on
Newton and Newton-like methods, which attempt to somehow linearize, or
approximately linearize, the residual function near the guess in order
to calculate an update.  In this case, the time integration scheme
requires both a nonlinear solver (to drive the residual to zero) and a
linear solver or approximate solver (to calculate the correction).

.. toctree::
   :maxdepth: 1

   time_integrators
   nonlinear_solvers
   linear_solvers
   preconditioners
   operators
