On use of VerboseObject:
----------------------------

The following levels are available to a user: `"none`", `"low`",
`"medium`", `"high`", and `"extreme`".  Please use the following
guidelines when using these levels in your code:

- `Teuchos::VERB_NONE` Really, don't ever use this.  Don't do output at all.

- `Teuchos::VERB_LOW` Should be used by time step-level work only.  It
  should print a message saying "start a timestep" and the result of
  that timestep (success/fail).  So really almost exclusively in a
  PK's AdvanceStep() and ValidStep() methods, if they exist.
  Evaluators should not use this level.  Note weak MPCs are the
  exception to this -- they shouldn't write this type of information
  because their coupled components should.  Strong MPCs on the other
  hand do not call their coupled components AdvanceStep() and
  therefore should.

- `Teuchos::VERB_MEDIUM` Should write things at the solver/iteration
  level.  So for instance, this should print out error norms at each
  iteration of the solver, but little else.  So a PK's ErrorNorm()
  should use this.  Evaluators should not use this level.

- `Teuchos::VERB_HIGH` Should write things at the physics levels.  PKs
  should write `"debug cell`" information, preconditioner information,
  etc.  Evaluators should write `"debug cell`" information on update.

- `Teuchos::VERB_EXTREME` Debugging output.  PKs should write
  something for every method called unless it is empty.  Evaluators
  should write their dependency graph checks for whether they must be
  updated or not.


Suggested use of `VerboseObject` looks like:

```
Teuchos::OSTab tab = vo_->getOSTab();
if (vo_->os_OK(Teuchos::VERB_LOW))
  *vo_->os_OK() << "my message" << std::endl;
```




