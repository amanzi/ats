/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* Simulation controller and top-level driver, defined in *.xml, For ONE E3SM Land Model (ELM) Timestep

     (1) ALL inputs for ELM-ATS interface, inc. *.xml, required by ATS would be in E3SM's input data directory
     under 'lnd/clm2/ats'.

     (2) ALL parameter list from *.xml ARE merely a data placeholder, by which ELM may reset or override,
     By Methods defined in this class. For an example, not matter what values of 'start time' and 'end time'
     in cycle_drvier, ELM will over-ride them each ELM timestep.

 Authors: F.-M. Yuan (yuanf@ornl.gov), Ethan Coon (coonet@ornl.gov)
*/

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

In the `"cycle driver`" sublist, the user specifies global control of the
simulation, including starting and ending times and restart options.

.. _coordinator-spec:
.. admonition:: coordinator-spec

    * `"start time`" ``[double]`` **0.** Specifies the start of time in model time.
    * `"start time units`" ``[string]`` **"s"** One of "s", "d", or "yr"

    ONE OF

    * `"end time`" ``[double]`` Specifies the end of the simulation in model time.
    * `"end time units`" ``[string]`` **"s"** One of `"s`", `"d`", or `"yr`"

    OR

    * `"end cycle`" ``[int]`` **optional** If provided, specifies the end of the
      simulation in timestep cycles.

      END

    * `"subcycled timestep`" ``[bool]`` **false**  If true, this coordinator creates
      a third State object to store intermediate solutions, allowing for failed
      steps.
    * `"restart from checkpoint file`" ``[string]`` **optional** If provided,
      specifies a path to the checkpoint file to continue a stopped simulation.
    * `"wallclock duration [hrs]`" ``[double]`` **optional** After this time, the
      simulation will checkpoint and end.
    * `"required times`" ``[io-event-spec]`` **optional** An IOEvent_ spec that
      sets a collection of times/cycles at which the simulation is guaranteed to
      hit exactly.  This is useful for situations such as where data is provided at
      a regular interval, and interpolation error related to that data is to be
      minimized.
    * `"PK tree`" ``[pk-typed-spec-list]`` List of length one, the top level
      PK_ spec.

Note: Either `"end cycle`" or `"end time`" are required, and if
both are present, the simulation will stop with whichever arrives
first.  An `"end cycle`" is commonly used to ensure that, in the case
of a time step crash, we do not continue on forever spewing output.

Example:

.. code-block:: xml

   <ParameterList name="cycle driver">
     <Parameter  name="end cycle" type="int" value="6000"/>
     <Parameter  name="start time" type="double" value="0."/>
     <Parameter  name="start time units" type="string" value="s"/>
     <Parameter  name="end time" type="double" value="1"/>
     <Parameter  name="end time units" type="string" value="yr"/>
     <ParameterList name="required times">
       <Parameter name="start period stop" type="Array(double)" value="{0,-1,86400}" />
     </ParameterList>
     <ParameterList name="PK tree">
       <ParameterList name="my richards pk">
         <Parameter name="PK type" type="string" value="richards" />
       </ParameterList>
     </ParameterList>
   </ParameterList>

*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "elm_ats_data.hh"

namespace ATS {

class elm_ats_plist {

  public:

    //
    elm_ats_plist(Teuchos::RCP<Teuchos::ParameterList> &plist);
    ~elm_ats_plist() = default;

    //
    void set_plist(elm_data &elmdata, const double start_ts, const double dt);

  private:

  // passing data via 'parameter list' (prior to model setup)
  void plist_general_mesh_reset(elm_data &elmdata, const bool elm_matched=false);
  void plist_materials_reset(elm_data &elmdata);
  void plist_cycle_driver_reset(const double start_ts, const double dt);

  Teuchos::RCP<Teuchos::ParameterList> parameter_list_;
  Teuchos::RCP<Teuchos::ParameterList> elm_drv_plist_;
  Teuchos::RCP<Teuchos::ParameterList> elm_pks_plist_;
  Teuchos::RCP<Teuchos::ParameterList> elm_state_plist_;
  Teuchos::RCP<Teuchos::ParameterList> subpk_plist_;
  Teuchos::RCP<Teuchos::ParameterList> srfpk_plist_;

}; // close of class elm_ats_plist

} //close of namespace ATS

