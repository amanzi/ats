Engineered Water Models 
------------------------

Gate Structure
^^^^^^^^^^^^^^^

`src/physics/ats/src/pks/flow/constitutive_relations/sources/surface_gate_structure_evaluator.hh <https://github.com/amanzi/ats/blob/master/src/pks/flow/constitutive_relations/sources/surface_gate_structure_evaluator.hh>`_

Simulates gravity-driven water movement/diversions through gate structure.

This evaluator models gate structure inside 2D flow area. The gate curve is kept general and can be provided by the user.
Gate structure can be used to move water two canals, two storage areas or canal to storage area.

`"evaluator type`" = `"gate structure`"

.. _evaluator-gate-structure-spec:
.. admonition:: evaluator-gate-structure-spec

   * `"gate intake region`" ``[str]`` Region of cells where gate flow is taken out.
   * `"storage area region`" ``[str]`` Region of cells where gate flow is introduced.
   * `"gate stage close`" ``[double]`` The water surface elevation in the storage area that should trigger the gate close.
   * `"is ponded depth function`" ``[bool]`` If true, the gate efficiency curve is a function of ponded depth in the intake region rather than stage.
   * `"function" ``[function-tabular]`` This is a function/table of head and flow. Here, head is the stage or ponded depth in the intake region (upstream) of the gate structure.

   KEYS:
   - `"cell volume`" **DOMAIN-cell_volume** 
   - `"ponded depth`" **DOMAIN-ponded_depth** 
   - `"potential`" **DOMAIN-pres_elev** stage or water surface elevation
   - `"elevation`" **DOMAIN-elevation** 
   - `"water content`" **DOMAIN-water_content** 
   - `"molar liquid density`" **DOMAIN-molar_density_liquid** 

Example:

.. code-block:: xml

      <ParameterList name="surface-gate_flow" type="ParameterList">
        <Parameter name="evaluator type" type="string" value="gate structure"/>
        <Parameter name="gate intake region" type="string" value="gate intake area"/>
        <Parameter name="storage area region" type="string" value="detention pond"/>
        <Parameter name="gate close stage" type="double" value="7.50"/>
        <Parameter name="is ponded depth function" type="bool" value="true"/>
        <ParameterList name="function" type="ParameterList">
          <ParameterList name="function-tabular" type="ParameterList">
            <Parameter name="x values" type="Array(double)" value="{0.0, 0.18, 0.1840954, 0.20202941, 0.21996341, 0.23789742, 0.25583143}"/>
            <Parameter name="y values" type="Array(double)" value="{0.0, 0.0, 0.01542766, 0.07531978, 0.1352119, 0.19510401, 0.25499613}"/>
            <Parameter name="forms" type="Array(string)" value="{linear, linear, linear, linear, linear, linear}"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>




Pump System
^^^^^^^^^^^^

`src/physics/ats/src/pks/flow/constitutive_relations/sources/surface_pump_system_evaluator.hh <https://github.com/amanzi/ats/blob/master/src/pks/flow/constitutive_relations/sources/surface_pump_system_evaluator.hh>`_

An evaluator for stage-based pump systems

This evaluator models stage-based pump station model inside 2D flow area. 
Pump stations can be used to move water between any combination of river reaches, storage areas or catchment regions. 
Based on pump on-off conditions and pump-operation curve, water is moved from pump-inlet to -outlet region instantly.

`"evaluator type`" = `"pump system`"

.. _evaluator-pump-system-spec:
.. admonition:: evaluator-pump-system-spec

   * `"pump inlet region`" ``[str]`` Region of cells where pump flow is taken out.
   * `"pump outlet region`" ``[str]`` Region of cells where pump flow is introduced (optional).
   * `"on off reference region`" ``[str]`` Region used to determine when the pump should turn on or off (optional). Defaults to "pump inlet region".
   * `"pump start at stage`" ``[double]`` The water surface elevation that should trigger the pump start.
   * `"pump stop at stage`" ``[double]`` The water surface elevation that should trigger the pump stop. 
   * `"maximum pumpline elevation`" ``[double]`` Allows the user to enter the highest elevation in the pump line. E.g., pumping water over top of a levee.
   * `"function`" ``[function-tabular]`` This is a function/table of head and pump flow. Here, head is head difference between outlet and inlet (or head by which the water is to be lifted).

   KEYS:
   - `"cell volume`" **DOMAIN-cell_volume** 
   - `"ponded depth`" **DOMAIN-ponded_depth** 
   - `"potential`" **DOMAIN-pres_elev** stage or water surface elevation
   - `"elevation`" **DOMAIN-elevation** 
   - `"water content`" **DOMAIN-water_content** 
   - `"molar liquid density`" **DOMAIN-molar_density_liquid** 
   - `"pump on`" **DOMAIN-pump_on_flag** status of pump

Example:

.. code-block:: xml

  <ParameterList name="surface-pump_flow" type="ParameterList">
    <Parameter name="evaluator type" type="string" value="pump system"/>
    <Parameter name="pump inlet region" type="string" value="outlet pump intake area"/>
    <Parameter name="pump start at stage" type="double" value="8.18"/>
    <Parameter name="pump stop at stage" type="double" value="8.16"/>
    <Parameter name="maximum pumpline elevation" type="double" value="8.4"/>
    <ParameterList name="function" type="ParameterList">
      <ParameterList name="function-tabular" type="ParameterList">
        <Parameter name="x values" type="Array(double)" value="{0.0, 0.01, 0.05, 0.1, 0.14, 0.1740954, 0.21, 0.225, 0.235, 0.24}"/>
        <Parameter name="y values" type="Array(double)" value="{0.3, 0.3, 0.28, 0.23, 0.18, 0.1352119, 0.076, 0.045, 0.02, 0.0}"/>
        <Parameter name="forms" type="Array(string)" value="{linear, linear, linear, linear, linear, linear, linear, linear, linear}"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>





Culverts
^^^^^^^^^

`src/physics/ats/src/pks/flow/constitutive_relations/sources/surface_culvert_evaluator.hh <https://github.com/amanzi/ats/blob/master/src/pks/flow/constitutive_relations/sources/surface_culvert_evaluator.hh>`_


Simulates water movement through culverts by instant transfer of water from inlet to outlet region.
Flow is calculated using standard culvert hydraulics, considering both inlet-controlled and outlet-controlled regimes.

Implements the following culvert hydraulics equations:

   - Inlet control:
     \f[
     Q_{inlet} = N_b C A \sqrt{2g h_i}
     \f]
     where:
       - \( N_b \) = number of barrels  
       - \( C \) = discharge coefficient  
       - \( A \) = culvert cross-sectional area  
       - \( h_i \) = head at culvert inlet  
       - \( g \) = gravity

   - Outlet control:
     \f[
     Q_{outlet} = N_b C A \sqrt{ \frac{2g h_o}{k} }
     \f]
     where:
       - \( h_o \) = head difference between inlet and outlet  
       - \( k = 1.5 + \frac{29 n^2 L}{R^{4/3}} \) (Manning-based resistance term)  
       - \( n \) = Manning's roughness  
       - \( L \) = culvert length  
       - \( R \) = hydraulic radius

   - Blended total discharge:
     \f[
     Q = \frac{Q_{inlet} Q_{outlet}}{\sqrt{Q_{inlet}^2 + Q_{outlet}^2 + \epsilon}}
     \f]
     where \( \epsilon \) is a small number to avoid divide-by-zero.

   The resulting \( Q \) is used to compute area-weighted water removal at the inlet and volume-weighted water addition at the outlet.

.. _evaluator-culvert-spec:
.. admonition:: evaluator-culvert-spec

   * `"culvert inlet region"`" ``[str]`` Region of cells where culvert flow is taken out.
   * `"culvert outlet region"`" ``[str]`` Region of cells where culvert flow is introduced.
   * `"number of barrels"`" ``[int]`` Number of culvert barrels, default is 1.
   * `"culvert length"`" ``[double]`` Length of the culvert in meters, default is 10.
   * `"culvert diameter"`" ``[double]`` Diameter of the culvert in meters, default is 1.
   * `"culvert roughness coefficient"`" ``[double]`` Manning's roughness coefficient for the culvert, default is 0.013.
   * `"culvert discharge coefficient"`" ``[double]`` Discharge coefficient for the culvert, default is 0.6.

   KEYS:
   - `"cell volume"`" **DOMAIN-cell_volume** 
   - `"ponded depth"`" **DOMAIN-ponded_depth** 
   - `"potential"`" **DOMAIN-pres_elev** (stage or water surface elevation)
   - `"elevation"`" **DOMAIN-elevation** 
   - `"water content"`" **DOMAIN-water_content** 
   - `"molar density liquid"`" **DOMAIN-molar_density_liquid**

Example:

.. code-block:: xml

      <ParameterList name="surface-culvert_flow" type="ParameterList">
        <Parameter name="evaluator type" type="string" value="culvert"/>
        <Parameter name="culvert inlet region" type="string" value="culvert inlet"/>
        <Parameter name="culvert outlet region" type="string" value="culvert outlet"/>
        <Parameter name="number of barrels" type="int" value="1"/>
        <Parameter name="culvert length" type="double" value="10.0"/>
        <Parameter name="culvert diameter" type="double" value="0.25"/>
        <Parameter name="culvert roughness coefficient" type="double" value="0.013"/>
        <Parameter name="culvert discharge coefficient" type="double" value="0.6"/>
      </ParameterList>





Impervious Interception
^^^^^^^^^^^^^^^^^^^^^^^^

`src/physics/ats/src/pks/flow/constitutive_relations/sources/impervious_interception_evaluator.hh <https://github.com/amanzi/ats/blob/master/src/pks/flow/constitutive_relations/sources/impervious_interception_evaluator.hh>`_

An evaluator for rerouting precip due to impervious surface.

This evaluator reroutes incoming precipitation, taking a portion of it (where
the portion is determined by the impervious area fraction) and moving it
(instantly) into the nearby stream network.

Note: this assumes that the runoff reciever is constant in time!

.. _impervious-interception-evaluator-spec:
.. admonition:: impervious-interception-evaluator-spec

   * `"maximum specific diversion rate [m s^-1]`" ``[double]`` **inf**
     Maximum rate of water removal through storm drains, etc, in units of m^3
     water per second per m^2 of _impervious_ area (specific area).

   KEYS:
   - `"impervious fraction`" **DOMAIN-impervious_fraction** The fraction of
     surface area that is impervious, this also defines the fraction of precip
     that is rerouted.
   - `"impervious runoff receiver`" **DOMAIN-impervious_runoff_receiver`" The
     Global ID of the cell that will recieve water from this cell.
   - `"incoming water source`" **DOMAIN-precipitation_rain** The source of
     water to be re-reouted -- this is typically rain, but might be
     canopy-throughfall_drainage_rain, and might be snow-melt, etc.
   - `"cell volume`" **DOMAIN-cell_volume**


