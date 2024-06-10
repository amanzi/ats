ATS Demonstrations
==================

This is a suite of demonstration problems which show the various
capabilities and typical workflows of ATS simulations.

A wide range of capabilities are shown here -- this is really meant to
be the first introduction to ATS and the primary introductory
documentation for ATS.  **Start here!**

Included in each demonstration folder is a collection of runs meant to
map out usage of the described physical capability.  These include
everything needed to get started using ATS.  Any given run contains
some of, and maybe all of, the following components:

* An input file: an `.xml` file read by ATS via `ats
  ../my_file.xml` **required**

* A mesh file for non-trivial meshes: typically an `.exo` file.  Most
  meshes have included a python script or ipython notebook used to
  generate that mesh, showing that workflow.

* Forcing data: typically an `.h5` file, and often generated from
  DayMet or other products.  Often a python script or ipython notebook
  is used to generate these files, and these are included to show that
  workflow as well.

* A jupyter notebook `.ipynb` which describes, documents, and plots
  the results from the given run.  Note these can be viewed, including
  plots, without even running the simulation.  This is where to
  understand what is being done in each run.

Running the demos
---------------------

Running all of the demos can take some time, but individual demos are
often fairly quick.  To run a given demo, make sure `ats` is compiled
and (preferably) in your path, or that you have a docker container ready to go, and then:

```python run_demos.py path_to_suite.cfg```

or

```python run_demos.py path_to_suite.cfg -t test_name```

or (with docker):

```python run_demos.py path_to_suite.cfg -e metsi/ats:latest -t test_name```

Note that some individual runs may depend upon results from previous
runs in that suite, and so all demos in that suite should be run.
This is particularly true for demos that show a full workflow, such as
`ecohydrology` or `arctic_hydrology`.

Also, if you want to run the demos without using these python scripts,
be sure to check out the FAQs on our `Wiki <https://github.com/amanzi/ats/wiki>`_


Visualizing the results
------------------------

Inside each subdirectory is a jupyter notebook.  Jupyter comes fairly
standard with most comprehensive python installations.  Anaconda
python is the one most ATS developers use, and its default
installation includes nearly all python packages used by ATS.


Demonstration Problems
----------------------

.. inclusion-marker


Richards Equation: Steady state
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. toctree::
   :maxdepth: 1
   
   ats_demos/01_richards_steadystate/richards_steadystate.ipynb

This shows examples of solving Richards equation to steadystate.
Typically this is used to establish a water column that satisfies
hydrostatic balance.


Richards Equation: Transient
>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. toctree::
   :maxdepth: 1

   ats_demos/02_richards/richards.ipynb

Transient problems show a variety of variably saturated cases, and
demonstrate seepage faces and other common boundary conditions for the
flow of water in a porous media.


Surface Water
>>>>>>>>>>>>>

.. toctree::
   :maxdepth: 1

   ats_demos/03_surface_water/surface_water.ipynb

Overland flow is solved through a diffusion wave equation.  This
demonstrates that as a standalone capability, solving surface water
problems to demonstrate common usages of boundary conditions and
forcing.


Integrated Hydrology
>>>>>>>>>>>>>>>>>>>>

.. toctree::
   :maxdepth: 1

   ats_demos/04_integrated_hydro/integrated_hydro.ipynb

Integrated hydrology brings the previous two examples together,
solving both surface and subsurface flow of water.


Ecohydrology
>>>>>>>>>>>>

.. toctree::
   :maxdepth: 1

   ats_demos/05_ecohydrology/ecohydrology.ipynb

Ecohydrogy brings in the effects of other ecological processes, here
loosely used to include all surface processes like evaporation,
transpiration, and canopy processes like interception and storage, and
even simplified biogeochemistry processes for a full carbon cycle.


Arctic Hydrology
>>>>>>>>>>>>>>>>

.. toctree::
   :maxdepth: 1

   ats_demos/06_arctic_hydrology/arctic_hydrology.ipynb

ATS was originally developed as an Arctic hydrology simulator.  It
includes state-of-the art constitutive models and numerical methods
for solving coupled freeze-thaw processes in Arctic environments.


Reactive Transport
>>>>>>>>>>>>>>>>>>

(Work in progress)

.. toctree::
   :maxdepth: 1

   ats_demos/07_reactive_transport/reactive_transport.ipynb 

ATS's sister code Amanzi was designed for solving problems of reactive
transport.  Interoperability of ATS and Amanzi allows ATS to leverage
this work to solve problems of nonreactive and reactive transport in
both the surface and subsurface, and even in frozen environments.
(See below for ATS reactive transport demos when used with integrated hydrology)

Integrated hydrology and Reactive Transport
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. toctree::
   :maxdepth: 1

   ats_demos/13_integrated_hydro_reactive_transport/integrated_hydro_reactive_transport.ipynb

ATS is unique in the its ability to simulate reactive transport in integrated
hydrology problems. In other words, it is capable of simulating muticomponent
reactive transport in both surface and subsurface compartments using a novel
coupling approach, and levering powerful, external geochemical engines
`Molins et al (2022) WRR <https://doi.org/10.1029/2022WR032074>`_
