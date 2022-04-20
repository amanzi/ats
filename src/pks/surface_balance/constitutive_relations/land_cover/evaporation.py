"""Downregulates evaporation from a potential.
"""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("saturation_gas", "sg"),
        ("porosity", "poro"),
        ("potential_evaporation", "pot_evap")]
params = [
    ("dess_dz", "double", "dessicated zone thickness [m]", 0.1),
    ("Clapp_Horn_b", "double", "Clapp and Hornberger b of surface soil [-]", 1.0),]

generate_evaluator("evaporation_downregulation", "SurfaceBalance",
                   "evaporation downregulation via soil resistance", "evaporation",
                   deps, params, doc=__doc__)
