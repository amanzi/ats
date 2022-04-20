"""Wilting factor.

Beta, or the water availability factor, or the plant wilting factor.

Beta =  (p_closed - p) / (p_closed - p_open)

where p is the capillary pressure or soil mafic potential, and closed
and open indicate the values at which stomates are fully open or fully
closed (the wilting point).

"""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("capillary_pressure_gas_liq", "pc"),]
params = [("pc_o", "double", "capillary pressure at fully open stomates [Pa]"),
          ("pc_c", "double", "capillary pressure at wilting point [Pa]")]

import sympy
z = sympy.var("pc")
pc_o_, pc_c_ = sympy.var("pc_o_,pc_c_")

# FIXME -- needs a max of 1, min of 0!
expression = (pc_c_ - pc) / (pc_c_ - pc_o_)

generate_evaluator("plant_wilting_factor", "Flow",
                   "plant wilting factor", "plant_wilting_factor",
                   deps, params, expression=expression, doc=__doc__)
