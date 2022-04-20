"""Rooting depth function.

Sets the root fraction as a function of depth,

F_root =  ( a*exp(-az) + b*exp(-bz) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.

"""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("depth", "z"),]
params = [("a", "double", "alpha", 7),
          ("b", "double", "beta", 1.75),
          ("z_max", "double", "max rooting depth [m]", 2.0)]

import sympy
z = sympy.var("z")
a_, b_, z_max_ = sympy.var("a_,b_,z_max_")
expression = sympy.Piecewise((0.5 * (a_ * sympy.exp(-a_*z) + b_ * sympy.exp(-b_*z)), z <= z_max_), (0., True))

generate_evaluator("rooting_depth_fraction", "Flow",
                   "rooting depth fraction", "rooting_depth_fraction",
                   deps, params, expression=expression, doc=__doc__)
