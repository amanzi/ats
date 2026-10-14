"""Richards water content evaluator: the standard form as a function of liquid saturation."""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("permeability","k"),
        ("mass_density_liquid", "rho"),
        ("viscosity_liquid", "mu"),
        ]
params = [("g", "double", "gravitational constant")]

import sympy
k, rho, mu = sympy.var("k,rho,mu")
g_ = sympy.var("g_")
expression = (k * rho * mu) / mu;

generate_evaluator("hydraulic_conductivity", "ecosim",
                   "hydraulic conductivity", "hydraulic_conductivity",
                   deps, params, expression=expression, doc=__doc__)
