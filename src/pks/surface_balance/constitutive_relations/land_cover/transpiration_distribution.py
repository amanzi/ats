"""Distributes transpiration based upon a rooting depth and a wilting-point water-potential factor.
"""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("pressure", "p"),
        ("rooting_depth_fraction", "f_root"),
        ("transpiration", "trans_total")]
params = [("wp_min", "double", "wilting point [Pa]", -2.0e6),]

generate_evaluator("transpiration_distribution", "Flow",
                   "transpiration distribution via rooting depth", "transpiration",
                   deps, params, doc=__doc__)
