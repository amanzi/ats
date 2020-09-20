"""Aspect modified shortwave radiation is determined by a factor which
is multiplied by the 'incoming radiation incident on a flat surface'
to determine the 'incoming radiation incident on a sloping surface of
a given aspect' as a function of latitude, slope, aspect, and Julian
day of the year, and time of day.

Note that some careful checking and experimentation has found that, in
general, the daily average incident radiation times the 12-noon aspect
modifier correlates reasonably well with the daily average of the
product of the hourly incident radiation and the hourly aspect
modifier.  It is notably better than the daily average radiation times
the daily average aspect modifier.

"""

import sys, os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], "tools", "evaluator_generator"))
from evaluator_generator import generate_evaluator

deps = [("slope", "slope"),
        ("aspect", "aspect"),
        ("incoming_shortwave_radiation", "qSWin"),
        ]
params = [("daily_avg", "bool", "estimate factor for a daily averaged radiation"),
          ("lat", "double", "domain-averaged latitude [degrees]"),
          ("doy0", "int", "day of year at time 0 [julian days]")]

generate_evaluator("incident_shortwave_radiation", "SurfaceBalance",
                   "incident shortwave radiation", "incident_shortwave_radiation",
                   deps, params, doc=__doc__)
                   


