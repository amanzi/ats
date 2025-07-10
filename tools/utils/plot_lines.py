import numpy as np
from matplotlib import pyplot as plt
import matplotlib.collections
import matplotlib.ticker

def plotLines(x, y, t, ax=None, colorbar=True, colorbar_ticks=True, colorbar_label=None, t_min=None, t_max=None, **kwargs):
    """Plots lines by color."""

    # create an axis if none
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    # set bounds
    if t_min is None:
        t_min = t[0]
    if t_max is None:
        t_max = t[-1]

    # plot the lines
    segments = [np.column_stack((x,y)) for (x,y) in zip(x,y)]
    lc = matplotlib.collections.LineCollection(segments, **kwargs)
    c = (t - t_min) / (t_max - t_min)
    lc.set_array(c)

    # add to axis
    ax.add_collection(lc)
    ax.autoscale()

    if colorbar:
        # create the colorbar
        cb = fig.colorbar(lc)

        # label it
        if colorbar_label is not None:
            cb.set_label(colorbar_label)

        if colorbar_ticks:
            # transform ticks from [0,1] to [times[0], times[-1]] to get
            # the correct labels
            def my_formatter(x, pos):
                xp = x * (t_max - t_min) + t_min
                if abs(xp - np.round(xp)) < 1.e-2:
                    return str(int(np.round(xp)))
                elif abs(xp - np.round(xp, 1)) < 1.e-2:
                    return str(np.round(xp, 1))
                return str(np.round(xp, 2))

            cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(my_formatter))
            
        else:
            cb.set_ticks(list())
    else:
        cb = None
    return ax, cb
    
