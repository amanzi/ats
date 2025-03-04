import numpy as np
from matplotlib import pyplot as plt
import matplotlib.collections

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
        axcb = fig.colorbar(lc)

        # label it
        if colorbar_label is not None:
            axcb.set_label(colorbar_label)

        # transform ticks from [0,1] to [times[0], times[-1]]
        if colorbar_ticks:
            xticks = axcb.get_ticks()
            new_ticks = [t * (t_max - t_min) + t_min for t in xticks]
            axcb.set_ticklabels([str(np.round(t, 2)) for t in new_ticks])
        else:
            axcb.set_ticks(list())
    else:
        axcb = None
    return ax, axcb
    
