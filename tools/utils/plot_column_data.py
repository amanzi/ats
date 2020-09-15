#/usr/bin/env python
"""Functions for plotting column runs of ATS."""
import sys,os
from matplotlib import pyplot as plt
import matplotlib.cm

import ats_xdmf, ats_units




epilog = """
Note this plots:
* subsurface quantities (var, z) with a line for each time
* surface quantities (t, var) with a single line

Layout
------

One way of specifying what to plot is to use the layout (-x) 
option.  This option takes a potentially 2D array of lists of variables
to plot in a 2D array of axes on the figure.  The best way to understand
how to use this option is via examples:

 --layout="[pressure, surface-pressure],[temperature, surface-temperature]"

   Makes a figure with 4 axes in a 2x2 block.  Note each list is a row.

 --layout="saturation_liquid + saturation_ice, temperature"

   Makes a figure with 2 axes, the first containing saturations and
   the second temperature.

 --layout=arctic

   Makes a plot for Arctic column runs.

Also, more default layouts might be reasonable to add.

Slicing
-------
Slicing in time (only for subsurface currently).  Valid choices for the
--time-slice option are:

 * A single int -- this specifies a single cycle
 * A list of comma-separated ints (no spaces) -- a list of cycles
 * A list of up to three colon-separated ints (no spaces) -- a numpy-style
   slice of cycles.

Note this also tries to guess reasonable values for axis extents.
There is no way of overriding the default values right now other than
to edit the code -- this might be an area for future development.
"""

def valid_layout(arg):
    if arg == 'arctic':
        arg = '[[saturation_liquid+saturation_ice, surface-ponded_depth],'+ \
               '[temperature, snow-depth]]'

    # whitespace never changes things -- remove it all
    arg = arg.replace(' ','')
    if not arg.startswith('[['):
        arg = '['+arg
    if not arg.startswith('[['):
        arg = '['+arg

    if not arg.endswith(']]'):
        arg = arg+']'
    if not arg.endswith(']]'):
        arg = arg+']'

    layout = list()
    assert(arg.startswith('['))
    assert(arg.endswith(']'))
    arg = arg[1:-1]
    rows = arg.split(']')
    assert(len(rows) > 1)

    empty_last = rows.pop()
    assert(empty_last == '' or empty_last == ',')
    for row in rows:
        if row.startswith(',['): # second or later row
            row = row[2:]
        elif row.startswith('['): # first row
            row = row[1:]
        else:
            assert(False)

        row_layout = list()
        columns = row.split(',')
        
        if columns[-1] == '':
            columns.pop()
        for column in columns:
            variables = column.split('+')
            row_layout.append(variables)
        layout.append(row_layout)

    n_cols = max(len(row) for row in layout)
    for row in layout:
        if len(row) != n_cols:
            row.append(list())
    return layout    


def domain_var(varname):
    vs = varname.split('-')
    if len(vs) == 2:
        return vs[0],vs[1]
    elif len(vs) == 1:
        return '', varname
    else:
        raise RuntimeError("Invalid variable name: {}".format(varname))

def label(varname):
    return '{} [{}]'.format(varname, ats_units.get_units(varname))
    
def annotate_surface(col, ax, time_unit):
    ax.set_ylabel(','.join([label(c) for c in col]))
    ax.set_xlabel('time [{}]'.format(time_unit))

def annotate_subsurface(col, ax):
    ax.set_ylabel('z [m]')
    ax.set_xlabel(','.join([label(c) for c in col]))
    
def annotate(layout, axs, time_unit):
    axs[0][0].legend()
    for i,row in enumerate(layout):
        for j,col in enumerate(row):
            if col[0].startswith('surface') or col[0].startswith('snow'):
                annotate_surface(col, axs[i][j], time_unit)
            else:
                annotate_subsurface(col, axs[i][j])

def plot_subsurface(vis, col, ax, label, color=None, cmap=None):
    if cmap is None:
        cmap = colors.alpha_cmap(color)
    
    z = vis.centroids[:,2]

    if len(vis.times) == 1:
        cm = colors.cm_mapper(vis.times[0]-1, vis.times[0], cmap)
    else:
        cm = colors.cm_mapper(vis.times[0], vis.times[-1], cmap)

    formats = ['-', '--', '-.']

    for varname, form in zip(col, formats):
        data = vis.getArray(varname)
        assert(len(data.shape) == 2)
        assert(data.shape[1] == len(vis.centroids))
        
        for i,time in enumerate(vis.times):
            mylabel = None
            if i == len(vis.times)-1:
                mylabel = label

            ax.plot(data[i,:], z, form, color=cm(time), label=mylabel)

def plot_surface(vis, col, ax, label, color):
    formats = ['-', '--', '-.']
    for varname, form in zip(col, formats):
        data = vis.getArray(varname)
        assert(data.shape[1] == 1)
        ax.plot(vis.times, data[:,0], form, color=color, label=label)

def plot(vis_objs, layout, axs, dirname, color):
    for i,row in enumerate(layout):
        for j,col in enumerate(row):
            ax = axs[i][j]
            domain, _ = domain_var(col[0])
            vis = vis_objs[domain]
            if domain == '':
                plot_subsurface(vis, col, ax, dirname, color)
            else:
                plot_surface(vis, col, ax, dirname, color)

def layout_flattener(layout):
    for row in layout:
        for col in row:
            for v in col:
                yield v

def directory(d):
    """Argparse validator for a directory."""
    if not os.path.isdir(d):
        raise argparse.ArgumentTypeError('Nonexistent directory "{}"'.format(d))
    return d    

def time_slice(sl):
    """Argparse validator for slices."""
    if ',' in sl:
        return map(int, sl.split(','))
    elif ':' in sl:
        def s_to_i(s):
            if s == '':
                return None
            else:
                return int(s)
            
        return slice(*map(s_to_i, sl.split(':')))
    else:
        return int(sl)
                

if __name__ == '__main__':
    import argparse
    import colors

    parser = argparse.ArgumentParser("Plot column runs of ATS", epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('directories', nargs='+', type=directory,
                        help='Directories of runs to plot.')
    parser.add_argument('--data-filename-format', type=str,
                        default='visdump_{}_data.h5',
                        help='Format string to generate a data hdf5 filename.')
    parser.add_argument('--mesh-filename-format', type=str,
                        default='visdump_{}_mesh.h5',
                        help='Format string to generate a mesh hdf5 filename.')
    parser.add_argument('--time-unit', type=str, default='d',
                        help='Unit for time plots.')

    parser.add_argument('--subsurface-time-slice', type=time_slice, 
                        help='Time slice information, see below')
    parser.add_argument('--surface-time-slice', type=time_slice, 
                        help='Time slice information, see below')
    
    parser.add_argument('-x', '--layout', type=valid_layout,
                        help='Layout in a complicated format -- see below.')
    parser.add_argument('-f', '--figsize', type=float, nargs=2, default=[5,3],
                        help='Figure size')
    parser.add_argument('-c', '--colors', type=str,
                        choices=['enumerated','sampled'], default='enumerated',
                        help='Choose colors from an enumerate set or a sampled colormap.')

    
    args = parser.parse_args()

    if args.colors == 'enumerated':
        color_list = colors.enumerated_colors(len(args.directories))
    else:
        color_list = colors.sampled_colors(len(args.directories), matplotlib.cm.jet)

    fig = plt.figure(figsize=args.figsize)
    axs = fig.subplots(len(args.layout), len(args.layout[0]), squeeze=False)

    domains = set([domain_var(v)[0] for v in layout_flattener(args.layout)])
    for dirname, color in zip(args.directories, color_list):
        vis_objs = dict()
        for domain in domains:
            vis = ats_xdmf.VisFile(dirname, domain,
                                   ats_xdmf.valid_data_filename(domain, args.data_filename_format),
                                   ats_xdmf.valid_mesh_filename(domain, args.mesh_filename_format),
                                   time_unit=args.time_unit)
            if domain == '':
                vis.loadMesh(columnar=True)
                if args.subsurface_time_slice is not None:
                    vis.filterIndices(args.subsurface_time_slice)
            else:
                vis.loadMesh()
                if args.surface_time_slice is not None:    
                    vis.filterIndices(args.surface_time_slice)

            vis_objs[domain] = vis

        plot(vis_objs, args.layout, axs, dirname, color)

    annotate(args.layout, axs, args.time_unit)        
    plt.show()
    sys.exit(0)
    
    
