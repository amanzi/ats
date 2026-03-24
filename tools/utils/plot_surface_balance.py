#!/usr/bin/env python
"""
Plot met data from an ATS input h5 file using default names.

This is currently only useful on a single, 1D column.
"""

import os,sys
import h5py
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm
import ats_xdmf
import itertools
import colors

def _get_surface_data(keys, dat, name):
    """Replacement for parse_ats.getSurfaceData using raw h5py handle."""
    if len(name.split(".")) != 3:
        name = name + ".cell.0"
    res = np.array([dat[name][key][0] for key in keys])
    if len(res.shape) == 2 and res.shape[1] == 1:
        res = res[:,0]
    return res


def get_filename_base(base, suffix):
    if suffix is not None:
        fname = "visdump_{0}_{1}_data.h5".format(base, suffix)
    else:
        fname = "visdump_{0}_data.h5".format(base)
    return fname

def get_filename(directory, base, suffix):
    return os.path.join(directory, get_filename_base(base, suffix))


def plot_surface((keys,times,dat), ax, color, varname, style='-', version='dev-new', domain_suffix=None,
                 negate=False, label=None, hackfactor=None, divide_by_cell_volume=False):
    if version == '0.86':
        if varname == 'pressure' or varname == 'temperature':
            if domain_suffix is None:
                prefix = "surface-"
            else:
                prefix = "surface_{0}-".format(domain_suffix)
        else:
            prefix = ""
    elif version == 'dev':
        if domain_suffix is None:
            prefix = "surface-"
        else:
            prefix = "surface_{0}-".format(domain_suffix)
    elif version == 'dev-new':
        if varname.startswith('snow'):
            varname = varname[5:]
            if domain_suffix is None:
                prefix = 'snow-'
            else:
                prefix = 'snow_{0}-'.format(domain_suffix)
        else:
            if domain_suffix is None:
                prefix = 'surface-'
            else:
                prefix = 'surface_{0}-'.format(domain_suffix)
    else:
        raise RuntimeError("Unrecognized version '%s'"%version)
    sd = _get_surface_data(keys, dat, prefix+varname)

    if divide_by_cell_volume:
        if domain_suffix is None:
            cv_key = 'surface-cell_volume'
        else:
            cv_key = 'surface_{0}-cell_volume'.format(domain_suffix)
        cv = _get_surface_data(keys, dat, cv_key)
        sd = sd / cv
    
    if negate:
        sd = -sd

    if sd[0] == 0.:
        sd[0] = np.nan

    if hackfactor is not None:
        sd = sd * hackfactor

    ax.plot(times*365.25, sd, style, color=color, label=label, markerfacecolor='none')
    ax.set_ylabel(varname)
    ax.set_xlabel("time [d]")


def plot_surface_balance(ktd, ax, color, marker=None, version='dev-new', label=None, hackfactor=None, domain_suffix=None):
    """Plot data provided as a (keys,times,dat) tuple on a set of axes."""
    if marker is not None:
        style = '-'+marker
        style2 = '--'+marker
    else:
        style = '-'
        style2 = '--'

    ax = ax.ravel()

    if ktd[2][1] is None:
        ktds = (ktd[0], ktd[1], ktd[2][0])
    else:
        ktds = (ktd[0], ktd[1], ktd[2][1])
    ktd = (ktd[0], ktd[1], ktd[2][0])

    plot_surface(ktd,ax[0], color, 'incoming_shortwave_radiation', style=style, version=version, label=label, domain_suffix=domain_suffix)
    try:
        plot_surface(ktd, ax[1], color, 'incoming_longwave_radiation', style=style, version=version, domain_suffix=domain_suffix)
    except KeyError:
        print "Cannot find longwave radiation (old runs may not have this).  Skipping..."

    if version == '0.86':
        plot_surface(ktd, ax[2], color, 'qE_lw_out', style=style, version=version, negate=True, domain_suffix=domain_suffix)
    else:
        plot_surface(ktd, ax[2], color, 'qE_lw_out', style=style, version=version, domain_suffix=domain_suffix)
    #plot_surface(ktd, ax[3], color, 'qE_latent_heat', style=style, version=version, domain_suffix=domain_suffix)
    plot_surface(ktd, ax[3], color, 'evaporative_flux', style=style, version=version, domain_suffix=domain_suffix)
    plot_surface(ktd, ax[4], color, 'qE_sensible_heat', style=style, version=version, domain_suffix=domain_suffix)

    if version == '0.86':
        plot_surface(ktd, ax[5], color, 'total_energy_source', style=style, hackfactor=hackfactor, divide_by_cell_volume=True, domain_suffix=domain_suffix)
    else:
        plot_surface(ktd, ax[5], color, 'total_energy_source', style=style, hackfactor=hackfactor, domain_suffix=domain_suffix)
        
    plot_surface(ktds, ax[6][0], color, 'snow_temperature', style=style, version=version, domain_suffix=domain_suffix)
    plot_surface(ktds, ax[6][1], colors.darken(color), 'snow_depth', style=style2, version=version, domain_suffix=domain_suffix)

    plot_surface(ktd, ax[7][0], color, 'temperature', style=style, version=version, domain_suffix=domain_suffix)
    plot_surface(ktd, ax[7][1], colors.darken(color), 'ponded_depth', style=style2, version=version, domain_suffix=domain_suffix)

    plot_surface(ktd, ax[8], color, 'unfrozen_fraction', style=style, version=version, domain_suffix=domain_suffix)


def get_version(dirname, domain_suffix):
    if os.path.isfile(get_filename(dirname, 'snow', domain_suffix)):
        version = 'dev-new'
    elif os.path.isfile(get_filename(dirname, 'snow', None)):
        version = 'dev-new'
        domain_suffix = None
    else:
        d = h5py.File(get_filename(dirname, 'surface', domain_suffix), 'r')
        if domain_suffix is None:
            precip_key = 'surface-precipitation_rain.cell.0'
        else:
            precip_key = 'surface_{0}-precipitation_rain.cell.0'.format(domain_suffix)
        if precip_key in d.keys():
            version = 'dev'
        else:
            version = '0.86'
        d.close()
    return version, domain_suffix

    
def get_axs():
    f,axs = plt.subplots(3,3, figsize=(14.4,8.3))
    axs = axs.ravel()
    axs[6] = (axs[6], axs[6].twinx())
    axs[7] = (axs[7], axs[7].twinx())
    return f,axs
        
if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("INPUT_DIRECTORIES", nargs="+", type=str,
                        help="List of directories to plot.")
    parser.add_argument("--colors", "-c", type=colors.float_list_type,
                        default=None,
                        help="List of color indices to use, of the form: --colors=[0,0.1,1], where the doubles are in the range (0,1) and are mapped to a color using the colormap.")
    parser.add_argument("--colormap", "-m", type=str,
                        default="jet",
                        help="Colormap used to pick a color.")
    parser.add_argument("--domain_suffix", type=str,
                        help="Suffix of surface and snow domains.")
    args = parser.parse_args()

    cm = colors.cm_mapper(0,1,args.colormap)
    fig, axs = get_axs()

    dirnames = args.INPUT_DIRECTORIES
    markers = reversed(['x', '+', 'o', '^', 'v', 's'][0:len(dirnames)])
    for i,(dirname,marker) in enumerate(zip(dirnames, itertools.cycle(markers))):
        version, domain_suffix = get_version(dirname, args.domain_suffix)
        print("Directory %s got version: %s and domain_suffix: %r"%(dirname,version, domain_suffix))
        hackfactor = None
        if dirname == "run-transient4":
            hackfactor = 2.0
            
        if args.colors is None:
            if len(dirnames) > 1:
                c = cm(float(i)/(len(dirnames)-1))
            else:
                c = 'b'
        else:
            if type(args.colors[i]) is float:
                c = cm(args.colors[i])
            else:
                c = args.colors[i]

        vf_surf = ats_xdmf.VisFile(dirname, filename=get_filename_base('surface', domain_suffix))
        ktd = (vf_surf.cycles, vf_surf.times, vf_surf.d)
        if version == 'dev-new':
            vf_snow = ats_xdmf.VisFile(dirname, filename=get_filename_base('snow', domain_suffix))
            ktds = vf_snow.d
        else:
            vf_snow = None
            ktds = None
        plot_surface_balance((ktd[0], ktd[1], (ktd[2], ktds)), axs, c, label=dirname, hackfactor=hackfactor, version=version, marker=marker, domain_suffix=domain_suffix)
        if vf_snow is not None:
            vf_snow.close()
        vf_surf.close()

    plt.tight_layout()
    axs[0].legend(bbox_to_anchor=(0., 1., 3.6, .05), loc=3,
               ncol=len(dirnames), mode="expand", borderaxespad=0.)
    plt.show()
