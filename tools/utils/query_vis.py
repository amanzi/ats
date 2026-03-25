#!/usr/bin/env python
"""Query ATS visualization output: print a summary of variables, cycles, and time range."""

import sys
import os
import argparse

import numpy as np

try:
    sys.path.insert(0, os.path.join(os.environ['ATS_SRC_DIR'], 'tools', 'utils'))
except KeyError:
    pass

import ats_xdmf


def _display_name(varname, domain):
    name = varname
    if domain and name.startswith(domain + '-'):
        name = name[len(domain)+1:]
    if name.endswith('.cell.0'):
        name = name[:-7]
    return name


def queryVisFiles(directory, domain, filename, time_unit=None):
    """Print a summary of ATS visualization output for one file.

    Parameters
    ----------
    directory : str
    domain : str or None
        ATS domain name, or None when opening by raw filename (filepath input).
    filename : str
        Basename of the data HDF5 file.
    time_unit : str or None
        Display unit.  None uses the native unit from the file.
    """
    vf = ats_xdmf.VisFile(directory, domain=domain, filename=filename,
                           output_time_unit=time_unit)

    n_cycles = len(vf.cycles)
    first_cycle = int(vf.cycles[0]) if n_cycles > 0 else None
    last_cycle  = int(vf.cycles[-1]) if n_cycles > 0 else None

    t_first = vf.times[0]  if n_cycles > 0 else None
    t_last  = vf.times[-1] if n_cycles > 0 else None

    diffs = np.diff(vf.times) if n_cycles > 1 else np.array([])

    n_cells = None
    all_vars = vf.variables()
    if all_vars and n_cycles > 0:
        try:
            n_cells = vf.d[all_vars[0]][str(vf.cycles[0])].shape[0]
        except (KeyError, IndexError):
            pass

    label = domain if domain is not None else filename
    print(f"Domain: {label}   ({filename})")
    for key, val in vf.d.attrs.items():
        val_str = val.decode() if isinstance(val, bytes) else str(val)
        print(f"  {key}: {val_str}")
    if n_cycles > 0:
        print(f"  Cycles:    {n_cycles}  ({first_cycle} \u2026 {last_cycle})")
        print(f"  Time:      {t_first:.3f} \u2026 {t_last:.3f} {vf.output_time_unit}")
        if diffs.size > 0:
            print(f"  Frequency: min {diffs.min():.3g} {vf.output_time_unit}  "
                  f"max {diffs.max():.3g} {vf.output_time_unit}  "
                  f"mean {diffs.mean():.3g} {vf.output_time_unit}")
        else:
            print(f"  Frequency: n/a (single cycle)")
    else:
        print("  (no cycles)")

    if n_cells is not None:
        print(f"  Cells:     {n_cells}")

    sorted_vars = sorted(all_vars, key=lambda v: _display_name(v, domain).lower())
    print(f"\n  Variables ({len(sorted_vars)}):")
    for v in sorted_vars:
        disp = _display_name(v, domain)
        if disp != v:
            print(f"    {disp}  ({v})")
        else:
            print(f"    {v}")

    vf.close()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('domain', metavar='DOMAIN_OR_FILE', nargs='?', default='*',
                        help='ATS domain name, path to a vis file (.h5 or .xmf), '
                             'or "*" to summarize all domains found (default: *). '
                             'A filepath implies its own directory; -d and -p are ignored.')
    parser.add_argument('-d', '--directory', dest='directory', default='.',
                        help='Directory containing visualization files '
                             '(default: current directory). Ignored when a filepath is given.')
    parser.add_argument('-p', '--prefix', dest='prefix', default='ats_vis',
                        help='Filename prefix (default: ats_vis). Ignored when a filepath is given.')
    parser.add_argument('--time-unit', dest='time_unit', default=None,
                        choices=['s', 'hr', 'd', 'yr', 'noleap'],
                        help='Time unit for display (default: native unit from file)')

    args = parser.parse_args()

    arg = args.domain
    is_path = (os.sep in arg or '/' in arg or
               arg.endswith('.h5') or arg.endswith('.xmf'))

    if is_path:
        directory, domain, filename = ats_xdmf.resolve_vis_input(arg)
        queryVisFiles(directory, domain, filename, time_unit=args.time_unit)
    elif arg == '*':
        domains = ats_xdmf.find_domains(args.directory)
        for i, domain in enumerate(domains):
            if i > 0:
                print()
            _, _, filename = ats_xdmf.resolve_vis_input(
                domain, directory=args.directory, prefix=args.prefix)
            queryVisFiles(args.directory, domain, filename, time_unit=args.time_unit)
    else:
        _, domain, filename = ats_xdmf.resolve_vis_input(
            arg, directory=args.directory, prefix=args.prefix)
        queryVisFiles(args.directory, domain, filename, time_unit=args.time_unit)


if __name__ == '__main__':
    main()
