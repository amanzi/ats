#!/usr/bin/env python
"""Subset ATS visualization output by time/cycle/index range and variable list.

Produces a small, self-contained directory suitable for downloading from HPC
and opening locally in VisIt.
"""

import sys
import os
import shutil
import argparse
import warnings
import xml.etree.ElementTree as ET

import numpy as np
import h5py

try:
    sys.path.insert(0, os.path.join(os.environ['ATS_SRC_DIR'], 'tools', 'utils'))
except KeyError:
    pass

import ats_xdmf


# ---------------------------------------------------------------------------
# Argument parsing helpers
# ---------------------------------------------------------------------------

def _parse_slice_or_list(s, cast):
    """Parse a slice or comma-separated list string, returning a slice or list.

    Slice syntax: 'START:STOP' or 'START:STOP:STEP' (any field may be empty)
    List syntax:  'v1,v2,v3'

    cast is applied to each non-empty value (e.g. float, int).
    """
    if ':' in s:
        parts = s.split(':')
        if len(parts) not in (2, 3):
            raise argparse.ArgumentTypeError(
                f"Cannot parse slice argument: {s!r}")
        vals = [cast(p) if p else None for p in parts]
        return slice(*vals)
    else:
        try:
            return [cast(v) for v in s.split(',')]
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Cannot parse list argument: {s!r}")


# ---------------------------------------------------------------------------
# Core subset function
# ---------------------------------------------------------------------------

def subsetVisFiles(in_directory, domain, out_directory,
                   times=None, time_unit='s', time_tolerance=1.0,
                   cycles=None, indices=None,
                   include_vars=None, exclude_vars=None,
                   dry_run=False):
    """Subset ATS visualization output files by time/cycle/index and variable.

    Parameters
    ----------
    in_directory : str
        Directory containing input visualization files.
    domain : str
        ATS domain name (e.g. 'surface', 'domain').  Used to construct
        filenames via ats_xdmf.valid_data_filename().
    out_directory : str
        Directory for output files (created if needed).  The same domain
        name is used for output filenames.
    times : slice or list or None
        Result of _parse_slice_or_list() for --times, or None.
    time_unit : str
        One of 's', 'hr', 'd', 'yr'.  Units for time filtering values.
    time_tolerance : float
        Tolerance for time matching, in time_unit.
    cycles : slice or list or None
        Result of _parse_slice_or_list() for --cycles, or None.
    indices : slice or list or None
        Result of _parse_slice_or_list() for --indices, or None.
    include_vars : list of str or None
        If given, only these variables are copied.
    exclude_vars : list of str or None
        If given, these variables are excluded.
    dry_run : bool
        If True, print what would be done and return without writing files.
    """
    h5_name   = ats_xdmf.valid_data_filename(domain)
    mesh_name = ats_xdmf.valid_mesh_filename(domain)

    in_h5 = os.path.join(in_directory, h5_name)
    if not os.path.isfile(in_h5):
        raise RuntimeError(f"No HDF5 data file found at {in_h5!r}.")

    vf = ats_xdmf.VisFile(in_directory, domain=domain,
                           output_time_unit=time_unit)

    # Apply cycle/time/index filter
    if times is not None:
        vf.filterTimes(times, time_unit=time_unit, tolerance=time_tolerance)
    elif cycles is not None:
        vf.filterCycles(cycles)
    elif indices is not None:
        vf.filterIndices(indices)
    # else: all cycles

    if not vf.cycles:
        raise RuntimeError("No cycles selected — check your filter arguments.")

    # Select variables from the h5 file
    selected_vars = vf.variables(names=include_vars, exclude=exclude_vars)

    # Dry run
    if dry_run:
        print(f"Selected {len(vf.cycles)} cycles:")
        for cycle, t in zip(vf.cycles, vf.times):
            print(f"  cycle {int(cycle):8d}  t = {t:.6g} {time_unit}")
        print(f"\nSelected {len(selected_vars)} variables:")
        for v in selected_vars:
            print(f"  {v}")
        vf.close()
        return

    # Create output directory
    os.makedirs(out_directory, exist_ok=True)

    # Write per-step XMFs for selected cycles
    for cycle in vf.cycles:
        N = int(cycle)
        in_xmf_path  = os.path.join(in_directory,  f'{h5_name}.{N}.xmf')
        out_xmf_path = os.path.join(out_directory, f'{h5_name}.{N}.xmf')

        ET.register_namespace('', 'http://www.w3.org/2001/XInclude')
        tree = ET.parse(in_xmf_path)
        root = tree.getroot()

        # Remove unselected Attribute elements
        for grid in root.iter('Grid'):
            to_remove = [
                e for e in list(grid)
                if e.tag == 'Attribute' and e.get('Name') not in selected_vars
            ]
            for e in to_remove:
                grid.remove(e)

        tree.write(out_xmf_path, xml_declaration=True, encoding='ASCII')

    # Write master VisIt.xmf
    _write_visit_xmf(out_directory, h5_name, vf.cycles)

    # Write subset H5 data file
    out_h5 = os.path.join(out_directory, h5_name)
    _write_subset_h5(in_h5, out_h5, vf.cycles, selected_vars)

    # Copy mesh H5
    in_mesh_h5  = os.path.join(in_directory,  mesh_name)
    out_mesh_h5 = os.path.join(out_directory, mesh_name)
    if os.path.isfile(in_mesh_h5):
        shutil.copyfile(in_mesh_h5, out_mesh_h5)
        print(f"Copied mesh: {out_mesh_h5}")
    else:
        warnings.warn(f"Mesh HDF5 not found: {in_mesh_h5}")

    # Copy mesh XMFs
    stem = mesh_name[:-3]  # strip '.h5'
    for suffix in [f'{mesh_name}.0.xmf', f'{stem}.VisIt.xmf']:
        in_mesh_xmf = os.path.join(in_directory, suffix)
        if os.path.isfile(in_mesh_xmf):
            out_mesh_xmf = os.path.join(out_directory, suffix)
            shutil.copyfile(in_mesh_xmf, out_mesh_xmf)
            print(f"Copied mesh XMF: {out_mesh_xmf}")

    vf.close()
    print(f"Done. Output in: {out_directory}")
    print(f"  {len(vf.cycles)} cycles, {len(selected_vars)} variables")
    print(f"  Open in VisIt: {os.path.join(out_directory, h5_name[:-3] + '.VisIt.xmf')}")


def _write_visit_xmf(out_directory, h5_name, cycles):
    """Write the master VisIt.xmf file with xi:include for each step XMF."""
    stem = h5_name[:-3]  # strip '.h5'
    out_path = os.path.join(out_directory, f'{stem}.VisIt.xmf')

    xi_ns = 'http://www.w3.org/2001/XInclude'
    xdmf = ET.Element('Xdmf')
    xdmf.set('xmlns:xi', xi_ns)
    xdmf.set('Version', '2.0')
    domain_elem = ET.SubElement(xdmf, 'Domain')
    coll = ET.SubElement(domain_elem, 'Grid')
    coll.set('GridType', 'Collection')
    coll.set('CollectionType', 'Temporal')
    coll.set('Name', 'Mesh')

    for cycle in cycles:
        href = f'{h5_name}.{int(cycle)}.xmf'
        xi_elem = ET.SubElement(coll, f'{{{xi_ns}}}include')
        xi_elem.set('href', href)

    ET.register_namespace('xi', xi_ns)
    tree = ET.ElementTree(xdmf)
    tree.write(out_path, xml_declaration=True, encoding='ASCII')
    print(f"Wrote VisIt XMF: {out_path}")


def _write_subset_h5(in_h5, out_h5, cycles, selected_vars):
    """Copy selected variables and cycles from in_h5 to out_h5."""
    with h5py.File(in_h5, 'r') as src, h5py.File(out_h5, 'w') as dst:
        # Copy the time unit attribute
        if 'time unit' in src.attrs:
            dst.attrs['time unit'] = src.attrs['time unit']
        for var in selected_vars:
            if var not in src:
                warnings.warn(f"Variable {var!r} not found in {in_h5}")
                continue
            grp = dst.create_group(var)
            for cycle in cycles:
                key = str(cycle)
                if key not in src[var]:
                    warnings.warn(
                        f"Cycle {key} not found for variable {var!r} in {in_h5}")
                    continue
                ds = grp.create_dataset(key, data=src[var][key][:])
                if 'Time' in src[var][key].attrs:
                    ds.attrs['Time'] = src[var][key].attrs['Time']
    print(f"Wrote data HDF5: {out_h5}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('domain', metavar='DOMAIN',
                        help='ATS domain name (e.g. "surface", "domain"), or '
                             '"*" to process all domains found in the input '
                             'directory.')
    parser.add_argument('-d', '--directory', dest='directory', default='.',
                        help='Directory containing input visualization files '
                             '(default: current directory)')
    parser.add_argument('--output', '-o', dest='output', default=None,
                        help='Output directory '
                             '(default: DIRECTORY/subset)')

    # Time/cycle/index selection (mutually exclusive)
    filter_group = parser.add_mutually_exclusive_group()
    filter_group.add_argument(
        '--times', dest='times', default=None,
        metavar='SLICE_OR_LIST',
        help='Time filter. Formats: '
             '"START:STOP" (range, inclusive), '
             '"START:STOP:INTERVAL" (evenly-spaced targets), '
             '"t1,t2,t3" (specific times). '
             'Units set by --time-unit.')
    filter_group.add_argument(
        '--cycles', dest='cycles', default=None,
        metavar='SLICE_OR_LIST',
        help='Cycle filter (ATS cycle numbers). Formats: '
             '"START:STOP" or "START:STOP:STEP" (inclusive both ends), '
             '"c1,c2,c3" (specific cycles).')
    filter_group.add_argument(
        '--indices', dest='indices', default=None,
        metavar='SLICE_OR_LIST',
        help='Index filter (0-based position into sorted cycle list). '
             'Formats: "START:STOP:STEP" (numpy-style, exclusive end), '
             '"i1,i2,i3". Supports negative indices.')

    parser.add_argument('--time-unit', dest='time_unit', default='s',
                        choices=['s', 'hr', 'd', 'yr', 'noleap'],
                        help='Time unit for --times values (default: s)')
    parser.add_argument('--time-tolerance', dest='time_tolerance',
                        type=float, default=1.0,
                        help='Tolerance for nearest-time matching, in '
                             '--time-unit units (default: 1.0 s)')

    # Variable selection (mutually exclusive)
    var_group = parser.add_mutually_exclusive_group()
    var_group.add_argument('--include', dest='include_vars',
                           action='append', metavar='VAR',
                           help='Include this variable in output (repeat for multiple).')
    var_group.add_argument('--exclude', dest='exclude_vars',
                           action='append', metavar='VAR',
                           help='Exclude this variable from output (repeat for multiple).')

    parser.add_argument('--dry-run', dest='dry_run', action='store_true',
                        default=False,
                        help='Print selected cycles and variables; do not '
                             'write any files.')

    args = parser.parse_args()

    in_directory = args.directory

    out_directory = args.output
    if out_directory is None:
        out_directory = os.path.join(in_directory, 'subset')

    # Parse filter arguments into slice objects or lists
    times_spec = None
    cycles_spec = None
    indices_spec = None

    if args.times is not None:
        times_spec = _parse_slice_or_list(args.times, float)
    elif args.cycles is not None:
        cycles_spec = _parse_slice_or_list(args.cycles, int)
    elif args.indices is not None:
        indices_spec = _parse_slice_or_list(args.indices, int)

    if args.domain == '*':
        domains = ats_xdmf.find_domains(in_directory)
        print(f"Found domains: {domains}")
    else:
        domains = [args.domain]

    for domain in domains:
        subsetVisFiles(
            in_directory=in_directory,
            domain=domain,
            out_directory=out_directory,
            times=times_spec,
            time_unit=args.time_unit,
            time_tolerance=args.time_tolerance,
            cycles=cycles_spec,
            indices=indices_spec,
            include_vars=args.include_vars,
            exclude_vars=args.exclude_vars,
            dry_run=args.dry_run,
        )


if __name__ == '__main__':
    main()
