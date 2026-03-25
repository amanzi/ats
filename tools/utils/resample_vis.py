#!/usr/bin/env python
"""Subset and/or combine ATS visualization output by time/cycle/index and variable.

With a single input directory, subsets that run's visualization output.
With multiple input directories (restart/continuation runs), combines them into
a single self-contained XDMF dataset after deduplicating restart seams, and then
subsets the combined dataset.

In all cases the time/cycle/index filter is applied to the merged dataset.

Produces a self-contained output directory suitable for downloading from HPC
and opening locally in VisIt.

Examples
--------
Subset a single run by time range::

  subset_vis.py surface run0/ -o out/ --times 0:2.5 --time-unit yr

Combine three directories into a single combined result::

  subset_vis.py domain run0/ run1/ run2/ -o combined/

Combine with glob expansion::

  subset_vis.py surface 'run*/' -o combined/ --sort

Variable subset with custom output prefix::

  subset_vis.py surface run0/ -o out/ \\
      --out-prefix ats_vis_ponded_depth --include ponded_depth

Downsample to every 10th visualized step::

  subset_vis.py surface run0/ run1/ -o out/ --indices 0:10 --dry-run
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
# Internal helpers
# ---------------------------------------------------------------------------

_TIME_IN_SECONDS = {
    'yr':     365.25 * 86400,
    'noleap': 365    * 86400,
    'd':      86400,
    'hr':     3600,
    's':      1,
}


def _default_time_tolerance(effective_unit):
    """Return 0.1 s expressed in effective_unit."""
    return 0.1 / _TIME_IN_SECONDS[effective_unit]


def _validate_time_ordering(vis_files, directories, time_unit):
    """Warn loudly if runs are not in chronological order."""
    for i in range(len(vis_files) - 1):
        t0 = vis_files[i].times[0]
        t1 = vis_files[i + 1].times[0]
        if t1 <= t0:
            warnings.warn(
                f"\n  *** TIME ORDERING WARNING ***\n"
                f"  Run '{directories[i+1]}' starts at t={t1:.6g} {time_unit}\n"
                f"  but run '{directories[i]}' starts at t={t0:.6g} {time_unit}.\n"
                f"  Directories may be out of order — combined output may be wrong.",
                stacklevel=3)


def _select_cycles(vis_files, directories, eps):
    """For each run, select cycles to include after deduplication.

    Parameters
    ----------
    vis_files : list of ats_xdmf.VisFile
    directories : list of str
    eps : float
        Overlap tolerance in the VisFile's output_time_unit.  Cycles from run i
        with time >= (start_time_of_run_i+1 - eps) are dropped.

    Returns
    -------
    list of (vf, directory, selected_cycle_strs, selected_times_array)
    """
    run_data = []
    for i, (vf, directory) in enumerate(zip(vis_files, directories)):
        cycles = list(vf.cycles)
        times  = np.array(vf.times)

        if i < len(vis_files) - 1:
            next_start = vis_files[i + 1].times[0]
            cutoff = next_start - eps
            mask = times < cutoff
            sel_cycles = [c for c, m in zip(cycles, mask) if m]
            sel_times  = times[mask]
        else:
            sel_cycles = cycles
            sel_times  = times

        if len(sel_cycles) == 0:
            warnings.warn(
                f"Run '{directory}': all cycles filtered out after deduplication. "
                "Check that directories are in chronological order and that "
                "--time-tolerance is appropriate.")
        else:
            run_data.append((vf, directory, sel_cycles, sel_times))

    return run_data


def _compute_selected_vars(vis_files, include_vars, exclude_vars):
    """Return the intersection of variable sets across runs, then apply filters."""
    var_sets = [set(vf.variables()) for vf in vis_files]
    common_vars = set.intersection(*var_sets)
    all_vars    = set.union(*var_sets)

    if common_vars != all_vars:
        for vf, vs in zip(vis_files, var_sets):
            missing = common_vars - vs
            if missing:
                warnings.warn(
                    f"Directory '{vf.directory}': missing variables present in "
                    f"other runs: {sorted(missing)}. These will be excluded.")

    filtered = vis_files[0].variables(names=include_vars, exclude=exclude_vars)
    selected = [v for v in filtered if v in common_vars]

    if not selected:
        raise RuntimeError(
            "No variables selected after filtering.  Check --include / --exclude.")

    return selected


def _apply_filter(combined, times_spec, cycles_spec, indices_spec,
                  time_tolerance, effective_unit):
    """Apply --times / --cycles / --indices filter to the flat combined list.

    Parameters
    ----------
    combined : list of (vf, directory, cycle_str, time_float)
    times_spec, cycles_spec, indices_spec : slice/list/None
    time_tolerance : float
    effective_unit : str

    Returns
    -------
    filtered : list of (vf, directory, cycle_str, time_float)
    """
    if indices_spec is not None:
        if isinstance(indices_spec, slice):
            return combined[indices_spec]
        else:
            n = len(combined)
            result = []
            for idx in indices_spec:
                if -n <= idx < n:
                    result.append(combined[idx])
                else:
                    warnings.warn(f"Index {idx} out of range (0..{n-1}); skipping.")
            return result

    if cycles_spec is not None:
        if isinstance(cycles_spec, slice):
            start = cycles_spec.start
            stop  = cycles_spec.stop
            step  = cycles_spec.step
            result = []
            for entry in combined:
                c = int(entry[2])
                lo = (start is None) or (c >= start)
                hi = (stop  is None) or (c <= stop)
                result.append(entry) if (lo and hi) else None
            if step is not None and step != 1:
                # Resample by step within the already-filtered list
                result = result[::step]
            return result
        else:
            cycle_set = set(cycles_spec)
            return [e for e in combined if int(e[2]) in cycle_set]

    if times_spec is not None:
        all_times = np.array([e[3] for e in combined])
        if isinstance(times_spec, slice):
            lo = times_spec.start
            hi = times_spec.stop
            interval = times_spec.stop - times_spec.start if (
                times_spec.start is not None and times_spec.stop is not None
                and times_spec.step is not None) else None

            if times_spec.step is None:
                # Simple range
                mask = np.ones(len(combined), dtype=bool)
                if lo is not None:
                    mask &= (all_times >= lo - time_tolerance)
                if hi is not None:
                    mask &= (all_times <= hi + time_tolerance)
                return [e for e, m in zip(combined, mask) if m]
            else:
                # Evenly-spaced targets: START:STOP:INTERVAL
                targets = np.arange(lo, hi + times_spec.step * 0.5,
                                    times_spec.step)
                selected_indices = set()
                for target in targets:
                    diffs = np.abs(all_times - target)
                    best = int(np.argmin(diffs))
                    if diffs[best] <= time_tolerance:
                        selected_indices.add(best)
                return [e for i, e in enumerate(combined) if i in selected_indices]
        else:
            # Explicit list of times
            selected_indices = set()
            for target in times_spec:
                diffs = np.abs(all_times - target)
                best = int(np.argmin(diffs))
                if diffs[best] <= time_tolerance:
                    selected_indices.add(best)
                else:
                    warnings.warn(
                        f"Requested time {target} {effective_unit} not found within "
                        f"tolerance {time_tolerance} {effective_unit}; skipping.")
            return [e for i, e in enumerate(combined) if i in selected_indices]

    # No filter — return all
    return combined


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _write_output_h5(combined, out_h5_path, selected_vars):
    """Write HDF5 with original cycle keys (no renumbering)."""
    with h5py.File(out_h5_path, 'w') as dst:
        # Merge file-level attributes from all source VisFiles (superset)
        seen_vfs = {}
        for vf, _, _, _ in combined:
            if id(vf) not in seen_vfs:
                seen_vfs[id(vf)] = vf
                for key, val in vf.d.attrs.items():
                    dst.attrs[key] = val

        for var in selected_vars:
            dst.create_group(var)

        for vf, directory, cycle_str, _ in combined:
            src = vf.d
            for var in selected_vars:
                if var not in src:
                    warnings.warn(f"Variable '{var}' missing in '{directory}'; skipping.")
                    continue
                if cycle_str not in src[var]:
                    warnings.warn(
                        f"Cycle key '{cycle_str}' missing for variable '{var}' "
                        f"in '{directory}'; skipping.")
                    continue
                ds = dst[var].create_dataset(cycle_str, data=src[var][cycle_str][:])
                if 'Time' in src[var][cycle_str].attrs:
                    ds.attrs['Time'] = src[var][cycle_str].attrs['Time']

    print(f"Wrote data HDF5: {out_h5_path}")


def _write_step_xmf(in_xmf_path, out_xmf_path, cycle_str,
                    in_h5_name, out_h5_name, selected_vars):
    """Write a per-step XMF, updating h5 filename and removing unselected vars."""
    if not os.path.isfile(in_xmf_path):
        warnings.warn(f"Per-step XMF not found: {in_xmf_path!r}; skipping.")
        return

    ET.register_namespace('xi', 'http://www.w3.org/2001/XInclude')
    tree = ET.parse(in_xmf_path)
    root = tree.getroot()

    for grid in root.iter('Grid'):
        to_remove = [
            e for e in list(grid)
            if e.tag == 'Attribute' and e.get('Name') not in selected_vars
        ]
        for e in to_remove:
            grid.remove(e)

    for item in root.iter('DataItem'):
        text = item.text or ''
        if in_h5_name in text:
            item.text = text.replace(in_h5_name, out_h5_name)

    tree.write(out_xmf_path, xml_declaration=True, encoding='ASCII')


def _write_visit_xmf(out_directory, out_h5_name, cycle_strs):
    """Write the master VisIt.xmf with xi:include for each step XMF."""
    stem     = out_h5_name[:-3]  # strip '.h5'
    out_path = os.path.join(out_directory, f'{stem}.VisIt.xmf')

    xi_ns = 'http://www.w3.org/2001/XInclude'
    ET.register_namespace('xi', xi_ns)

    xdmf = ET.Element('Xdmf')
    xdmf.set('Version', '2.0')
    xdmf.set('xmlns:xi', xi_ns)
    domain_elem = ET.SubElement(xdmf, 'Domain')
    coll = ET.SubElement(domain_elem, 'Grid')
    coll.set('GridType', 'Collection')
    coll.set('CollectionType', 'Temporal')
    coll.set('Name', 'Mesh')

    for cycle_str in cycle_strs:
        href     = f'{out_h5_name}.{int(cycle_str)}.xmf'
        xi_elem  = ET.SubElement(coll, f'{{{xi_ns}}}include')
        xi_elem.set('href', href)
        xi_elem.set('xpointer', 'xpointer(//Xdmf/Domain/Grid)')

    tree = ET.ElementTree(xdmf)
    tree.write(out_path, xml_declaration=True, encoding='ASCII')
    print(f"Wrote VisIt XMF: {out_path}")


def _copy_mesh(directories, out_directory, domain):
    """Copy mesh H5 and XMF files from the first directory that has them."""
    mesh_name = ats_xdmf.valid_mesh_filename(domain)
    stem = mesh_name[:-3]  # strip '.h5'

    mesh_src_dir = None
    for d in directories:
        if os.path.isfile(os.path.join(d, mesh_name)):
            mesh_src_dir = d
            break

    if mesh_src_dir is None:
        warnings.warn(
            f"Mesh HDF5 '{mesh_name}' not found in any input directory. "
            "Output will not render correctly in VisIt.")
        return

    for fname in [mesh_name, f'{mesh_name}.0.xmf', f'{stem}.VisIt.xmf']:
        src_path = os.path.join(mesh_src_dir, fname)
        dst_path = os.path.join(out_directory, fname)
        if os.path.isfile(src_path):
            shutil.copyfile(src_path, dst_path)
            print(f"Copied mesh file: {dst_path}")
        else:
            warnings.warn(f"Mesh file not found, skipping: {src_path}")


# ---------------------------------------------------------------------------
# Core public function
# ---------------------------------------------------------------------------

def subsetVisFiles(directories, domain, out_directory, out_stem,
                   in_prefix='ats_vis',
                   times=None, time_unit=None, time_tolerance=None,
                   cycles=None, indices=None,
                   include_vars=None, exclude_vars=None,
                   dry_run=False):
    """Subset and/or combine ATS visualization output files.

    Parameters
    ----------
    directories : list of str
        Input directory or directories in chronological order.
        A single-element list is the subset case.
    domain : str
        ATS domain name (e.g. 'surface', 'domain').
    out_directory : str
        Directory for output files (created if needed).
    out_stem : str
        Output filename stem (everything before ``_data.h5``).
    in_prefix : str
        Input filename prefix.  Default ``'ats_vis'``.
    times : slice or list or None
        Time filter (result of _parse_slice_or_list with float cast), or None.
    time_unit : str or None
        Time unit for display and filtering.  None uses the native unit from file.
    time_tolerance : float or None
        Tolerance for time matching and restart deduplication, in time_unit.
        Default: 0.1 s expressed in the effective unit.
    cycles : slice or list or None
        Cycle filter, or None.
    indices : slice or list or None
        Index filter (0-based into the combined timeline), or None.
    include_vars : list of str or None
        If given, only these variables are copied.
    exclude_vars : list of str or None
        If given, these variables are excluded.
    dry_run : bool
        If True, print what would be done and return without writing files.
    """
    if not directories:
        raise ValueError("Must provide at least one directory.")

    for d in directories:
        if not os.path.isdir(d):
            raise RuntimeError(f"Input directory not found: {d!r}")

    if not dry_run:
        if os.path.abspath(out_directory) in [os.path.abspath(d) for d in directories]:
            raise RuntimeError(
                f"Output directory '{out_directory}' is the same as one of the "
                "input directories.  Choose a different output path.")

    _, _, in_filename = ats_xdmf.resolve_vis_input(
        domain, directory=directories[0], prefix=in_prefix)
    out_h5_name = f'{out_stem}_data.h5'

    vis_files = []
    try:
        # Open first file to resolve effective time unit
        vf0 = ats_xdmf.VisFile(directories[0], domain=domain,
                                filename=in_filename,
                                output_time_unit=time_unit)
        effective_unit = vf0.output_time_unit
        vis_files = [vf0] + [
            ats_xdmf.VisFile(d, domain=domain, filename=in_filename,
                              output_time_unit=effective_unit)
            for d in directories[1:]
        ]

        if time_tolerance is None:
            time_tolerance = _default_time_tolerance(effective_unit)

        if len(vis_files) > 1:
            _validate_time_ordering(vis_files, directories, effective_unit)

        selected_vars = _compute_selected_vars(vis_files, include_vars, exclude_vars)

        # Deduplicate restart seams
        run_data = _select_cycles(vis_files, directories, eps=time_tolerance)
        if not run_data:
            raise RuntimeError("No cycles selected after deduplication.")

        # Build flat combined timeline
        combined = []
        for vf, directory, sel_cycles, sel_times in run_data:
            for c, t in zip(sel_cycles, sel_times):
                combined.append((vf, directory, c, t))

        # Apply user filter
        combined = _apply_filter(combined, times, cycles, indices,
                                  time_tolerance, effective_unit)

        if not combined:
            raise RuntimeError("No cycles selected — check your filter arguments.")

        # Dry-run output
        if dry_run:
            n_dirs = len({e[1] for e in combined})
            by_dir = {}
            for vf, directory, cycle_str, t in combined:
                by_dir.setdefault(directory, []).append((cycle_str, t))

            print(f"{n_dirs} director{'y' if n_dirs == 1 else 'ies'}, "
                  f"{len(combined)} cycles selected "
                  f"(after deduplication and filtering):")
            for directory, entries in by_dir.items():
                cycles_in_dir = [e[0] for e in entries]
                times_in_dir  = [e[1] for e in entries]
                print(f"  {directory}  {len(entries)} cycles  "
                      f"cycle {int(cycles_in_dir[0])} … {int(cycles_in_dir[-1])}  "
                      f"t = {times_in_dir[0]:.3f} … {times_in_dir[-1]:.3f} "
                      f"{effective_unit}")
            print(f"\n{len(selected_vars)} variables selected:")
            for v in selected_vars:
                print(f"  {v}")
            return

        os.makedirs(out_directory, exist_ok=True)

        # Write per-step XMFs
        for vf, directory, cycle_str, _ in combined:
            N = int(cycle_str)
            in_xmf_path  = os.path.join(directory,     f'{in_filename}.{N}.xmf')
            out_xmf_path = os.path.join(out_directory, f'{out_h5_name}.{N}.xmf')
            _write_step_xmf(in_xmf_path, out_xmf_path, cycle_str,
                             in_filename, out_h5_name, selected_vars)

        # Write master VisIt.xmf
        cycle_strs = [e[2] for e in combined]
        _write_visit_xmf(out_directory, out_h5_name, cycle_strs)

        # Write HDF5 data file
        out_h5 = os.path.join(out_directory, out_h5_name)
        _write_output_h5(combined, out_h5, selected_vars)

        # Copy mesh files
        _copy_mesh(directories, out_directory, domain)

        n_dirs_used = len({e[1] for e in combined})
        print(f"\nDone. {len(combined)} cycles from "
              f"{n_dirs_used} director{'y' if n_dirs_used == 1 else 'ies'}.")
        print(f"Output directory: {out_directory}")
        print(f"Open in VisIt:    "
              f"{os.path.join(out_directory, out_stem + '.VisIt.xmf')}")

    finally:
        for vf in vis_files:
            try:
                vf.close()
            except Exception:
                pass


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('domain', metavar='DOMAIN',
                        help='ATS domain name (e.g. "surface", "domain"), or '
                             '"*" to process all domains found in the first '
                             'input directory.')
    parser.add_argument('directories', metavar='DIRECTORY', nargs='+',
                        help='One or more input directories in chronological order. '
                             'A single directory subsets that run; multiple directories '
                             'are combined after restart deduplication. '
                             'Glob patterns (e.g. \'run*/\') are accepted if quoted '
                             'and are automatically natural-sorted.')

    parser.add_argument('-p', '--prefix', dest='prefix', default='ats_vis',
                        help='Input filename prefix (default: ats_vis)')
    parser.add_argument('-o', '--output-directory', dest='output', default=None,
                        required=False,
                        help='Output directory (default: DIRECTORY/subset for a '
                             'single input, required for multiple inputs).')
    parser.add_argument('--out-prefix', dest='out_prefix', default=None,
                        help='Output filename prefix (default: same as -p). '
                             'Output stem is {out_prefix}_{domain}.')

    # Time/cycle/index selection (mutually exclusive)
    filter_group = parser.add_mutually_exclusive_group()
    filter_group.add_argument(
        '--times', dest='times', default=None,
        metavar='SLICE_OR_LIST',
        help='Time filter. Formats: '
             '"START:STOP" (range), '
             '"START:STOP:INTERVAL" (evenly-spaced targets), '
             '"t1,t2,t3" (specific times). '
             'Units set by --time-unit.')
    filter_group.add_argument(
        '--cycles', dest='cycles', default=None,
        metavar='SLICE_OR_LIST',
        help='Cycle filter. Formats: '
             '"START:STOP" or "START:STOP:STEP" (inclusive both ends), '
             '"c1,c2,c3" (specific cycles).')
    filter_group.add_argument(
        '--indices', dest='indices', default=None,
        metavar='SLICE_OR_LIST',
        help='Index filter (0-based into the combined timeline). '
             'Formats: "START:STOP:STEP" (numpy-style, exclusive end), '
             '"i1,i2,i3". Supports negative indices.')

    parser.add_argument('--time-unit', dest='time_unit', default=None,
                        choices=['s', 'hr', 'd', 'yr', 'noleap'],
                        help='Time unit for display and --times values '
                             '(default: native unit from file).')
    parser.add_argument('--time-tolerance', dest='time_tolerance',
                        type=float, default=None,
                        help='Tolerance for time matching and restart boundary '
                             'deduplication, in --time-unit units. '
                             'Default: 0.1 s expressed in the effective unit '
                             '(~1e-9 yr, ~3e-8 d, 2.8e-5 hr, 0.1 s).')

    # Variable selection (mutually exclusive)
    var_group = parser.add_mutually_exclusive_group()
    var_group.add_argument('--include', dest='include_vars',
                           action='append', metavar='VAR',
                           help='Include this variable (repeat for multiple).')
    var_group.add_argument('--exclude', dest='exclude_vars',
                           action='append', metavar='VAR',
                           help='Exclude this variable (repeat for multiple).')

    parser.add_argument('--sort', dest='sort', action='store_true', default=False,
                        help='Natural-sort the directory list. '
                             'Applied automatically when glob patterns are given.')
    parser.add_argument('--dry-run', dest='dry_run', action='store_true',
                        default=False,
                        help='Print selected cycles and variables; do not write files.')

    args = parser.parse_args()

    import glob as _glob
    directories = []
    is_glob = False
    for d in args.directories:
        if any(c in d for c in ('*', '?', '[', ']')):
            expanded = sorted(_glob.glob(d), key=ats_xdmf.natural_sort_key)
            if not expanded:
                raise RuntimeError(f"Glob pattern '{d}' matched no directories.")
            directories.extend(expanded)
            is_glob = True
        else:
            directories.append(d)

    if is_glob or args.sort:
        directories = sorted(directories, key=ats_xdmf.natural_sort_key)
        print(f"Directories (natural-sorted): {directories}")

    # Default output directory: DIRECTORY/subset for single input
    if args.output is None:
        if len(directories) > 1:
            parser.error(
                "argument -o/--output-directory is required when multiple "
                "input directories are given.")
        out_directory = os.path.join(directories[0], 'subset')
    else:
        out_directory = args.output

    out_prefix = args.out_prefix or args.prefix

    # Parse filter arguments
    times_spec = cycles_spec = indices_spec = None
    if args.times is not None:
        times_spec = _parse_slice_or_list(args.times, float)
    elif args.cycles is not None:
        cycles_spec = _parse_slice_or_list(args.cycles, int)
    elif args.indices is not None:
        indices_spec = _parse_slice_or_list(args.indices, int)

    if args.domain == '*':
        domains = ats_xdmf.find_domains(directories[0])
        print(f"Found domains: {domains}")
    else:
        domains = [args.domain]

    for domain in domains:
        out_stem = out_prefix if domain == 'domain' else f'{out_prefix}_{domain}'
        subsetVisFiles(
            directories=directories,
            domain=domain,
            out_directory=out_directory,
            out_stem=out_stem,
            in_prefix=args.prefix,
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
