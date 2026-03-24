#!/usr/bin/env python
"""Combine ATS visualization output from restarted/continuation runs into a single
self-contained XDMF dataset.

Takes N input directories (in chronological order) and produces one output
directory containing:
  - A combined HDF5 data file with renumbered timestep keys (0, 1, 2, ...)
  - One per-step XMF per selected cycle
  - A master VisIt.xmf
  - The mesh H5 and XMF files (copied from the first run that has them)

Overlapping cycles at restart boundaries are deduplicated: for each run except
the last, any cycle whose time >= the start time of the next run is dropped.

Examples
--------
Combine three restart directories for the "domain" domain::

  combine_vis.py domain run0/ run1/ run2/ --output combined/

Same for the surface domain, with a 1-hour overlap tolerance::

  combine_vis.py surface run0/ run1/ run2/ --output combined/ --eps 3600 --time-unit s
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
# Helpers
# ---------------------------------------------------------------------------

def _validate_time_ordering(vis_files, directories, time_unit):
    """Warn loudly if runs are not in chronological order.

    The only condition we can guarantee is that start[i+1] > start[i].
    """
    for i in range(len(vis_files) - 1):
        t0 = vis_files[i].times[0]
        t1 = vis_files[i + 1].times[0]
        if t1 <= t0:
            warnings.warn(
                f"\n  *** TIME ORDERING WARNING ***\n"
                f"  Run '{directories[i+1]}' starts at t={t1:.6g} {time_unit}\n"
                f"  but run '{directories[i]}' starts at t={t0:.6g} {time_unit}.\n"
                f"  Directories may be out of order — combined output may be wrong.",
                stacklevel=2)


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
    list of (vf, directory, selected_cycle_strs, selected_times)
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
                "Check that directories are in chronological order and that --eps is appropriate.")
        else:
            run_data.append((vf, directory, sel_cycles, sel_times))

    return run_data


def _compute_selected_vars(vis_files, include_vars, exclude_vars):
    """Return the intersection of variable sets across runs, then apply filters.

    Parameters
    ----------
    vis_files : list of ats_xdmf.VisFile
    include_vars : list of str or None
    exclude_vars : list of str or None

    Returns
    -------
    list of str
    """
    var_sets = [set(vf.variables()) for vf in vis_files]
    common_vars = set.intersection(*var_sets)

    # Warn about any per-run differences
    all_vars = set.union(*var_sets)
    if common_vars != all_vars:
        for i, (vf, vs) in enumerate(zip(vis_files, var_sets)):
            missing = common_vars - vs
            if missing:
                warnings.warn(
                    f"Run {i}: missing variables present in other runs: {sorted(missing)}. "
                    "These will be excluded from the combined output.")

    # Apply user filters using the first VisFile's matching logic, then restrict
    # to common_vars.
    filtered = vis_files[0].variables(names=include_vars, exclude=exclude_vars)
    selected = [v for v in filtered if v in common_vars]

    if not selected:
        raise RuntimeError(
            "No variables selected after filtering.  Check --include / --exclude.")

    return selected


def _write_combined_h5(run_data, out_h5_path, selected_vars):
    """Write combined HDF5 with sequentially renumbered cycle keys.

    Parameters
    ----------
    run_data : list of (vf, directory, sel_cycles, sel_times)
    out_h5_path : str
    selected_vars : list of str

    Returns
    -------
    steps : list of (run_idx, old_key_str, new_key_int, h5_native_time)
        One entry per written cycle step.
    """
    steps = []
    new_key = 0

    with h5py.File(out_h5_path, 'w') as dst:
        # Propagate time unit from first run
        first_vf = run_data[0][0]
        if 'time unit' in first_vf.d.attrs:
            dst.attrs['time unit'] = first_vf.d.attrs['time unit']

        # Create variable groups up front
        for var in selected_vars:
            dst.create_group(var)

        for run_idx, (vf, directory, sel_cycles, sel_times) in enumerate(run_data):
            src = vf.d
            for old_key_str in sel_cycles:
                new_key_str = str(new_key)
                h5_native_time = None

                for var in selected_vars:
                    if var not in src:
                        warnings.warn(f"Variable '{var}' missing in '{directory}'; skipping.")
                        continue
                    if old_key_str not in src[var]:
                        warnings.warn(
                            f"Cycle key '{old_key_str}' missing for variable '{var}' "
                            f"in '{directory}'; skipping.")
                        continue

                    ds = dst[var].create_dataset(new_key_str, data=src[var][old_key_str][:])
                    if 'Time' in src[var][old_key_str].attrs:
                        t_val = src[var][old_key_str].attrs['Time']
                        ds.attrs['Time'] = t_val
                        if h5_native_time is None:
                            h5_native_time = float(t_val)

                steps.append((run_idx, old_key_str, new_key, h5_native_time))
                new_key += 1

    return steps


def _write_step_xmf(in_xmf_path, out_xmf_path, old_key_str, new_key_int,
                    h5_name, h5_native_time, selected_vars):
    """Write a per-step XMF with updated key references.

    Rewrites DataItem text from ``h5name:VARNAME/{old_key}`` to
    ``h5name:VARNAME/{new_key}``.  Removes Attribute elements for variables
    not in selected_vars.  Updates the ``<Time Value="..."/>`` element.

    Parameters
    ----------
    in_xmf_path : str
        Source per-step XMF from the input run directory.
    out_xmf_path : str
        Destination path in the output directory.
    old_key_str : str
        Original HDF5 dataset key (cycle number string, e.g. '20').
    new_key_int : int
        New sequential key.
    h5_name : str
        Data filename, e.g. 'ats_vis_data.h5'.
    h5_native_time : float or None
        Time value (in H5 native units) to write into ``<Time Value="..."/>``.
    selected_vars : list of str or None
        If None, keep all Attribute elements.
    """
    if not os.path.isfile(in_xmf_path):
        warnings.warn(
            f"Per-step XMF not found: {in_xmf_path!r}. "
            "Writing a minimal placeholder XMF.")
        _write_minimal_step_xmf(out_xmf_path, h5_name, new_key_int,
                                 h5_native_time, selected_vars)
        return

    # Register the XInclude namespace so it round-trips properly
    ET.register_namespace('xi', 'http://www.w3.org/2001/XInclude')

    tree = ET.parse(in_xmf_path)
    root = tree.getroot()

    old_suffix = '/' + old_key_str
    new_suffix = '/' + str(new_key_int)

    for grid in root.iter('Grid'):
        # Update <Time Value="..."/>
        time_elem = grid.find('Time')
        if time_elem is not None and h5_native_time is not None:
            time_elem.set('Value', f'{h5_native_time:.17e}')

        # Update Attribute elements
        to_remove = []
        for attr in list(grid):
            if attr.tag != 'Attribute':
                continue
            var_name = attr.get('Name')
            if selected_vars is not None and var_name not in selected_vars:
                to_remove.append(attr)
                continue
            di = attr.find('DataItem')
            if di is not None and di.text:
                text = di.text.strip()
                # text is like: ats_vis_data.h5:VARNAME/OLD_KEY
                # Update h5 filename to the combined file (same base name, same domain)
                parts = text.split(':', 1)
                if len(parts) == 2:
                    new_text = h5_name + ':' + parts[1]
                else:
                    new_text = text
                # Replace the trailing /OLD_KEY with /NEW_KEY
                if new_text.endswith(old_suffix):
                    new_text = new_text[:-len(old_suffix)] + new_suffix
                # Preserve leading/trailing whitespace pattern
                leading  = len(di.text) - len(di.text.lstrip())
                trailing = len(di.text) - len(di.text.rstrip())
                di.text  = di.text[:leading] + new_text + di.text[len(di.text) - trailing:]
        for elem in to_remove:
            grid.remove(elem)

    tree.write(out_xmf_path, xml_declaration=True, encoding='ASCII')


def _write_minimal_step_xmf(out_xmf_path, h5_name, new_key_int,
                              h5_native_time, selected_vars):
    """Write a skeleton per-step XMF when the source XMF is missing."""
    xdmf = ET.Element('Xdmf')
    xdmf.set('Version', '2.0')
    domain_elem = ET.SubElement(xdmf, 'Domain')
    grid = ET.SubElement(domain_elem, 'Grid')
    grid.set('GridType', 'Uniform')
    if h5_native_time is not None:
        time_elem = ET.SubElement(grid, 'Time')
        time_elem.set('Value', f'{h5_native_time:.17e}')
    tree = ET.ElementTree(xdmf)
    tree.write(out_xmf_path, xml_declaration=True, encoding='ASCII')


def _write_visit_xmf(out_directory, h5_name, n_total_steps):
    """Write the master VisIt.xmf with xi:include for each step XMF.

    Parameters
    ----------
    out_directory : str
    h5_name : str
        Data filename, e.g. 'ats_vis_data.h5'.
    n_total_steps : int
        Number of steps (0 … n_total_steps-1).
    """
    stem = h5_name[:-3]  # strip '.h5'
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

    for new_key in range(n_total_steps):
        href = f'{h5_name}.{new_key}.xmf'
        xi_elem = ET.SubElement(coll, f'{{{xi_ns}}}include')
        xi_elem.set('href', href)
        xi_elem.set('xpointer', 'xpointer(//Xdmf/Domain/Grid)')

    tree = ET.ElementTree(xdmf)
    tree.write(out_path, xml_declaration=True, encoding='ASCII')
    print(f"Wrote VisIt XMF: {out_path}")


def _copy_mesh(directories, out_directory, domain):
    """Copy mesh H5 and XMF files from the first run that has them.

    Parameters
    ----------
    directories : list of str
    out_directory : str
    domain : str
    """
    mesh_name = ats_xdmf.valid_mesh_filename(domain)
    stem = mesh_name[:-3]  # strip '.h5'

    # Find the first directory that has the mesh H5
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

    for fname in [mesh_name,
                  f'{mesh_name}.0.xmf',
                  f'{stem}.VisIt.xmf']:
        src_path = os.path.join(mesh_src_dir, fname)
        dst_path = os.path.join(out_directory, fname)
        if os.path.isfile(src_path):
            shutil.copyfile(src_path, dst_path)
            print(f"Copied mesh file: {dst_path}")
        else:
            warnings.warn(f"Mesh file not found, skipping: {src_path}")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def combineVisFiles(directories, domain, out_directory,
                    eps=1.0, time_unit='d',
                    include_vars=None, exclude_vars=None,
                    dry_run=False):
    """Combine ATS visualization output from multiple restarted runs.

    Parameters
    ----------
    directories : list of str
        Input run directories in chronological order.
    domain : str
        ATS domain name (e.g. 'domain', 'surface').
    out_directory : str
        Output directory (created if needed).
    eps : float
        Overlap tolerance in *time_unit*.  Cycles from run i with time >=
        (first_time_of_run_i+1 - eps) are dropped.  Default 1.0.
    time_unit : str
        Time unit for *eps* and for printing.  One of 's', 'hr', 'd', 'yr',
        'noleap'.  Default 'd'.
    include_vars : list of str or None
        If given, only these variables are written.
    exclude_vars : list of str or None
        If given, these variables are excluded.
    dry_run : bool
        If True, print what would be written and return without writing files.
    """
    if not directories:
        raise ValueError("Must provide at least one directory.")

    # Validate inputs
    for d in directories:
        if not os.path.isdir(d):
            raise RuntimeError(f"Input directory not found: {d!r}")

    if os.path.abspath(out_directory) in [os.path.abspath(d) for d in directories]:
        raise RuntimeError(
            f"Output directory '{out_directory}' is the same as one of the input "
            "directories.  Choose a different output path.")

    h5_name = ats_xdmf.valid_data_filename(domain)

    # Open VisFiles
    vis_files = []
    try:
        for d in directories:
            vis_files.append(
                ats_xdmf.VisFile(d, domain=domain, output_time_unit=time_unit))

        # Validate time ordering
        _validate_time_ordering(vis_files, directories, time_unit)

        # Compute selected variables (intersection + user filter)
        selected_vars = _compute_selected_vars(vis_files, include_vars, exclude_vars)

        # Deduplicate cycles across runs
        run_data = _select_cycles(vis_files, directories, eps)
        if not run_data:
            raise RuntimeError("No cycles selected across all runs.")

        total_cycles = sum(len(sel) for _, _, sel, _ in run_data)

        if dry_run:
            print(f"Combined output would contain {total_cycles} cycles "
                  f"from {len(run_data)} run(s):")
            new_key = 0
            for (vf, directory, sel_cycles, sel_times) in run_data:
                print(f"\n  Run: {directory}")
                for old_key, t in zip(sel_cycles, sel_times):
                    print(f"    new key {new_key:6d}  (old {old_key:>8})  "
                          f"t = {t:.6g} {time_unit}")
                    new_key += 1
            print(f"\n{len(selected_vars)} variables selected:")
            for v in selected_vars:
                print(f"  {v}")
            return

        # Write output
        os.makedirs(out_directory, exist_ok=True)

        # Write combined HDF5
        out_h5 = os.path.join(out_directory, h5_name)
        print(f"Writing combined HDF5: {out_h5}")
        steps = _write_combined_h5(run_data, out_h5, selected_vars)

        # Write per-step XMFs
        for (run_idx, old_key_str, new_key_int, h5_native_time) in steps:
            directory = run_data[run_idx][1]
            in_xmf  = os.path.join(directory,    f'{h5_name}.{old_key_str}.xmf')
            out_xmf = os.path.join(out_directory, f'{h5_name}.{new_key_int}.xmf')
            _write_step_xmf(in_xmf, out_xmf, old_key_str, new_key_int,
                             h5_name, h5_native_time, selected_vars)

        # Write master VisIt.xmf
        _write_visit_xmf(out_directory, h5_name, len(steps))

        # Copy mesh files
        _copy_mesh(directories, out_directory, domain)

        print(f"\nDone. Combined {len(steps)} cycles from {len(run_data)} run(s).")
        print(f"Output directory: {out_directory}")
        stem = h5_name[:-3]
        print(f"Open in VisIt:    {os.path.join(out_directory, stem + '.VisIt.xmf')}")

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
                        help='ATS domain name (e.g. "domain", "surface"), or '
                             '"*" to process all domains found in the first '
                             'input directory.')
    parser.add_argument('directories', metavar='DIRECTORY', nargs='+',
                        help='Input run directories in chronological order.  '
                             'Glob patterns (e.g. \'run*\') are accepted if quoted '
                             'and are automatically natural-sorted.  Explicit lists '
                             'are used as-is unless --sort is given.')
    parser.add_argument('--output', '-o', dest='output', required=True,
                        help='Output directory (created if needed).')
    parser.add_argument('--time-unit', dest='time_unit', default='d',
                        choices=['s', 'hr', 'd', 'yr', 'noleap'],
                        help='Time unit for --eps and printed times (default: d).')
    parser.add_argument('--eps', dest='eps', type=float, default=1.0,
                        help='Overlap tolerance in --time-unit units.  Cycles '
                             'from run i with time >= (start_of_run_i+1 - eps) '
                             'are dropped.  Default: 1.0.')

    var_group = parser.add_mutually_exclusive_group()
    var_group.add_argument('--include', dest='include_vars',
                           action='append', metavar='VAR',
                           help='Include this variable in output (repeat for multiple).')
    var_group.add_argument('--exclude', dest='exclude_vars',
                           action='append', metavar='VAR',
                           help='Exclude this variable from output (repeat for multiple).')

    parser.add_argument('--sort', dest='sort', action='store_true', default=False,
                        help='Natural-sort the directory list before combining.  '
                             'Applied automatically when a glob pattern is given.')
    parser.add_argument('--dry-run', dest='dry_run', action='store_true',
                        default=False,
                        help='Print what would be written; do not create files.')

    args = parser.parse_args()

    # Expand any glob patterns and natural-sort; trust explicit lists unless --sort
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

    if args.domain == '*':
        domains = ats_xdmf.find_domains(directories[0])
        print(f"Found domains: {domains}")
    else:
        domains = [args.domain]

    for domain in domains:
        combineVisFiles(
            directories=directories,
            domain=domain,
            out_directory=args.output,
            eps=args.eps,
            time_unit=args.time_unit,
            include_vars=args.include_vars,
            exclude_vars=args.exclude_vars,
            dry_run=args.dry_run,
        )


if __name__ == '__main__':
    main()
