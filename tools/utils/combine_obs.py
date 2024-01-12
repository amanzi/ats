#!/usr/bin/env python
# -------------------------------------------------------------
# file: combine_xmf.py
#
# Combines xmf files from restarted runs into one VisIt readable file.
#

import sys,os
import pandas as pd


def loadFile(filename):
    print(f'Loading {filename}')
    return pd.read_csv(filename, comment='#')

def readAndCombine(input_files, eps=1.e-5):
    dat = [loadFile(f) for f in input_files]
    dat_merged = []

    time_col = dat[0].columns[0]
    for i in range(len(dat)-1):
        new_time_start = dat[i+1][time_col][0]
        print(f'Truncating: {input_files[i]}')
        print(f'  new time start = {new_time_start}')
        dat_new = dat[i][dat[i][time_col] < new_time_start-eps]
        print(f'  old extent = {dat[i][time_col][0]}, {dat[i][time_col][len(dat[i])-1]}')
        print(f'  new extent = {dat_new[time_col][0]}, {dat_new[time_col][len(dat_new)-1]}')
        dat_merged.append(dat_new)

    print(f'Not Truncating: {input_files[-1]}')
    print(f'  extent = {dat[-1][time_col][0]}, {dat[-1][time_col][len(dat[-1])-1]}')
    dat_merged.append(dat[-1])
        
    merged = pd.concat(dat_merged, ignore_index=True)
    return merged


def readHeader(input_file):
    with open(os.path.join(input_file), 'r') as fid:
        line = fid.readline()
        header = []
        while line.strip().startswith("#"):
            header.append(line)
            line = fid.readline()
    return ''.join(header)


def write(output_file, header, data):
    column_name_header = ','.join([f'"{n}"' for n in list(data.columns)])

    with open(output_file, 'w') as fid:
        fid.write(header)
        fid.write(column_name_header+'\n')

    data.to_csv(output_file, index=False, header=None, mode='a',
                float_format='%1.8e', chunksize=1000)
    


epilog = \
"""This can be used in two ways.  The simpler is to provide a list of
FILES, and an output filename.  The more convenient version when
running a series of ATS restarts in different directories is to
provide a list of DIRECTORIES, the name of the observation file
(e.g. water_balance.dat) created by each run in those directories, and
a new directory as output that will be used a meta-run collecting
output from all runs.

Examples:

$> combine_obs.py -o combined.dat run1.dat run2.dat run3.dat
$> combine_obs.py -o combined/ -f water_balance.dat run0/ run1/ run2/

"""

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Combine observation files from a sequence of restarted runs into a single CSV.",
                                     epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("FILES_OR_DIRECTORIES", nargs="+", type=str,
                        help="List of directories to combine.")

    parser.add_argument("--filename", "-f", type=str,
                        help="File name to be combined, if using directories.")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Filename or directory in which to place the merged file.")

    args = parser.parse_args()
    if len(args.FILES_OR_DIRECTORIES) < 1:
        raise RuntimeError("Specify nonzero length list of directories.")

    # are we in directory or file mode?  Parse to generate list of input files and output file
    if os.path.isfile(args.FILES_OR_DIRECTORIES[0]):
        input_files = args.FILES_OR_DIRECTORIES
        output_file = args.output
        if output_file is None:
            raise RuntimeError("If using with FILES, must provide the --output filename")
    else:
        if args.filename is None:
            raise RuntimeError("If using with DIRECTORIES, must provide the --filename")

        input_files = [os.path.join(dirname, args.filename) for dirname in args.FILES_OR_DIRECTORIES]

        output_dir = args.output
        if output_dir is None:
            output_dir = '.'
        output_file = os.path.join(output_dir, args.filename)

    # get the header
    header = readHeader(input_files[0])

    # create the combined dataframe
    data = readAndCombine(input_files)

    # write the file
    write(output_file, header, data)





