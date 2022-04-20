#!/usr/bin/env python
# -------------------------------------------------------------
# file: combine_xmf.py
#
# Combines xmf files from restarted runs into one VisIt readable file.
#

import sys,os
import pandas as pd


def loadFile(filename, directory):
    ff = os.path.join(directory, filename)
    print(f'Loading {ff}')
    return pd.read_csv(ff, comment='#')

def combine(filename, directory_list, eps=1.e-5):
    dat = [loadFile(filename, directory) for directory in directory_list]
    dat_merged = []
    time_col = dat[0].columns[0]
    for i in range(len(dat)-1):
        new_time_start = dat[i+1][time_col][0]
        print(f'Truncating: {directory_list[i]}')
        print(f'  new time start = {new_time_start}')
        dat_new = dat[i][dat[i][time_col] < new_time_start-eps]
        print(f'  old extent = {dat[i][time_col][0]}, {dat[i][time_col][len(dat[i])-1]}')
        print(f'  new extent = {dat_new[time_col][0]}, {dat_new[time_col][len(dat_new)-1]}')
        dat_merged.append(dat_new)
        
    print(f'Not Truncating: {directory_list[-1]}')
    print(f'  extent = {dat[-1][time_col][0]}, {dat[-1][time_col][len(dat[-1])-1]}')
    dat_merged.append(dat[-1])

    merged = pd.concat(dat_merged, ignore_index=True)
    return merged

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Combine *.VisIt.xmf files into a single VisIt readable file.")
    parser.add_argument("DIRECTORIES", nargs="+", type=str,
                        help="List of directories to combine.")

    parser.add_argument("--filename", "-f", type=str,
                        help="File name to be combined")
    parser.add_argument("--output-directory", "-o", type=str, default='.',
                        help="Directory in which to place the merged file.")

    args = parser.parse_args()
    if len(args.DIRECTORIES) < 1:
        parser.print_help()
        raise RuntimeError("Specify nonzero length list of directories.")

    # create the directory objects
    data = combine(args.filename, args.DIRECTORIES)

    # get the header
    with open(os.path.join(args.DIRECTORIES[0],args.filename), 'r') as fid:
        line = fid.readline()
        header = []
        while line.strip().startswith("#"):
            header.append(line)
            line = fid.readline()
    headerlines = ''.join(header)

    # pandas seems to be failing on writes to file, write by line
    column_name_header = ','.join([f'"{n}"' for n in list(data.columns)])

    # write the file
    outfile = os.path.join(args.output_directory, args.filename)
    with open(outfile, 'w') as fid:
        fid.write(headerlines)
        fid.write(column_name_header+'\n')

    data.to_csv(outfile, index=False, header=None, mode='a', float_format='%1.8e', chunksize=1000, line_terminator='\n')





