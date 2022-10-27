import h5py
import numpy as np
from argparse import ArgumentParser

def rh2vp(filename_rh_h5, filename_vp_h5):
    data = {}
    attrs = {}
    with h5py.File(filename_rh_h5, 'r') as f:
        for k, v in f.items():
            data[k] = v[:]
        for k, v in f.attrs.items():
            attrs[k] = v
            
    vp_sat = 611.2*np.exp(17.67*(data['air temperature [K]']-273.15)/(data['air temperature [K]']-273.15+243.5))
    data['vapor pressure air [Pa]'] = data['relative humidity [-]']*vp_sat
    data.pop('relative humidity [-]')
    
    with h5py.File(filename_vp_h5, 'w') as f:
        for k, v in data.items():
            f.create_dataset(k, data=v)
        for k, v in attrs.items():
            f.attrs[k] = v
            
if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('filename_rh_h5', type=str, 
                        help='filename of ats forcing.h5 using relative humidity')
    parser.add_argument('-f', '--filename_vp_h5', type=str,
                        help='filename of ats forcing.h5 using vapor pressure')
    parser.add_argument('--inplace', action='store_true', help='Convert rh to vp inplace')
    
    args =  parser.parse_args()
    
    if args.inplace:
        rh2vp(args.filename_rh_h5, args.filename_rh_h5)
    else:
        rh2vp(args.filename_rh_h5, args.filename_vp_h5)
