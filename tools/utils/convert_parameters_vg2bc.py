"""Computes Brooks-Corey water retention curve parameters given a set
of van Genuchten parameters, using the method of Lenhard et al. (1989)
or Ma et al. (1999) method 2.
"""

from argparse import ArgumentParser
from plot_wrm import *
import numpy as np
from matplotlib import pyplot as plt


def get_bc_param_from_vg(alpha, n):
    lbda = (n-1) * (1-0.5**(n/(n-1)))
    a = 0.72 - 0.35 * np.exp(-n**4)
    b = a**(1/lbda)
    c = (a**(n/(1-n))-1)**(1/n)
    p_sat = b / alpha * c
    return p_sat, 1/lbda
    

def get_bc_pc_by_sat(s, sr, p_sat, clapp_horn_b):
    if abs(1.-s) <= 1.e-6:
        return 0.
    elif s <= sr:
        return np.inf
    else:
        return p_sat * ((s - sr)/(1 - sr))**(-clapp_horn_b)
    
    
def plot_vg_and_bc(alpha, n, sr, smoothing_interval_sat=0.05):
    vg = VanGenuchten(alpha=alpha, n=n, sr=sr, smoothing_interval_sat=smoothing_interval_sat)
    sat = np.linspace(sr, 1., 100)
    kr_vg = np.array([vg.k_relative(s) for s in sat])
    pc_vg = np.array([vg.capillaryPressure(s) for s in sat])
    p_sat, clapp_horn_b = get_bc_param_from_vg(alpha, n)
    kr_bc = ((sat - sr)/(1 - sr))**(2*clapp_horn_b+3)
    pc_bc = np.array([get_bc_pc_by_sat(s, sr, p_sat, clapp_horn_b) for s in sat])
    
    fig, ax = plt.subplots(1, 2, figsize=(8, 5))
    plt.subplots_adjust(left=0.09, right=0.97, top=0.92, wspace=0.3)
    ax[0].plot(sat, kr_vg)
    ax[0].plot(sat, kr_bc)
    ax[0].set_xlabel('Saturation')
    ax[0].set_ylabel('Relative permeability')
    ax[1].plot(sat, pc_vg, label='van Genuchten, \n  alpha={:.2g}, \
               \n  n={:.2f}'.format(alpha, n))
    ax[1].plot(sat, pc_bc, label='Brooks-Corey, \n  lambda={:.2f}, \
               \n  clapp_horn_b={:.2f}, \n  p_sat={:.2f} Pa'.format(1/clapp_horn_b, clapp_horn_b, p_sat))
               
    ax[1].set_yscale('log')
    ax[1].set_xlabel('Saturation')
    ax[1].set_ylabel('Matric suction (Pa)')
    ax[1].legend();
    
    
def get_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('alpha', type=float, help='van Genuchten alpha [Pa^-1]')
    parser.add_argument('n', type=float, help='van Genuchten n')
    parser.add_argument('--sr', type=float, help='residual saturation')
    parser.add_argument('--smooth_sat', type=float, default=0, help='smoothing_interval_sat')
    parser.add_argument('--plot', action='store_true', help='plot van Genuchten vs. Brooks-Corey')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    p_sat, clapp_horn_b = get_bc_param_from_vg(args.alpha, args.n)
    print('Brooks-Corey parameters equivalent to van Genuchten: ' +
          '\nlambda={:.2f}, clapp_horn_b={:.2f}, p_sat={:.2f} Pa'\
         .format(1/clapp_horn_b, clapp_horn_b, p_sat))
    if args.plot:
        plot_vg_and_bc(args.alpha, args.n, args.sr, args.smooth_sat)
        plt.show()
