from matplotlib import pyplot as plt
import numpy as np


class Spline(object):
    """Forms a cublic spline on an interval given values and derivatives at the endpoints of that interval."""
    def __init__(self, x1, y1, dy1, x2, y2, dy2):
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.dy1 = dy1
        self.dy2 = dy2

    def T(self, x):
        return (x - self.x1) / (self.x2 - self.x1)

    def Value(self, x):
        t = self.T(x)
        return (1-t)**2 * ((1+2*t) * self.y1 + t * (self.x2 - self.x1) * self.dy1) \
          + t**2 * ((3-2*t) * self.y2 + (t-1) * (self.x2 - self.x1) * self.dy2)
        
    def Derivative(self, x):
        t = self.T(x)
        dtdx = 1./(self.x2 - self.x1)
        dydt = (6*t**2 - 6*t)* self.y1 \
          + (3*t**2 - 4*t + 1) * (self.x2 - self.x1) * self.dy1 \
          + (-6*t**2 + 6*t) * self.y2 \
          + (3*t**2 - 2*t) * (self.x2 - self.x1) * self.dy2
        return dydt * dtdx


class VanGenuchten(object):
    def __init__( self, alpha, n, sr, l=0.5, smoothing_interval_sat=0.0, smoothing_interval_p=0.0 ):
        self._alpha = alpha
        self._n = n
        self._sr = sr
        self._l = l
        self._m = 1 - 1.0/n 

        # smoothing for sat
        self._s0 = 1.0 - smoothing_interval_sat
        if self._s0 < 1.:
            self._spline = spline.Spline(self._s0, self.k_relative(self._s0), self.d_k_relative(self._s0),
                                         1.0, 1.0, 0.)
            
        # smoothing for pc
        self._pc0 = smoothing_interval_p
        if self._pc0 > 0.:
            self._spline_sat = spline.Spline(0., 1., 0., self._pc0, self.saturation(self._pc0), self.d_saturation(self._pc0))

    def capillaryPressure( self, s ):
        if s <= self._sr:
            return np.inf
        if s >= 1.:
            return 0.
        se = (s - self._sr) / (1.0 - self._sr)
        if (se < 1.e-8):
            return pow(se, -1.0/(self._m * self._n)) / self._alpha
        else:
            return (pow(pow(se, -1.0/self._m) - 1.0, 1/self._n)) / self._alpha

    def saturation( self, pc ):
        if pc <= 0.0:
            return 1.0
        elif pc < self._pc0:
            return self._spline_sat.Value(pc)
        else:
            se = pow(1.0 + pow(self._alpha*pc, self._n), -self._m)
            return  se * (1.0 - self._sr) + self._sr

    def k_relative( self, s ):
        if s >= 1.:
            return 1.
        elif s <= self._sr:
            return 0.
        elif s <= self._s0:
            se = (s - self._sr) / (1.0-self._sr)
            return (se**self._l) * pow( 1.0 - pow( 1.0 - pow(se,1.0/self._m),self._m), 2)
        else:
            return self._spline.Value(s)

    def label( self ):
        return "VG: a=%1.2e [1/Pa], n=%1.2g, sr=%1.2g, smooth=%g"%(self._alpha, self._n, self._sr, 1-self._s0)

    def short_label( self ):
        return "VG: a=%1.2e [1/Pa], n=%1.2g, sr=%1.2g"%(self._alpha, self._n, self._sr)
        
    
    def plot(self, ax=None, color='b', format='-', label=None, y_units='Pa'):
        pc = np.linspace(0,7, 1000)
        pc = 10**pc

        if label is None:
            label = self.short_label()

        if ax is None:
            fig,ax = plt.subplots(1,1,squeeze=True)

        s = np.array([self.saturation(apc) for apc in pc])
            
        if y_units == 'hPa':
            pc = pc / 100.
        elif y_units == 'm':
            pc = pc / 1000 / 9.81
        elif y_units == 'cm':
            pc = pc / 1000 / 9.81 * 100
        elif y_units == 'Pa':
            pass
        else:
            raise ValueError("Invalid units for yaxis, must be one of [Pa, m, cm, hPa]")
    
        ax.semilogy(s, pc, color=color, label=label)
        ax.set_xlabel("saturation [-]")
        ax.set_ylabel("capillary pressure [{}]".format(y_units))
        return ax

    def plot_kr(self, ax=None, color='b', format='-', label=None):
        if ax is None: 
            fig,ax = plt.subplots(1,1,squeeze=True)

        if label is None:
            label = self.short_label()

        pc = np.linspace(0,7, 1000)
        pc = 10**pc

        sat = np.array([self.saturation(apc) for apc in pc])
        kr = np.array([self.k_relative(s) for s in sat])
        ax.plot(sat, kr, color=color, label=label)

       
if __name__ == "__main__":
    import sys
    import argparse
    import shlex
    import colors

    parser = argparse.ArgumentParser('plot WRM curves')

    import wrm_vangenuchten
    def option_to_wrm(s):
        print("got: {}".format(s))
        try:
            s = shlex.split(s)
            print("s = {}".format(s))
            assert(len(s) == 4 or len(s) == 3)

            alpha, n, sr = map(float, s[0:3])
            print("WRM:")
            print("  alpha = {}".format(alpha))
            print("  n = {}".format(n))
            print("  sr = {}".format(sr))

            if (len(s) == 4):
                label = s[3]
            else:
                label = None
            print("  label = {}".format(label))
        except:
            raise argparse.ArgumentTypeError("WRM must be van Genucten parameters (alpha, n, sr, label)")
        else:
            return label, VanGenuchten(alpha=alpha, n=n, sr=sr)
            
    parser.add_argument('--wrm', type=option_to_wrm, action='append', help='WRM parameters, "alpha n sr [label]"')
    parser.add_argument('--y-units', type=str, choices=['Pa','m','hPa','cm'], default='Pa', help='units of the y-axis, in log space')

    args = parser.parse_args()
    color_list = colors.enumerated_colors(len(args.wrm))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for (label,wrm), color in zip(args.wrm, color_list):
        wrm.plot(ax, color, y_units=args.y_units, label=label)
    ax.legend()
    plt.show()
    sys.exit(0)
        
        
    
