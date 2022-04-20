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
            self._spline = Spline(self._s0, self.k_relative(self._s0), self.d_k_relative(self._s0),
                                         1.0, 1.0, 0.)
            
        # smoothing for pc
        self._pc0 = smoothing_interval_p
        if self._pc0 > 0.:
            self._spline_sat = Spline(0., 1., 0., self._pc0, self.saturation(self._pc0), self.d_saturation(self._pc0))

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

    def d_k_relative( self, s ):
        if s >= 1.:
            return 0
        elif s <= self._sr:
            return 0.
        elif s <= self._s0 + 1.e-6:
            se = (s - self._sr)/(1-self._sr);
            x = pow(se, 1.0 / self._m);

            if (abs(1.0 - x) < 1.e-10):
                return 0.0;
            y = pow(1.0 - x, self._m);

            dkdse = (1.0 - y) * (self._l * (1.0 - y) + 2 * x * y / (1.0 - x)) * pow(se, self._l - 1.0);
            return dkdse / (1 - self._sr);
        else:
            return self._spline.Derivative(s)

    def label( self ):
        return "VG: a=%1.2e [1/Pa], n=%1.2g, sr=%1.2g, smooth=%g"%(self._alpha, self._n, self._sr, 1-self._s0)

    def short_label( self ):
        return "VG: a=%1.2e [1/Pa], n=%1.2g, sr=%1.2g"%(self._alpha, self._n, self._sr)


class WiltingPointLimiter(object):
    def __init__( self, pc_open, pc_closed ):
        self._pc_open = pc_open
        self._pc_closed = pc_closed
        assert(pc_open < pc_closed)

    def saturation( self, pc ):
        """This is the wilting point coefficient."""
        if pc > self._pc_closed:
            return 0.
        elif pc < self._pc_open:
            return 1.
        else:
            return (self._pc_closed - pc) / (self._pc_closed - self._pc_open)

    def label( self ):
        return "WP: pc_open=%1.2e [Pa], pc_closed=%1.2e"%(self._pc_open, self._pc_closed)

    def short_label( self ):
        return self.label()
    
    
def plot(wrm, ax=None, color='b', format='-', label=None, y_units='Pa'):
    pc = np.linspace(0,7, 1000)
    pc = 10**pc

    if label is None:
        label = wrm.short_label()

    if ax is None:
        fig,ax = plt.subplots(1,1,squeeze=True)
    else:
        try:
            ax1 = ax[1]
        except TypeError:
            ax1 = None
        else:
            ax = ax[0]

    s = np.array([wrm.saturation(apc) for apc in pc])
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

    ax.semilogy(s, pc, format, color=color, label=label)
    ax.set_xlabel("saturation [-]")
    ax.set_ylabel("capillary pressure [{}]".format(y_units))

    if ax1 is not None:
        ax1.plot(s, pc/1000/9.8, format, color=color)
        ax1.set_xlabel('saturation [-]')
        ax1.set_ylabel('elevation above water table [m]')
        ax1.set_ylim([0,5])
        ax = [ax, ax1]

    return ax

def plot_kr(wrm, ax=None, color='b', format='-', label=None):
    if ax is None: 
        fig,ax = plt.subplots(1,1,squeeze=True)

    if label is None:
        label = wrm.short_label()

    pc = np.linspace(0,7, 1000)
    pc = 10**pc

    sat = np.array([wrm.saturation(apc) for apc in pc])
    kr = np.array([wrm.k_relative(s) for s in sat])
    ax.plot(sat, kr, color=color, label=label)

       
if __name__ == "__main__":
    import sys
    import argparse
    import shlex
    import colors

    parser = argparse.ArgumentParser('plot WRM curves')

    def option_to_wrm(s):
        print("got: {}".format(s))
        try:
            s = shlex.split(s)
            print("s = {}".format(s))
            assert(3 <= len(s) <= 5)

            alpha, n, sr = map(float, s[0:3])

            if len(s) > 3:
                label = s[3]
            else:
                label = None

            if len(s) > 4:
                smooth_int_sat = float(s[4])
            else:
                smooth_int_sat = 0.
                
            print("WRM:")
            print(f"  alpha = {alpha}")
            print(f"  n = {n}")
            print(f"  sr = {sr}")
            print(f"  smoothing_interval_sat = {smooth_int_sat}")
            print(f"  label = {label}")
            
        except:
            raise argparse.ArgumentTypeError("WRM must be van Genucten parameters (alpha, n, sr, label, smoothing_interval_sat)")
        else:
            return label, VanGenuchten(alpha=alpha, n=n, sr=sr, smoothing_interval_sat=smooth_int_sat)

    def option_to_wp(s):
        print("got: {}".format(s))
        try:
            s = shlex.split(s)
            print("s = {}".format(s))
            assert(2 <= len(s) <= 3)

            p_open, p_closed = map(float, s[0:2])

            if len(s) > 2:
                label = s[2]
            else:
                label = None

            print("Wilting Point limiter:")
            print(f"  p_open = {p_open}")
            print(f"  p_closed = {p_closed}")
            print(f"  label = {label}")
            
        except:
            raise argparse.ArgumentTypeError("Wilting point limiter must be parameters  (pc_open, pc_closed, label)")
        else:
            return label, WiltingPointLimiter(pc_open, pc_closed)

        
    parser.add_argument('--wrm', type=option_to_wrm, action='append', help='WRM parameters, "alpha n sr [label [smoothing_interval_sat]]"')
    parser.add_argument('--wp', type=option_to_wp, action='append', help='wilting point capillary pressure parameters, "stomata_pc_open stomata_pc_closed [label]"')
    parser.add_argument('--y-units', type=str, choices=['Pa','m','hPa','cm'], default='Pa', help='units of the y-axis, in log space')
    parser.add_argument('--kr', action='store_true', help='Plot relative permeability curve')

    args = parser.parse_args()
    if args.wrm is None:
        args.wrm = []
    if args.wp is None:
        args.wp = []
    
    color_list = colors.enumerated_colors(len(args.wrm)+len(args.wp))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if args.kr:
        for (label,wrm), color in zip(args.wrm, color_list):
            plot_kr(wrm, ax, color, label=label)
    else:
        for (label,wrm), color in zip(args.wrm, color_list):
            plot(wrm, ax, color, y_units=args.y_units, label=label)
        for (label,wp), color in zip(args.wrm, color_list[len(args.wrm):]):
            plot(wp, ax, color, y_units=args.y_units, label=label)
    ax.legend()
    plt.show()
    sys.exit(0)
        
        
    
