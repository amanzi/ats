from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize
import pytest

class AreaFraction:
    def __init__(self, del_max, del_ex, snow_transition=0.02, rho=1000.):
        self.del_max = del_max
        self.del_ex = del_ex
        self.snow_transition = snow_transition
        self.rho = rho

    def f(self, delta):
        if (delta >= self.del_max):
            return delta - self.del_ex
        else:
            d_on_dmax = delta / self.del_max
            return pow(d_on_dmax, 2) * (2*self.del_max - 3*self.del_ex) \
                + pow(d_on_dmax, 3) * (2*self.del_ex - self.del_max)
    def f_prime(self, delta):
        if (delta >= self.del_max):
            return 1
        else:
            d_on_dmax = delta / self.del_max
            return 2 * d_on_dmax * (2*self.del_max - 3*self.del_ex) / self.del_max \
                + 3 * pow(d_on_dmax, 2) * (2*self.del_ex - self.del_max) / self.del_max

    def f_inv(self, phi):
        def fun(delta):
            return self.f(delta) - phi
        
        return scipy.optimize.brentq(fun, 0, self.del_max)

    def area_fracs_analytic2(self, pressure, snow_height):
        p_atm = 101325.
        rhogz = self.rho * 9.81

        if pressure > p_atm:
            pd = (pressure - p_atm) / rhogz
        else:
            pd = 0.
        sd = snow_height

        # determine the land
        full = self.del_max - self.del_ex
        if pd + sd >= self.del_max:
            area_frac_land = 0
        else:
            area_frac_land = 1 - self.f_prime(pd + sd)

        # now partition the snow
        potential_snow = 1 - area_frac_land
        vpd = self.f(pd)
        vsd = self.f(pd+sd) - vpd
        if (vsd >= potential_snow * self.snow_transition):
            area_frac_snow = potential_snow
            area_frac_water = 0.
        else:
            area_frac_snow = vsd / self.snow_transition
            area_frac_water = 1-area_frac_snow-area_frac_land

        return area_frac_land, area_frac_water, area_frac_snow
            
        
    def area_fracs_analytic(self, pressure, snow_depth):
        p_atm = 101325.
        rhogz = self.rho * 9.81
        pd = (pressure - p_atm) / rhogz

        if pd <= 0.:
            if snow_depth <= 0.:
                return 1,0,0
            elif snow_depth >= self.del_max - self.del_ex:
                return 0,0,1
            else:
                del_snow = self.f_inv(snow_depth)
                snow_frac = self.f_prime(del_snow)
                return 1-snow_frac, 0, snow_frac
        elif snow_depth <= 0.:
            water_frac = self.f_prime(pd)
            return 1-water_frac, water_frac, 0
        
        # water
        vpd = 0.
        if (pd <= 0.):
            vpd = 0.
            water_frac = 0.
        elif pd < self.del_max:
            vpd = self.f(pd)
            water_frac = self.f_prime(pd)
        else:
            vpd = pd - self.del_ex
            water_frac = 1.
        init_water_frac = water_frac

        #land
        if vpd + snow_depth >= self.del_max - self.del_ex:
            land_frac = 0.
            total_depth = self.del_max - self.del_ex + (vpd + snow_depth - (self.del_max-self.del_ex))
        else:
            total_depth = self.f_inv(vpd + snow_depth)
            land_frac = 1.0 - self.f_prime(total_depth)

        # snow
        if snow_depth <= self.snow_transition * (1 - land_frac):
            snow_frac = snow_depth / self.snow_transition
            water_frac = 1 - land_frac - snow_frac
        else:
            snow_frac = 1 - land_frac
            water_frac = 0

        try:
            assert(-1.e-10 <= land_frac <= 1+1.e-10)
            assert(-1.e-10 <= water_frac <= 1+1.e-10)
            assert(-1.e-10 <= snow_frac <= 1+1.e-10)
            assert(abs(1 - (land_frac + snow_frac + water_frac)) < 1.e-10)
        except AssertionError as err:
            print("broken at:")
            print('sd, pd = ',snow_depth, pd)
            print('dm, de = ', self.del_max, self.del_ex)
            print('vpd =',vpd)
            print('init wa =', init_water_frac)
            print('total d =', total_depth)
            print(land_frac, water_frac, snow_frac)
            raise err
        return land_frac, water_frac, snow_frac


            
            
            
        
    def area_fracs(self, pressure, snow_depth):
        p_atm = 101325.
        rhogz = self.rho * 9.81
        pd = (pressure - p_atm) / rhogz

        if (abs(snow_depth - 0.007) < .0011 and abs(pd - 0.04) < 1.e-5):
            debug = True
            print('----')
            print("debugging: sd, pd = ", snow_depth, pd)
        else:
            debug = False

        # water
        vpd = 0.
        if (pd <= 0.):
            water_frac = 0.
            vpd = 0.
        elif pd < self.del_max:
            vpd = self.f(pd)
            water_frac = self.f_prime(pd)
        else:
            water_frac = 1.
            vpd = pd - self.del_ex

        if debug:
            print("vpd = ", vpd)
            print("water frac = ", water_frac)

        #snow
        # volume available for snow before it covers everything is self.del_max - self.del_ex - volumetric pd
        fully_covered = max(self.del_max - self.del_ex - vpd, 0.)
        if debug:
            print("fully covered = ", self.del_max - self.del_ex - vpd)
            
        if (snow_depth < water_frac*self.snow_transition):
            if debug:
                print("NOT WATER COVERED")
            # not fully covered, and within threshold of the water-covered part
            # simply transition from 0 to water-covered fraction
            snow_frac = snow_depth / self.snow_transition
            assert(snow_frac <=1)

        elif (snow_depth >= fully_covered):
            if debug:
                print("FULLY COVERED")
            # snow and water overtop the topography, everything is flat, simply transition from 0 to 1
            if (snow_depth >= self.snow_transition):
                snow_frac = 1.
            else:
                snow_frac = snow_depth / self.snow_transition
                assert(snow_frac <=1)

            
        else:
            if debug:
                print("WATER COVERED")
            # not fully covered, but more than just the water part
            # linearly transition between the above two cases
            y0 = water_frac
            x0 = self.snow_transition * water_frac
            y1 = 1.
            x1 = fully_covered
            assert(x1 - x0 > 0)
            snow_frac = y0 + (y1 - y0)/(x1-x0) * (snow_depth-x0)
            assert(snow_frac <=1)

        if debug:
            print('snow frac =', snow_frac)
        # subtract off snow from water, assume snow covers water-covered areas first (which are low-lying)
        if (snow_frac > water_frac):
            water_frac = 0.
        else:
            water_frac -= snow_frac

        # remainder goes to land
        land_frac = max(1. - water_frac - snow_frac,0.)

        try:
            assert(-1.e-10 <= land_frac <= 1+1.e-10)
            assert(-1.e-10 <= water_frac <= 1+1.e-10)
            assert(-1.e-10 <= snow_frac <= 1+1.e-10)
            assert(abs(1 - (land_frac + snow_frac + water_frac)) < 1.e-10)
        except AssertionError as err:
            print("broken at:")
            print(snow_depth, pd)
            print(land_frac, water_frac, snow_frac)
            raise err
        return land_frac, water_frac, snow_frac
    


def get_af(del_max, del_ex):
    return AreaFraction(del_max, del_ex)

def run_test(del_max, del_ex, fig, axs, bottom=False, alg='smoothed'):
    af = get_af(del_max, del_ex)

    rounding = 3
    
    pd = np.arange(-.01, 0.1301, 0.001)
    pres = pd * 1000. * 9.81 + 101325.0
    snow_depth = np.arange(0, 0.1301, 0.001)

    #    if alg is 'analytic':
    area_fracs = np.array([[af.area_fracs_analytic2(p, sd) for sd in snow_depth] for p in pres])
    #else:
    #    area_fracs = np.array([[af.area_fracs(p, sd) for sd in snow_depth] for p in pres])

    def decorate(ax, title, left=False, bottom=False):
        ax.set_title(title)
        if left:
            ax.set_ylabel('snow depth [m]')
        if bottom:
            ax.set_xlabel('ponded depth [m]')

        zero = np.where(np.abs(pd) < 1.e-8)[0]
        zero_arr = np.array([0.,])
        others = np.linspace(zero, len(pd)-1, 4)
        xticks = np.round(np.concatenate([zero_arr, others]))
        xticklabels = [np.round(pd[int(t)],2) for t in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)

        yticks = np.round(np.linspace(0, len(snow_depth)-1, 4))
        yticklabels = [np.round(snow_depth[int(t)],2) for t in yticks]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        
    
    axs[0].imshow(area_fracs[:,:,0].transpose(), origin='lower left', vmin=0, vmax=1, cmap='jet')
    decorate(axs[0], 'land frac', True, bottom)

    axs[1].imshow(area_fracs[:,:,1].transpose(), origin='lower left', vmin=0, vmax=1, cmap='jet')
    decorate(axs[1], 'water frac', False, bottom)

    mp = axs[2].imshow(area_fracs[:,:,2].transpose(), origin='lower left', vmin=0, vmax=1, cmap='jet')
    decorate(axs[2], 'snow frac', False, bottom)

    text_x = 0.05
    bounds = axs[0].get_position().bounds
    text_y = bounds[1] + 0.5*bounds[3]
    fig.text(text_x, text_y, 'del_max = %g, del_ex = %g'%(del_max, del_ex), rotation='vertical', ha='center', va='center')
    
    return mp
    
def test0():
    fig, axs = plt.subplots(2,3,squeeze=False)
    cax = fig.add_axes([0.92,0.1,0.03,0.8])

    mp = run_test(0, 0, fig, axs[0])
    mp = run_test(0.1, 0.05, fig, axs[1], True)
    plt.colorbar(mp, cax=cax)


def test1():
    fig, axs = plt.subplots(2,3,squeeze=False)
    cax = fig.add_axes([0.92,0.1,0.03,0.8])

    mp = run_test(0, 0, fig, axs[0], alg='analytic')
    mp = run_test(0.1, 0.05, fig, axs[1], True, alg='analytic')
    plt.colorbar(mp, cax=cax)
    plt.show()
    
