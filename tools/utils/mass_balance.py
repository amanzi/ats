"""Mass Balance 

This is a notebook demonstrating how mass balance calculations are
best done for ATS simulations.

Note that all of these are approximate, and should not expect to
perfectly conserve mass.  The reason for this is that ATS will be
simulating at a higher time resolution than your vis, and therefore
the integrals implied here are not exact.  Similarly, to get an even
reasonable approximation, make sure that, for the time period you are
interested in, you use vis that is on the same or higher resolution
than your forcing data.

Author: Ethan Coon (coonet@ornl.gov)

## Surface Mass Balance

The surface mass balance is given by a balance of precipitation,
infiltration, runoff, and evaporation:

```d Theta / dt = P*n*A - I - R - E*n*A```

where:

* ```Theta``` is ```surface-water_content [mol]```,
* ```P``` is ```surface-precipitation_rain``` plus ```surface-precipitation_snow [m/s]```,
* ```I``` is infiltration, the negative of ```surface-surface_subsurface_flux [mol/s]```
* ```R``` runoff ```[mol/s]```,
* ```E``` evaporation (which is actually evaporation - condensation), the _negative_ of ```surface-evaporative_flux```, (the negative is a consistency bug in the code, fix me!) ```[m/s]```
* ```n``` is the ```surface-molar_density_liquid``` of water ```[mol/m^3]```,
* and ```A``` is the ```surface-cell_volume [m^2]```.  

Note that snow precipitation can be directly added to rain
precipitation, as snow is assumed to be in Snow-Water-Equivalent (SWE)
units.

Note that some quantities here are extensive flow rates (infiltration,
runoff, water content) while others are intensive fluxes (precip,
evaporation).  In an ideal world, these would be changed in the code
to all be consistent (fix me, contributions welcome!)


## Subsurface Mass Balance

The subsurface mass balance is given by:
    
```d Theta / dt = I - S - T*V```

where 

* ```Theta``` is ```water_content [mol]```,
* ```I``` is infiltration, the negative of ```surface-surface_subsurface_flux [mol/s]```
* ```S``` seepage ```[mol/s]```, 
* ```T``` is ```transpiration [mol/m^3/s]```
* and ```V``` is ```cell_volume [m^3]```


## Simulation instrumentation

There are two ways to measure the needed quanties discussed above: by post-processing the visualization files and by adding observations.  The advantage of post-processing the vis files is that you needn't have thought of this before running your run.  The advantage of using observations is that you can get away with higher resolution (temporally) because the integrations are performed in the code, so you don't need to save huge quantities of data.  Visualization files are not _entirely_ sufficient, however, as runoff/seepage are not post-processable.  

Therefore there are two strategies discussed here.  

* For a run with vis and no observations, runoff/seepage can be inferred by solving the above equations for ```R``` and ```S```, respectively, and post-processing everything else.
* For a run with a full set of observations, check that mass is conserved (modulo time sampling error) and understand the balance of terms.

### Visualization

Throughout this workbook we will assume that visualization has been added to the input file, and looks something like this:

```xml
  <ParameterList name="visualization">
    <ParameterList name="domain" type="ParameterList">
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1}" />
    </ParameterList>
    <ParameterList name="surface" type="ParameterList">
      <Parameter name="times start period stop" type="Array(double)" value="{0.0,86400.0,-1}" />
    </ParameterList>
  </ParameterList>
```

This generates daily vis, and is consistent with a daily prescribed precipitation forcing.


"""


import numpy as np
from matplotlib import pyplot as plt

import ats_xdmf # from $ATS_SRC_DIR/tools/utils/


names_dev = {"ponded depth":"surface-ponded_depth.cell.0",
         "surface pressure":"surface-pressure.cell.0",
         "pressure":"pressure.cell.0",
         "rain":"surface-precipitation_rain.cell.0",
         "snow":"surface-precipitation_snow.cell.0",       
         "water content":"water_content.cell.0",
         "surface water content":"surface-water_content.cell.0",
         "saturation":"saturation_liquid.cell.0",
         "saturation gas":"saturation_gas.cell.0",
         "saturation ice":"saturation_ice.cell.0",
         "cell volume":"cell_volume.cell.0",
         "surface cell volume":"surface-cell_volume.cell.0",
         "surface mass density":"surface-mass_density_liquid.cell.0",
         "exfiltration":"surface-surface_subsurface_flux.cell.0",
         "evaporation":"surface-evaporative_flux.cell.0",
         "transpiration":"transpiration.cell.0",
         "surface density": "surface-molar_density_liquid.cell.0",
         "density": "molar_density_liquid.cell.0",
        }


names_086 = {"ponded depth":"ponded_depth.cell.0",
         "surface pressure":"surface-pressure.cell.0",
         "pressure":"pressure.cell.0",
         "rain":"precipitation_rain.cell.0",
         "snow":"precipitation_snow.cell.0",       
         "water content":"water_content.cell.0",
         "surface water content":"surface-water_content.cell.0",
         "saturation":"saturation_liquid.cell.0",
         "saturation gas":"saturation_gas.cell.0",
         "saturation ice":"saturation_ice.cell.0",
         "cell volume":"cell_volume.cell.0",
         "surface cell volume":"surface-cell_volume.cell.0",
         "surface mass density":"surface-mass_density_liquid.cell.0",
         "exfiltration":"surface_subsurface_flux.cell.0",
         "evaporation":"evaporative_flux.cell.0",
         "surface density": "surface-molar_density_liquid.cell.0",
         "density": "molar_density_liquid.cell.0",
        }


class MassBalanceFromVis(object):
    """Simulation class encapsulating reading output for common plotting work."""
    _names = names_dev

    @classmethod
    def set_names_dev(cls):
        """Change this class to use dev variable names.

        Usage:
          MassBalanceFromVis.set_names_dev()
          sim = MassBalanceFromVis(...)
        """
        cls._names = names_dev

    @classmethod
    def set_names_086(cls):
        """Change this class to use dev variable names.

        Usage:
          MassBalanceFromVis.set_names_dev()
          sim = MassBalanceFromVis(...)
        """
        cls._names = names_086
    

    def __init__(self, dirname, names=None, gravity=9.80665, typical_density=55000.):
        """Create a simulation class.
        
        Usage:
          with MassBalanceFromVis(dirname) as sim:
              ...
        
        Arguments:
          dirname        | The directory where output is.
          names          | A dictionary mapping common names to variable names in the 
                         | visualization files.  Defaults are provided, and most 
                         | simulations need not override them. (default=None)
          gravity        | Magnitude of gravity (default=9.80665 is the standard ATS value used)
          typical_density| A typical molar density of liquid, used to print things in [m] instead of [mol]
        
        Note that this should typically be used in 'with ... as' context, as this opens file 
        resources which should be closed on exit.  Alternatively, close() can be called directly
        by the user.
        """
        self.dirname = dirname
        if names is not None:
            self._names.update(names)
        self.gravity = gravity
        self.p_atm = 101325.0
        
        # load the vis files
        self.vis = ats_xdmf.VisFile(dirname, time_unit='d')
        self.vis.loadMesh()
        
        self.vis_surf = ats_xdmf.VisFile(dirname, domain='surface', time_unit='d')
        self.vis_surf.loadMesh()

        self.length = min(len(self.vis.cycles), len(self.vis_surf.cycles))
        self.times = self.vis.times[0:self.length]

        # save cell volume instead of loading it repeatedly
        self.volumes = self.vis.get("cell_volume", self.vis.cycles[0])
        self.volume = self.volumes.sum()
        self.areas = self["surface cell volume",0]
        self.surface_area = self.areas.sum()
        
        self.typical_density = typical_density
        

    def __getitem__(self, name_and_index):
        """Reads data from a simulation, based on indices into time/time_s.
        
        Example:
          s = Sim(dirname)
          s["saturation", -1]  # returns saturation at the final timestep
          
        name_and_index is a tuple of:

        name           | Common name of the variable, must be in the names dict.
        i              | Index into timesteps, i.e. returns the value at s.times[i].
        """
        name, i = name_and_index
        # ponded depth is a special case as it is important, but doesn't always exist.
        # instead derive it if needed
        if name == "ponded depth":
            internal_name = self._names[name]
            try:
                return self.vis_surf.get(internal_name, self.vis_surf.cycles[i])
            except KeyError:
                pres = self["surface pressure", i]
                return np.maximum((pres - self.p_atm) / self.density / self.gravity, 0.)
        
        try:
            internal_name = self._names[name]
        except KeyError:
            raise KeyError('Unknown common name: "%s", call self.common_names() to see what'%name+
                           'are available, or add to the names dict.')
        else:
            if self.isSurface(name):
                return self.vis_surf.get(internal_name, self.vis_surf.cycles[i])
            else:
                return self.vis.get(internal_name, self.vis.cycles[i])
            
    def __enter__(self):
        """Supports 'with MassBalanceFromVis(...) as...' semantics."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Supports 'with MasBalanceFromVis(...) as ...' semantics."""
        self.close()
    
    def close(self):
        """Close all file resources."""
        self.vis.close()
        self.vis_surf.close()

    def print_common_names(self):
        """Prints a list of recognized variable names"""
        print("Recognized names are:")
        for n in self._names.keys():
            print(" ", n)
            
    def __len__(self):
        """Provides length of the datasets.  Syntatic sugar."""
        return self.length
    
    def isSurface(self, name):
        """Is this variable a surface domain quantity?"""
        if self._names[name].startswith('surface'):
            return True
        # 0.86!
        elif self._names[name].startswith('precipitation') or self._names[name].startswith('evaporative'):
            return True
        
        else:
            return False
    
    def integrate(self,name,i, scale_by_density=False):
        """Spatially integrates intensive data over the domain at timestep i.
        
        If scale_by_density is set, this additional multiplies by the molar density.
        """
        if not scale_by_density:
            if self.isSurface(name):
                return (self[name,i]*self.areas).sum()
            else:
                return (self[name,i]*self.volumes).sum()
        else:
            if self.isSurface(name):
                return (self[name,i]*self.areas*self['surface density',i]).sum()
            else:
                return (self[name,i]*self.volumes*self['density',i]).sum()
            
    def sum(self, name, i):
        """Spatially integrates extensive data over the domain at timestep i."""
        return self[name,i].sum()

    def surfaceWC(self,i):
        """surface WC in [mol] at timestep i"""
        return self.sum('surface water content',i)
    
    def evaporation(self, i):
        """evaporation in [mol/s]:
        
        Note that some simulations may not include surface energy balance!
        """
        try:
            return self.integrate('evaporation',i,True)
        except KeyError:
            return 0.0
        
    def precipitation(self, i):
        """precipitation in [mol/s]"""
        rain = self.integrate('rain',i,True)
        try:
            snow = self.integrate('snow',i,True)
        except KeyError:
            snow = 0.
        return rain + snow
    
    def infiltration(self, i):
        """infiltration in [mol/s]"""
        return -self.sum('exfiltration',i)
    
    def runoff(self, i, log=False):
        """runoff [mol/s] in the interval (i,i+1).  
        
        Calculated via a surface water mass balance:
        
        dWC/dt = precip - infiltration - runoff
        
        As a result this is an approximation, and may be wrong for a variety of reasons:
          - assumes visualization is at least as frequent as precip data
          - assumes infiltration is relatively smooth
        
        If any of these caveats affects you, you might find that introducing a 
        runoff observation is a better choice than using this function.
        """
        # precip in mol / s
        precip = (self.precipitation(i+1) + self.precipitation(i)) / 2.0
        
        # infiltration in mol / s
        infiltration = (self.infiltration(i+1) + self.infiltration(i)) / 2.0
        
        # evap in mol/s
        evap = (self.evaporation(i+1) + self.evaporation(i)) / 2.0
        
        # dwc/dt in mol / s
        dwc = (self.surfaceWC(i+1) - self.surfaceWC(i)) \
                 / ((self.times[i+1]-self.times[i]) * 86400.)

        # runoff in mol / s
        runoff = precip - infiltration - evap - dwc
        
        if log:
            print("Surface mass balance (surface area [m^2] =", self.surface_area, "):")
            print("  precip [m/day]   = ", precip / self.surface_area / self.typical_density * 86400)
            print("  infilt [m/day]   = ", infiltration / self.surface_area / self.typical_density * 86400)
            print("  evap   [m/day]   = ", evaporation / self.surface_area / self.typical_density * 86400)
            print("  dWC/dt [m/day]   = ", dwc / self.surface_area / self.typical_density * 86400)
            print("  runoff [m^3/s]   = ", runoff / self.typical_density)
        return runoff
        
    def WC(self,i=None):
        """subsurface WC in [mol] at timestep i"""
        return self.sum('water content', i)
                
    def saturatedWC(self,i):
        """saturated zone WC in [mol] at timestep i"""
        return np.where(self['saturation gas',i] < 0.001,
                        self['water content',i], 0.0).sum()
    
    def unsaturatedWC(self,i):
        """unsaturated zone WC in [mol]"""
        return np.where(self['saturation gas',i] >= 0.001,
                        self['water content',i], 0.0).sum()
    
    def transpiration(self, i):
        """transpiration in [mol/s] is in the negative direction (it is a sink)"""
        try:
            return -self.integrate('transpiration',i)
        except KeyError:
            return 0.0
    
    def seepage(self, i, log=False):        
        """seepage [mol/s] in the interval (i,i+1).  
        
        Calculated via a subsurface water mass balance:
        
        dWC/dt = infiltration - transpiration - seepage
        
        As a result this is an approximation, and may be wrong for a variety of reasons:
          - assumes visualization is at least as frequent as precip data
          - assumes infiltration is relatively smooth
        
        If any of these caveats affects you, you might find that introducing a 
        runoff observation is a better choice than using this function.
        """

        # infiltration in mol / s
        infiltration = (self.infiltration(i+1) + self.infiltration(i)) / 2.0
        
        # transpiration in mol/s
        trans = (self.transpiration(i+1) + self.transpiration(i)) / 2.0
        
        # dwc/dt in mol / s
        dwc = (self.WC(i+1) - self.WC(i)) \
                 / ((self.times[i+1]-self.times[i]) * 86400.)

        # seepage in mol / s
        seepage = infiltration - trans - dwc
        
        if log:
            print("Subsurface mass balance (surface area [m^2] =", self.surface_area, "):")
            print("  infiltration [m/day]  = ", infiltration / self.surface_area / self.typical_density * 86400)
            print("  transpiration [m/day] = ", trans / self.surface_area / self.typical_density * 86400)
            print("  dWC/dt [m/day]   = ", dwc / self.surface_area / self.typical_density * 86400)
            print("  seepage [m^3/s]   = ", seepage / self.typical_density)
        return seepage
        
    def vectorize(self, method):
        """Calculates at all timesteps and returns an array"""
        return np.array([getattr(self,method)(i) for i in range(len(self))])
    
    def half_times_vectorize(self, method):
        """Calculates at all half-timesteps and returns an array"""
        return np.array([getattr(self,method)(i) for i in range(len(self)-1)])
    
    def half_times(self):
        return (self.times[1:] + self.times[:-1])/2.0



# plotting functions
def get_axs(figsize=None):
    """Gets a list of axes for use in mass balance plots."""
    if figsize is None:
        figsize = (10,12)
    fig, axs = plt.subplots(3,2, figsize=figsize)
    return axs.flatten()

def plot(sim, axs, color='b', symbol=None, label=None, derived_runoff=True, derived_seepage=True):
    """Plots the mass balance for a given simulation."""
    if symbol is None:
        style = '-'
        dash = '--'
    else:
        style = '-'+symbol
        dash = '--'+symbol
    
    # units to convert from mol/s to m/day
    conv = 86400.0 / sim.surface_area / sim.typical_density
    precip = sim.vectorize('precipitation')*conv
    axs[0].plot(sim.times, precip, style, color=color, label=label)
    axs[1].plot(sim.times, (sim.vectorize('transpiration') + sim.vectorize('evaporation'))               
                *conv, dash, color=color)

    to_m = 1.0 / sim.surface_area / sim.typical_density
    axs[2].plot(sim.times, sim.vectorize('surfaceWC')*to_m, style, color=color)
    axs[3].plot(sim.times, sim.vectorize('WC')*to_m, style, color=color)

    if derived_runoff:
        half_times = sim.half_times()
        axs[4].plot(half_times, sim.half_times_vectorize('runoff')/sim.typical_density, style, color=color)
        axs[4].plot(half_times, (precip[1:] + precip[:-1])/2. * sim.surface_area / 86400.0, 'k--')
    else:
        axs[4].plot(sim.times, sim.vectorize('runoff')/sim.typical_density, style, color=color)
        axs[4].plot(sim.times, precip * sim.surface_area / 86400.0, 'k--')

    if derived_seepage:
        half_times = sim.half_times()
        axs[5].plot(half_times, sim.half_times_vectorize('seepage')/sim.typical_density, style, color=color)
    else:
        axs[5].plot(sim.times, sim.vectorize('seepage')/sim.typical_density, style, color=color)

    
def decorate(axs):
    """Decorates the axes with labels"""
    for ax in axs:
        ax.set_xlabel("time [days]")
    axs[0].set_ylabel("precip [m/dy]")
    axs[1].set_ylabel("E+T [m/dy]")
    axs[2].set_ylabel("surface WC [m]")
    axs[3].set_ylabel("subsurface WC [m]")
    axs[4].set_ylabel("runoff [m^3/s]")
    axs[5].set_ylabel("seepage [m^3/s]")
    
    #    axs[1].set_ylim(0.015,0.025)
    #    axs[4].set_ylim(0,2)

def legend(axs):
    axs[0].legend()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot mass balance of a set of runs.")
    parser.add_argument('directories', metavar='dirs', type=str, nargs='+',
                        help="list directories to plot")
    parser.add_argument('--dev', action="store_true", default=False,
                        help="use ats-dev variable names")
    parser.add_argument('--names', type=str, default=None,
                        help="string containing a pythonic dictionary of variable names")
    parser.add_argument('--symbol', '-s', default=None, type=str,
                        help="symbol to place on the line")
    args = parser.parse_args()

    # axes
    axs = get_axs()

    # color iterator
    import itertools
    colors = itertools.cycle( (c for c in "bgrcmykw") )

    # parse optional names
    if args.names is not None:
        import ast
        names = dict(ast.literal_eval(args.names))
    else:
        names = None
    
    if args.dev:
        MassBalanceFromVis.set_names_dev()

    # loop and plot
    for d in args.directories:
        sim = MassBalanceFromVis(d, names=names)
        color = colors.next()
        plot(sim, axs, color, symbol=args.symbol, label=d)

    decorate(axs)

    if len(args.directories) > 1:
        legend(axs)
        
    plt.show()
    
    
