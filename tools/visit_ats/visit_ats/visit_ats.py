import sys, os
import visit_rcParams as vrc
import visit_colormaps as vc
import visit_time as vt
import visit as v

rcParams = dict( surface_linewidth=3,
                  snow_linewidth=4
                  )

class Operator:
    def __init__(self, oname, otype, oatts):
        self.oname = oname
        self.otype = otype
        self.oatts = oatts

class Plot:
    def __init__(self, pname, ptype, patts, varname=None):
        """Plot class"""
        self.pname = pname
        self.varname = varname
        self.ptype = ptype
        self.patts = patts
        self.operators = []
        self.annot = v.GetAnnotationObject(pname)
        self.annot.fontFamily = vrc.getDefaultFont()
        self.annot.fontHeight = vrc.rcParams['legend.fontheight']
        self.annot.managePosition = 0
        self.annot.xScale = vrc.rcParams['legend.scale'][0]
        self.annot.yScale = vrc.rcParams['legend.scale'][1]
        self.annot.position = vrc.rcParams['legend.position']
        self.annot.drawMinMax = vrc.rcParams['legend.minmax']

        if varname is not None:
            self.annot.drawTitle = 0
            self.title = v.CreateAnnotationObject("Text2D")
            self.title.text = varname
            self.title.fontFamily = vrc.getDefaultFont()
            self.title.height = vrc.rcParams['legend.title.fontheight']
            self.title.position = vrc.rcParams['legend.title.position']

    def getLimits(self):
        assert self.ptype == 'Pseudocolor'
        min = None
        max = None
        if self.patts.minFlag:
            min = self.patts.min
        if self.patts.maxFlag:
            min = self.patts.max
        return min,max

class VisItWindow:
    """Class for a window"""
    
    class Slice:
        """Helper class for slicing into 2D"""
        def __init__(self, point=None, normal=None):
            if point is None:
                point = (0,0,0)
            if normal is None:
                normal = (0,-1,0)

            assert type(point) is tuple
            assert len(point) == 3
            assert type(normal) is tuple
            assert len(normal) == 3
            self.point = point
            self.normal = normal

        def toAttributes(self):
            s = v.SliceAttributes()
            s.originType = s.Point
            s.originPoint = self.point
            s.normal = self.normal
            s.axisType = s.Arbitrary
            return s

    def __init__(self, index):
        self.i = index
        self.annot = vrc.getAnnotationAttributes()
        
        self.setDimension(3)
        self.plots = []
        self.nonplots = [] # other objects like meshes
        self._slice = None
        self.exaggeration = None
        
    def setDimension(self, dim):
        """Sets the dimension, which is used in controlling the view"""
        self.dim = dim
        if dim == 2:
            self.view = v.GetView2D()
        elif dim == 3:
            self.view = v.GetView3D()
        else:
            raise RuntimeError("Invalid dimension %s"%str(dim))

    def slice(self, point=None, normal=None):
        """Registers a slice -- this is not immediately added"""
        self._slice = self.Slice(point, normal)
        self.setDimension(2)

    def exaggerateVertical(self, factor):
        """Registers an exxageration -- this is not immediately added"""
        self.exaggeration = factor

    def _exaggerateVertical(self):
        if self.dim == 3:
            self.view.axis3DScaleFlag = 1
            self.view.axis3DScales = (self.view.axis3DScales[0],
                                      self.view.axis3DScales[1],
                                      self.exaggeration)
        else:
            for i,plot in enumerate(self.plots):
                done = False
                for op in plot.operators:
                    if "exaggerate_vertical" == op.oname:
                        done = True
                if not done:
                    print("transforming plot %d..."%i)
                    tr = v.TransformAttributes()
                    tr.doScale = 1
                    tr.scaleY = self.exaggeration
                    v.SetActivePlots(i)
                    v.AddOperator("Transform")
                    v.SetOperatorOptions(tr)
                    plot.operators.append(Operator("exaggerate_vertical", "Transform", tr))

    def createMesh(self, color='w', opacity=0.15, silo=False):
        _colors = dict(w=(255,255,255,255),
                       k=(0,0,0,255),
                       gray=(175,175,175),
                       )

        if silo:
            v.AddPlot('Mesh', "mesh")
        else:
            v.AddPlot('Mesh', "Mesh")
        ma = v.MeshAttributes()
        ma.legendFlag = 0
        ma.meshColor = _colors[color]
        ma.meshColorSource = ma.MeshCustom
        if (opacity < 1.):
            ma.opaqueMode = ma.On
            ma.opacity = opacity

        v.SetPlotOptions(ma)

        pname = v.GetPlotList().GetPlots(v.GetNumPlots()-1).plotName
        if silo:
            plot = Plot(pname, 'mesh', ma)
        else:
            plot = Plot(pname, 'Mesh', ma)
        self.nonplots.append(plot)
        return plot

        
    def createPseudocolor(self, varname, display_name=None, cmap=None, 
                          limits=None, linewidth=None, legend=True, alpha=False):
        """Generic creation of pseudocolor"""
        if display_name is None:
            display_name = vrc.renameScalar(varname)

        if "temperature" in display_name:
            display_name = display_name.replace("[K]", "[C]")
            print("defining alias: %s = %s"%(display_name, varname))
            v.DefineScalarExpression(display_name, "<%s> - 273.15"%varname)
        elif display_name != varname:
            print("defining alias: %s = %s"%(display_name, varname))
            v.DefineScalarExpression(display_name, '<'+varname+'>')

        v.AddPlot('Pseudocolor', display_name)
        pa = v.PseudocolorAttributes()

        # limits
        if limits is None:
            limits = vrc.getLimits(varname)

        if limits is not None:
            min = limits[0]
            max = limits[1]
            if min is not None:
                pa.minFlag = 1
                pa.min = min
            if max is not None:
                pa.maxFlag = 1
                pa.max = max

        # opacity
        if alpha:
            pa.opacity = 0
            pa.opacityType = pa.ColorTable

        # colormap
        if cmap is not None:
            reverse = cmap.endswith("_r")
            if reverse:
                cmap = cmap.strip("_r")
                pa.invertColorTable = 1
            pa.colorTableName = cmap

        # linewidth for 2D
        if linewidth is None:
            linewidth = vrc.rcParams['pseudocolor.linewidth']
        pa.lineWidth = linewidth

        # turn off legend for 2D surf
        if not legend:
            pa.legendFlag = 0

        v.SetActivePlots(len(self.plots)+1)
        v.SetPlotOptions(pa)
        pname = v.GetPlotList().GetPlots(v.GetNumPlots()-1).plotName
        if legend:
            plot = Plot(pname, 'Pseudocolor', pa, display_name)
        else:
            plot = Plot(pname, 'Pseudocolor', pa)
        self.plots.append(plot)
        return plot

    def createContour(self, varname, value, color=None, linewidth=None):
        """Generic creation of a single contour without a legend"""
        
        v.AddPlot('Contour', varname)
        ca = v.ContourAttributes()

        ca.contourMethod = ca.Value
        ca.contourValue = (value,)
        ca.colorType = ca.ColorBySingleColor

        if color is None:
            color = vrc.rcParams['contour.color']
        if type(color) is str:
            color = vc.common_colors[color]
        ca.singleColor = color

        if linewidth is None:
            linewidth = vrc.rcParams['contour.linewidth']
        ca.lineWidth = linewidth

        # turn off legend for 2D surf
        ca.legendFlag = 0

        v.SetPlotOptions(ca)
        pname = v.GetPlotList().GetPlots(v.GetNumPlots()-1).plotName
        plot = Plot(pname, 'Contour', ca)
        self.plots.append(plot)
        return plot

    def draw(self):
        print("drawing window %d of dimension %d"%(self.i,self.dim))
        v.SetActiveWindow(self.i)
        v.SetAnnotationAttributes(self.annot)

        if self.dim == 2:
            # add the slice
            assert self._slice is not None
            for i,plot in enumerate(self.plots):
                sliced = False
                for op in plot.operators:
                    if "slice" == op.oname:
                        sliced = True

                if not sliced:
                    print("slicing plot %d..."%i)
                    v.SetActivePlots(i)
                    v.AddOperator("Slice")
                    sa = self._slice.toAttributes()
                    v.SetOperatorOptions(sa)
                    plot.operators.append(Operator("slice", "Slice", sa))
                

        if self.exaggeration is not None:
            print("exaggerating...")
            self._exaggerateVertical()

        # set the plot options
        for i, plot in enumerate(self.plots):
            print("setting plot options for plot %i..."%i)
            v.SetActivePlots(i)
            v.SetPlotOptions(plot.patts)
        
        # set the view
        print("setting the view...")
        if self.dim == 2:
            v.SetView2D(self.view)
        else:
            v.SetView3D(self.view)

        print("drawing...")
        v.DrawPlots()
            
   
 
class Vis:
    """Container class for windows, also manages sources and correlations"""
    def __init__(self, directory, hostname="localhost", n_windows=1):
        self.directory = directory
        self.hostname = hostname
        self.windows = []
        self._active_window = 0

        for i in range(n_windows):
            self.addWindow()

        self.setActiveWindow(1)

        # if self.hostname != "localhost":
        #     args = []
        #     v.OpenMDServer(self.hostname,args)
        #     v.OpenComputeEngine(self.hostname,args)

    def loadSources(self, prefix="visdump_data", 
                    surface_prefix="visdump_surface_data", filetype="xdmf"):
        """Loads source files for subsurface and potentially surface."""

        if prefix is None:
            self.subsurface_src = None
        else:
            if filetype == "xdmf":
                self.subsurface_src = ":".join((self.hostname, os.path.join(self.directory, "%s.VisIt.xmf"%prefix)))
            elif filetype == "silo":
                self.subsurface_src = ":".join((self.hostname, os.path.join(self.directory, "%s.silo"%prefix)))

        if surface_prefix is None:
            self.surface_src = None
        else:
            if filetype == "xdmf":
                self.surface_src = ":".join((self.hostname, os.path.join(self.directory, "%s.VisIt.xmf"%surface_prefix)))
            elif filetype == "silo":
                self.surface_src = ":".join((self.hostname, os.path.join(self.directory, "%s.silo"%surface_prefix)))

        # open the subsurface database
        if self.subsurface_src is not None:
            v.OpenDatabase(self.subsurface_src)

        if surface_prefix is not None:
            # open the surface database
            v.OpenDatabase(self.surface_src)

            if prefix is not None:
                # create the database correlation
                v.CreateDatabaseCorrelation("my_correlation", (self.subsurface_src, self.surface_src), 0)
                v.SetActiveTimeSlider("my_correlation")

            # create vector expressions for ponded depth and snow depth
            v.DefineVectorExpression("ponded_depth_displace", "{0,0,ponded_depth.cell.0}")
            v.DefineVectorExpression("snow_displace", "{0,0,snow_depth.cell.0+ponded_depth.cell.0}")

    def loadSourcesList(self, srclist):
        """A generic set of sources."""
        self.src = [":".join((self.hostname, os.path.join(self.directory, s))) for s in srclist]
        for s in self.src:
            v.OpenDatabase(s)
            
    def unloadSources(self):
        v.DeleteAllPlots()
        if self.subsurface_src is not None:
            v.CloseDatabase(self.subsurface_src)
        if self.surface_src is not None:
            v.CloseDatabase(self.surface_src)
            
    def addWindow(self):
        """Adds a window to VisIt and makes it active"""
        if len(self.windows) != 0:
            v.AddWindow()

        win = VisItWindow(len(self.windows)+1)
        self.windows.append(win)
        self.setActiveWindow(len(self.windows))
        v.ToggleLockTime()

    def getActiveWindow(self):
        return self.windows[self._active_window-1]
        
    def setActiveWindow(self, i):
        assert 0 < i <= len(self.windows)
        self._active_window = i
        v.SetActiveWindow(i)

    def activateSurface(self):
        v.ActivateDatabase(self.surface_src)

    def activateSubsurface(self):
        v.ActivateDatabase(self.subsurface_src)

    def activateSource(self, src):
        v.ActivateDatabase(src)


class ATSVis(Vis):
    def __init__(self, *args, **kwargs):
        Vis.__init__(self, *args, **kwargs)
        self.time_annot = None
    
    def createSubsurfacePseudocolor(self, varname, limits=None, cmap=None, window=None):
        """Simplified interface to create standard pseudocolors"""
        self.activateSubsurface()

        if window is not None:
            self.setActiveWindow(window)
        win = self.getActiveWindow()
            
        return win.createPseudocolor(varname, limits=limits, cmap=cmap,
                                     legend=True)
    
    def createSurfacePseudocolor(self, varname, limits=None, cmap=None, window=None,
                                 displace=True, alpha=False, legend=False):
        """Simplified interface to create standard pseudocolors on the surface."""
        self.activateSurface()

        if window is not None:
            self.setActiveWindow(window)
        win = self.getActiveWindow()

        pcolor = win.createPseudocolor(varname, limits=limits, cmap=cmap, 
                                        legend=legend, alpha=alpha)
        if displace:
            # deform by surface vector
            v.AddOperator("Displace")
            da = v.DisplaceAttributes()
            da.variable = "ponded_depth_displace"
            v.SetOperatorOptions(da)
            pcolor.operators.append(Operator("displace", "Displace", da))
        return pcolor

    def createSnowPseudocolor(self, varname, limits=None, cmap=None, window=None, legend=False):
        """Simplified interface to create standard pseudocolors on the snow surface."""
        self.activateSurface()

        if window is not None:
            self.setActiveWindow(window)
        win = self.getActiveWindow()

        if cmap is None:
            cmap = "hot"
        
        pcolor = win.createPseudocolor(varname, limits=limits, cmap=cmap, 
                                       legend=legend)
                                    
        # deform by surface vector
        v.AddOperator("Displace")
        da = v.DisplaceAttributes()
        da.variable = "snow_displace"
        v.SetOperatorOptions(da)
        pcolor.operators.append(Operator("displace", "Displace", da))
        return pcolor

    def createSubsurfaceContour(self, varname, value, window=None, color=None, 
                                linewidth=None):
        self.activateSubsurface()
        
        if window is not None:
            self.setActiveWindow(window)
        win = self.getActiveWindow()
            
        return win.createContour(varname, value, color=color, linewidth=linewidth)
    
    def plotPressure(self, window=None, limits=None):
        """Adds a plot of subsurface pressure"""
        return self.createSubsurfacePseudocolor("pressure.cell.0", limits=limits, window=window)

    def plotSurfacePressure(self, window=None, limits=None):
        """Adds a plot of surface pressure"""
        return self.createSurfacePseudocolor("surface_pressure.cell.0", limits=limits, window=window)


    def plotLiquidSaturation(self, window=None, limits=None, cmap="saturation_liquid_r"):
        """Adds a plot of subsurface pressure"""

        return self.createSubsurfacePseudocolor("saturation_liquid.cell.0", limits=limits,
                                                cmap=cmap, window=window)

    def plotGasSaturation(self, window=None, limits=None):
        """Adds a plot of subsurface pressure"""
        return self.createSubsurfacePseudocolor("saturation_gas.cell.0", limits=limits,
                                                window=window)

    def plotIceSaturation(self, window=None, limits=None, cmap="saturation_ice_r"):
        """Adds a plot of subsurface pressure"""
        return self.createSubsurfacePseudocolor("saturation_ice.cell.0", limits=limits,
                                                cmap=cmap, window=window)

    def plotTemperature(self, window=None, limits=None):
        """Adds a plot of subsurface temperature"""
        # create the colormap
        cmap = None
        if limits is None:
            limits = vrc.rcParams['var.limits']["temperature"]

        if limits is not None:
            if limits[0] != None and limits[1] != None:
                vc.createTemperaturesColorMap(limits[0], limits[1], 0., "temperature")
                cmap = "temperature"

        return self.createSubsurfacePseudocolor("temperature.cell.0", limits=limits,
                                                cmap=cmap, window=window)

    def plotSurfaceTemperature(self, window=None, limits=None):
        """Adds a plot of surface temperature"""
        # create the colormap
        cmap = None
        if limits is None:
            limits = vrc.rcParams['var.limits']["temperature"]
        
        if limits != None:
            if limits[0] != None and limits[1] != None:
                vc.createTemperaturesColorMap(limits[0], limits[1], 0., "surface_temperature")
                cmap = "surface_temperature"

        return self.createSurfacePseudocolor("surface_temperature.cell.0", limits=limits,
                                             cmap=cmap, window=window)

    def plotSnowTemperature(self, window=None, limits=None):
        """Adds a plot of snow temperature"""
        # create the colormap
        cmap = None
        if limits is None:
            limits = vrc.rcParams['var.limits']["temperature"]

        if limits != None:
            if limits[0] != None and limits[1] != None:
                vc.createTemperaturesColorMap(limits[0], limits[1], 0., "snow_temperature")
                cmap = "snow_temperature"

        return self.createSnowPseudocolor("snow_temperature.cell.0", limits=limits,
                                          cmap=cmap, window=window)

    def plotPondedDepth(self, **kwargs):
        """Adds a plot of surface ponded depth"""
        if 'domain_name' in kwargs.keys():
            varname = kwargs['domain_name']+"-ponded_depth.cell.0"
            kwargs.pop("domain_name")
        else:
            varname = "surface-ponded_depth.cell.0"
        return self.createSurfacePseudocolor(varname, **kwargs)

    def plotSnowDepth(self, **kwargs):
        """Adds a plot of snow depth"""
        return self.createSnowPseudocolor("snow_depth.cell.0", **kwargs)

    def _getIndexByTime(self, time):
        pass

    def plotALD(self, yr_start, yr_end):
        """Adds a plot of contours of ALD from year start to year end, inclusive"""
        # find the time slice that starts the full period
        # for yr in range(yr_start, yr_end+1):
        #     # find the time slice that starts the year
        pass        
    
    def writeTime(self, round=None):
        self.setActiveWindow(vrc.rcParams['time.window'])

        if self.time_annot is None:
            ta = v.CreateAnnotationObject("Text2D")
            ta.position = vrc.rcParams['time.location']
            ta.height = vrc.rcParams['time.fontheight']
            ta.fontFamily = vrc.getDefaultFont()
            self.time_annot = ta

        if round is None:
            round = vrc.rcParams['time.round']
        self.time_annot.text = vt.visitTime(round)

    def plotSurfaceMesh(self, color='w', opacity=.15):
        """Simplified interface to create standard pseudocolors on the surface."""
        self.activateSurface()
        win = self.getActiveWindow()
        mesh = win.createMesh(color,opacity)
        return mesh

    def plotSubsurfaceMesh(self, color='w', opacity=.15):
        """Simplified interface to create standard pseudocolors on the surface."""
        self.activateSubsurface()
        win = self.getActiveWindow()
        mesh = win.createMesh(color,opacity)
        return mesh
    

    def draw(self):
        """Draw the plots"""
        for win in self.windows:
            win.draw()

        # leave with 1 as the active window, which seems to be
        # required to get saving to work
        self.setActiveWindow(1)

    def update(self):
        """Any changes not made by SetTimeSliderState()"""
        if self.time_annot is not None:
            self.writeTime()

                     

    
        
        
    
                    
    
            
