import visit as v

common_colors = dict(w=(255,255,255,255),
                     k=(0,0,0,255),
                     r=(255,0,0,255),
                     g=(0,255,0,255),
                     b=(0,0,255,255)
                     )


_added_temps = []
def createTemperaturesColorMap(Tmin, Tmax, zero=0., name="temperature"):
    """Creates a ColorMap with warm colors above zero, and cold colors below."""
    # color control point format: (position, (r,g,b,a)) for (r,g,b,a) in [0,255]
    if name in _added_temps:
        print("Already added a temperature color map of name %s"%name)
        return

    if Tmin >= Tmax:
        raise RuntimeError("Tmin >= Tmax")
    
    warms = [(247,   0, 255, 255),
             (255,   0,   9, 255),
             (255, 119,   0, 255),
             (255, 247,   0, 255)]
    warms.reverse()
    cools = [(  9, 255,   0, 255),
    #             (  0, 255, 119, 255),
             (  0, 255, 247, 255),
             (  0, 137, 255, 255),
    #             (  0,   9, 255, 255),
             (119,   0, 255, 255)]
    cools.reverse()

    def _interp(mymin, mymax, colors):
        dp = (mymax - mymin) / (len(colors) - 1)
        pos = [mymin + i*dp for i in range(len(cools))]
        pos[-1] = mymax # roundoff issues
        return zip(pos, colors)
    
    # set up the positions
    if Tmax <= zero:
        cmap = _interp(0,1,cools)
    elif Tmin >= zero:
        cmap = _interp(0,1,warms)
    else:
        myzero = (zero - Tmin) / (Tmax - Tmin)
        cmap = _interp(0, myzero-0.00001, cools) \
              + _interp(myzero, 1, warms)

    clist = v.ColorControlPointList()
    for val in cmap:
        cp = v.ColorControlPoint()
        cp.position = val[0]
        cp.colors = val[1]
        clist.AddControlPoints(cp)
    print("Colormap Temp", name, [val[0] for val in cmap])

    v.AddColorTable(name, clist)
    _added_temps.append(name)


_sat_ice_added = False
def createSaturationIceColorMap():
    global _sat_ice_added
    if _sat_ice_added:
        return
    
    # cmap = [(0,      (125, 145, 198, 255)),
    #         (0.2584, (129, 203, 255, 255)),
    #         (0.7116, (255, 255, 163, 157)),
    #         (1.0,    (161, 128,  77,  74))]
    cmap = [(0,      (115, 145, 198, 255)),
            (0.2084, (129, 203, 255, 255)),
            (0.4584, (124, 186,  53, 255)),
            (0.7116, (255, 255, 163, 157)),
            (1.0,    (161, 128,  77,  74))]
    clist = v.ColorControlPointList()
    for val in cmap:
        cp = v.ColorControlPoint()
        cp.position = val[0]
        cp.colors = val[1]
        clist.AddControlPoints(cp)

    v.AddColorTable("saturation_ice", clist)
    _sat_ice_added = True


_sat_liquid_added = False
def createSaturationLiquidColorMap():
    global _sat_liquid_added
    if _sat_liquid_added:
        return
    
    cmap = [(0,      (  0,   0, 255, 255)),
            (0.2584, (  0, 255, 255, 255)),
            (0.7116, (255, 255,   0, 170)),
            (1.0,    ( 93,  64,   0,  98))]
    clist = v.ColorControlPointList()
    for val in cmap:
        cp = v.ColorControlPoint()
        cp.position = val[0]
        cp.colors = val[1]
        clist.AddControlPoints(cp)

    v.AddColorTable("saturation_liquid", clist)
    _sat_liquid_added = True


_ponded_water_added = False
def createPondedWaterColorMap():
    global _ponded_water_added    
    if _ponded_water_added:
        return
    
    cmap = [(0,      (  0,   0, 255, 255)),
            (0.2584, (  0, 255, 255, 255)),
            (0.7116, (255, 255,   0, 170)),
            (1.0,    ( 93,  64,   0,  98))]
    clist = v.ColorControlPointList()
    for val in cmap:
        cp = v.ColorControlPoint()
        cp.position = val[0]
        cp.colors = val[1]
        clist.AddControlPoints(cp)

    v.AddColorTable("ponded_water", clist)
    _ponded_water_added = True
    

createSaturationIceColorMap()
createSaturationLiquidColorMap()
createPondedWaterColorMap()
