"""A single source for all units in known ATS variables (WIP!)"""

native_units = { 'pressure' : 'Pa',
                 'temperature' : 'K',
                 'ponded_depth' : 'm',
                 'saturation_liquid' : '-',
                 'saturation_gas' : '-',
                 'saturation_ice' : '-',
                 'snow-depth' : 'm',
}


def get_units(varname):
    if varname in native_units:
        return native_units[varname]

    if len(varname.split('-')) == 2:
        return get_units(varname[1])


