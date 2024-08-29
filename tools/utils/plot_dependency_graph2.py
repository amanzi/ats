import networkx as nx
import subprocess

_colors = {'independent':'C4973B',
           'secondary':'368CC9',#'032CFC',
           'primary':'EB3B14',
           'aliased':'5885af',
           'time_advanced':'832DD2'
           }

_names =  { 'pressure': ("p","p"),
            'temperature': ("T","T"),
            'saturation_liquid': ("s_l", "sl"),
            'enthalpy': ("H", "H"),
            'porosity': (r"\phi", "0"),
            'permeability': (r"\kappa","K"),
            'density': (r"\rho", "r"),
            'mol_frac_gas': (r"\omega", "2"),
            'viscosity_liquid': (r"\mu", "3"),
            'unfrozen_fraction': (r"\eta", "4"),
            'water_content': (r"\Theta", "5"),
            'energy': ("E", "E"),
            'internal_energy': ("u", "u"),
            'air_temperature': (r"T^{air}", "Tair"),
            'base_porosity': (r"\phi_0", "00"),
            'capillary_pressure_gas_liq': (r"p_c^{gl}", "pcgl"),
            'capillary_pressure_liq_ice': (r"p_c^{li}", "pcli"),
            'cell_volume': (r"d\Omega","6"),
            'elevation': ("z", "z"),
            'slope_magnitude': ("S","S"),
            'incoming_shortwave_radiation': (r"Q^{sw-in}_E", "QswinE"),
            'manning_coefficient': (r"m_n", "mn"),
            'overland_conductivity': (r"\kappa_s", "Ks"),
            'ponded_depth': ("h","h"),
            'surface-conducted_energy_source': (r"Q^{cond}_E", "Qcond"),
            'surface-mass_source': ("Q","Q"),
            'surface-total_energy_source': (r"Q^{tot}_E", "QEtot"),
            'thermal_conductivity': (r"\kappa^T", "KT"),
            'vapor_pressure': ("vp", "vp"),
            'wind_speed': (r"U_{wind}", "Uwind"),
            'snow-depth': (r"h_{snow}", "hsnow"),
            'relative_humidity': (r"RH", "RH"),
            'surface-precipitation_rain': (r"P^r", "Pr"),
            'snow-precipitation': (r"P^s", "Ps"),
            'pres_elev': (r"h+z","hpz"),
            'relative_permeability': (r"k_r", "kr"),
            'dwater_content_dtime': (r"\frac{d\Theta}{dt}", "dWC"),
           }

document_preamble = \
r"""
\documentclass{{report}}
\usepackage{{tikz}}
\usepackage[active,tightpage]{{preview}}
\PreviewEnvironment{{tikzpicture}}
\setlength\PreviewBorder{{5pt}}
\usepackage{{subcaption}}
\usepackage{{xcolor}}
"""

document_beginend = \
r"""
\begin{{document}}
{content}
\end{{document}}
"""

def splitKeyTag(keytag):
    if '@' not in keytag:
        return keytag, None
    else:
        split = keytag.split('@')
        assert len(split) == 2
        return split[0], split[1]

def getKeyTag(key, tag, tex=False):
    if tag is None:
        return key
    else:
        if tex:
            return key+f"\\;\;_{{@{tag}}}"
        else:
            return key+'@'+tag

def getNames(keytag):
    """Given an ATS-style key@tag, returns two strings, one for use as a node name, the other as the label."""
    key, tag = splitKeyTag(keytag)
    if key in _names:
        k1,k2 = _names[key]
    else:
        k1,k2 = key,key
        k1.replace("_", "\\_")

    return getKeyTag(k1, tag, True), getKeyTag(k2, tag)

def texify(name):
    if name is None:
        return None
    return rf"${name}$"

def writeLatex(G, output_name='output',
               layout='igraph', layout_scaling=(3,2),
               scale=1,
               fontsize='normalsize'):
    """Plots a networkx graph.

    Expects a networkx graph that includes vertex attributes:
    - label: the variable name
    - vtype: the variable type (one of _colors.keys())
    """
    if not output_name.endswith('.tex'):
        output_name = output_name + '.tex'

    # use igraph to get a layout and then tikz position
    if layout == 'igraph':
        import igraph
        IG = igraph.Graph.from_networkx(G)
        layout = IG.layout('rt')
        layout.mirror(1)
        if layout_scaling is not None:
            layout.scale(*layout_scaling)
        
        position = dict(zip(G.nodes, layout))
    elif layout == 'graphviz' or layout == 'graphviz: dot':
        position = nx.nx_agraph.graphviz_layout(G, prog='dot')

    # get node drawing options
    node_opts = dict()
    for n in G.nodes:
        vtype = G.nodes[n]['vtype']
        node_opts[n] = f"rectangle, shading=axis, left color={vtype}!60!white, right color = {vtype}!30!white,"\
                       f"shading angle=-45, anchor=north,"\
                       f"style={{draw={vtype}, align=center}},"\
                       f"minimum height=2em, minimum width=4em, font=\\{fontsize},"\
                       f"scale={scale}"

    # get the document context
    document_colors = ["\\definecolor{{"+key+"}}{{HTML}}{{"+value+"}}" for (key,value) in _colors.items()]
    document_wrapper = "\n".join([document_preamble,]+document_colors+[document_beginend,])

    nx.write_latex(G, output_name, as_document=True,
               pos=position,
               node_options=node_opts,
               node_label='label',
               document_wrapper=document_wrapper
               )
    return output_name


def writePDF(G, *args, **kwargs):
    output_name = writeLatex(G, *args, **kwargs)
    subprocess.run(['pdflatex', output_name])
    subprocess.run(['convert', output_name[0:-3]+'pdf', output_name[0:-3]+'png'])



    
    
    
