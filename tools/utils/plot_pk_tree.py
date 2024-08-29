import networkx as nx
import subprocess

_colors = {'physical':'ff8811',
           'mpc':'9dd9d2',
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
                       f"scale={scale}, rounded corners"

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
