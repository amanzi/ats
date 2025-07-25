from plot_dependency_graph2 import writePDF, texify, getNames # in $ATS_SRC_DIR/tools/utils
import networkx as nx

G = nx.DiGraph()

WC = getNames('water_content')
sat = getNames('saturation_liquid')
dens = getNames('density')
poro = getNames('porosity')
p = getNames('pressure')
T = getNames('temperature')

G.add_node(WC[1], vtype='secondary', label=texify(WC[0]))
G.add_node(sat[1], vtype='secondary', label=texify(sat[0]))
G.add_node(dens[1], vtype='secondary', label=texify(dens[0]))
G.add_node(poro[1], vtype='independent', label=texify(poro[0]))
G.add_node(p[1], vtype='primary', label=texify(p[0]))
G.add_node(T[1], vtype='independent', label=texify(T[0]))
G.add_edge(WC[1], sat[1])
G.add_edge(WC[1], dens[1])
G.add_edge(WC[1], poro[1])
G.add_edge(sat[1], p[1])
G.add_edge(dens[1], p[1])
G.add_edge(dens[1], T[1])

writePDF(G, 'dag1', scale=2, layout_scaling=(6,4))
