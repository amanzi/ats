from plot_dependency_graph2 import writePDF, texify, getNames # in $ATS_SRC_DIR/tools/utils
import networkx as nx


G = nx.DiGraph()

WC = getNames('water_content@NEXT')
WC2 = getNames('water_content@CURRENT')
sat = getNames('saturation_liquid@NEXT')
dens = getNames('density@NEXT')
poro = getNames('porosity@NEXT')
p = getNames('pressure@NEXT')
T = getNames('temperature@NEXT')
dWC = getNames('dwater_content_dtime@NEXT')

G.add_node(WC2[1], vtype='time_advanced', label=texify(WC2[0]))
G.add_node(dWC[1], vtype='secondary', label=texify(dWC[0]))

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

G.add_edge(dWC[1], WC[1])
G.add_edge(dWC[1], WC2[1])
G.add_edge(WC2[1], WC[1])

writePDF(G, 'dag2', scale=2, layout_scaling=(6,4))
