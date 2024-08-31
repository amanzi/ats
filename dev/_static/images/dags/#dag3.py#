from plot_dependency_graph2 import writePDF, texify, getNames # in $ATS_SRC_DIR/tools/utils
import networkx as nx


G = nx.DiGraph()

WC = getNames('water_content@NEXT')
WC2 = getNames('water_content@CURRENT')
sat = getNames('saturation_liquid@NEXT')
sat2 = getNames('saturation_liquid@CURRENT')
dens = getNames('density@NEXT')
dens2 = getNames('density@CURRENT')
poro = getNames('porosity@NEXT')
poro2 = getNames('porosity@CURRENT')
p = getNames('pressure@NEXT')
p2 = getNames('pressure@CURRENT')
T = getNames('temperature@NEXT')
T2 = getNames('temperature@CURRENT')
dWC = getNames('dwater_content_dtime@NEXT')


G.add_node(WC2[1], vtype='secondary', label=texify(WC2[0]))
G.add_node(sat2[1], vtype='secondary', label=texify(sat2[0]))
G.add_node(dens2[1], vtype='secondary', label=texify(dens2[0]))
G.add_node(poro2[1], vtype='independent', label=texify(poro2[0]))
G.add_node(p2[1], vtype='primary', label=texify(p2[0]))
G.add_node(T2[1], vtype='independent', label=texify(T2[0]))
G.add_edge(WC2[1], sat2[1])
G.add_edge(WC2[1], dens2[1])
G.add_edge(WC2[1], poro2[1])
G.add_edge(sat2[1], p2[1])
G.add_edge(dens2[1], p2[1])
G.add_edge(dens2[1], T2[1])

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



G.add_node(dWC[1], vtype='secondary', label=texify(dWC[0]))
G.add_edge(dWC[1], WC[1])
G.add_edge(dWC[1], WC2[1])

writePDF(G, 'dag3', scale=2, layout_scaling=(6,4))
