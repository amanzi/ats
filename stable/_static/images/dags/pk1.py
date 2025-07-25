from plot_pk_tree import writePDF # in $ATS_SRC_DIR/tools/utils
import networkx as nx

G = nx.DiGraph()
G.add_node("coupled water", vtype="mpc", label="Coupled Water \\\\ \\hrulefill \\\\ MPC: Custom")
G.add_node("richards", vtype="physical", label="Richards Flow")
G.add_node("overland", vtype="physical", label="Overland Flow")
G.add_edge("coupled water", "richards")
G.add_edge("coupled water", "overland")

writePDF(G, 'pk1', scale=2, layout_scaling=(6,4))

