from plot_pk_tree import writePDF # in $ATS_SRC_DIR/tools/utils
import networkx as nx

G = nx.DiGraph()
G.add_node("permafrost", vtype="mpc", label="Permafrost \\\\ \\hrulefill \\\\ MPC: Custom")
G.add_node("richards", vtype="physical", label="Richards Flow")
G.add_node("sub_energy", vtype="physical", label="Subsurface Energy\\\\Transport")
G.add_node("overland", vtype="physical", label="Overland Flow")
G.add_node("surf_energy", vtype="physical", label="Surface Energy\\\\Transport")
G.add_edge("permafrost", "richards")
G.add_edge("permafrost", "overland")
G.add_edge("permafrost", "sub_energy")
G.add_edge("permafrost", "surf_energy")

G.add_node("coupler", vtype="mpc", label="Top Level Coupler \\\\ \\hrulefill \\\\ MPC: Weak")
G.add_node("snow", vtype="physical", label="Snow Water\\\\Conservation")
G.add_edge("coupler", "permafrost")
G.add_edge("coupler", "snow")


writePDF(G, 'pk2', scale=2, layout_scaling=(6,4))

