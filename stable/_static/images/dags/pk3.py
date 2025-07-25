from plot_pk_tree import writePDF # in $ATS_SRC_DIR/tools/utils
import networkx as nx

G = nx.DiGraph()

G.add_node("tlc", vtype="mpc", label="Top Level Coupler \\\\ \\hrulefill \\\\ MPC: Custom")
G.add_node("surf_star", vtype="mpc", label="Surface Lateral \\\\ \\hrulefill \\\\ MPC: Custom")
G.add_node("surf_flow", vtype="physical", label="Surface Flow")
G.add_node("surf_energy", vtype="physical", label="Surface Energy \\\\ Transport")
G.add_edge("tlc", "surf_star")
G.add_edge("surf_star", "surf_flow")
G.add_edge("surf_star", "surf_energy")

G.add_node("column_coupler", vtype="mpc", label="Column Coupler \\\\ \\hrulefill \\\\ MPC: Subcycling")
G.add_edge("tlc", "column_coupler")


G.add_node("permafrost1", vtype="mpc", label="")
G.add_node("richards1", vtype="physical", label="")
G.add_node("sub_energy1", vtype="physical", label="")
G.add_node("overland1", vtype="physical", label="")
G.add_node("surf_energy1", vtype="physical", label="")
G.add_edge("permafrost1", "richards1")
G.add_edge("permafrost1", "overland1")
G.add_edge("permafrost1", "sub_energy1")
G.add_edge("permafrost1", "surf_energy1")
G.add_node("coupler1", vtype="mpc", label="Column 0")
G.add_node("snow1", vtype="physical", label="")
G.add_edge("coupler1", "permafrost1")
G.add_edge("coupler1", "snow1")
G.add_edge("column_coupler", "coupler1")

G.add_node("permafrost2", vtype="mpc", label="")
G.add_node("richards2", vtype="physical", label="")
G.add_node("sub_energy2", vtype="physical", label="")
G.add_node("overland2", vtype="physical", label="")
G.add_node("surf_energy2", vtype="physical", label="")
G.add_edge("permafrost2", "richards2")
G.add_edge("permafrost2", "overland2")
G.add_edge("permafrost2", "sub_energy2")
G.add_edge("permafrost2", "surf_energy2")
G.add_node("coupler2", vtype="mpc", label="Column N")
G.add_node("snow2", vtype="physical", label="")
G.add_edge("coupler2", "permafrost2")
G.add_edge("coupler2", "snow2")
G.add_edge("column_coupler", "coupler2")

writePDF(G, 'pk3', scale=2, layout_scaling=(6,4))

# NOTE: add the followign line to pk3.tex:
#      \draw[loosely dashed, line width=5] (coupler1) to (coupler2);
# just before the line, \end{tikzpicture}
