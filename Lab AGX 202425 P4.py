# -*- coding: utf-8 -*-
import os
import networkx as nx  

# ------- IMPLEMENT HERE ANY AUXILIARY FUNCTIONS NEEDED ------- #
def annotate_nodes_for_gephi(G):
    for node in G.nodes():
        in_deg = G.in_degree(node)
        out_deg = G.out_degree(node)
        if out_deg > 0 and in_deg > 0:
            G.nodes[node]['type'] = 'TF-TG'
        elif in_deg > 0:
            G.nodes[node]['type'] = 'TG'
        else:
            G.nodes[node]['type'] = 'TF'
    return G
# --------------- END OF AUXILIARY FUNCTIONS ------------------ #


if __name__ == "__main__":
    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #
    path = os.path.dirname(os.path.abspath(__file__))
    G = nx.read_graphml(os.path.join(path, "output_graphs\Ecoli_TRN.graphml"))

    G = annotate_nodes_for_gephi(G)

    nx.write_graphml(G, os.path.join(path,"output_graphs\Ecoli_selfloop_TRN.graphml"))
    # ------------------- END OF MAIN ------------------------ #