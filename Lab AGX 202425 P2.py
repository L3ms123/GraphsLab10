# -*- coding: utf-8 -*-

from collections import defaultdict, Counter
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from random import randint, sample
from scipy import stats
from itertools import combinations

# ------- IMPLEMENT HERE ANY AUXILIARY FUNCTIONS NEEDED ------- #

# --------------- END OF AUXILIARY FUNCTIONS ------------------ #


def degree_distribution(G : nx.DiGraph, degree_type : str, \
                        node_type : str, bestN : int) -> tuple:
    '''Compute the in-, out- or general degree distribution.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - param: degree_type : str
        Type of degree to compute ('in', 'out', 'general')
    - param: node_type : str
        Type of nodes on which we report degrees ('TG', 'TG', 'all')

    - return: list
        list with degree frequency distribution
    - return: dict
        dictionary with node names as keys and their degrees as values
        for the top N nodes in the degree distribution
    '''
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    # Validate inputs
    if degree_type not in {'in', 'out', 'general'}:
        raise ValueError("degree_type must be one of: 'in', 'out', 'general'")
    if node_type not in {'TF', 'TG', 'all'}:
        raise ValueError("node_type must be one of: 'TF', 'TG', 'all'")

    # Filter nodes based on node_type
    if node_type == 'all':
        nodes_to_consider = G.nodes()
    else:
        nodes_to_consider = [n for n, attr in G.nodes(data=True) if attr.get('ntype') == node_type]

    if degree_type == 'in':  degrees = G.in_degree(nodes_to_consider)
    elif degree_type=='out':  degrees = G.out_degree(nodes_to_consider)
    elif degree_type=='general':  degrees = G.degree(nodes_to_consider)


    degree_dict = dict(degrees)

    # Degree frequency distribution
    degree_counts = Counter(degree_dict.values())
    max_degree = max(degree_counts) if degree_counts else 0
    degree_freq_list = [degree_counts.get(i, 0) for i in range(max_degree + 1)]

    # Top N nodes by degree
    topN_nodes = dict(sorted(degree_dict.items(), key=lambda item: item[1], reverse=True)[:bestN])

    return degree_freq_list, topN_nodes
    # ----------------- END OF FUNCTION --------------------- #


def largest_CC_graph(G : nx.Graph) -> nx.Graph:
    '''Generate the largest connected component graph.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - return: Networkx graph
        The graph corresponding to the largest connected component in G
    '''
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    ccs = list(nx.connected_components(G))
    largest_cc_nodes = max(ccs, key=len)
    largest_cc_subgraph = G.subgraph(largest_cc_nodes).copy()
    return largest_cc_subgraph
    # ----------------- END OF FUNCTION --------------------- #


def average_distance(G : nx.Graph, iterations : int) -> float:
    '''Estimate the average distance in the input graph.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - param: iterations : int
        Number of iterations to perform
    - return: float
        The estimated average distance in G
    '''
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    nodes = list(G.nodes())
    total = 0
    valid_samples = 0

    for _ in range(iterations):
        u, v = sample(nodes, 2)
        if nx.has_path(G, u, v):
            dist = nx.shortest_path_length(G, source=u, target=v)
            total += dist
            valid_samples += 1

    if valid_samples == 0:    return float('inf')

    return total / valid_samples
    # ----------------- END OF FUNCTION --------------------- #


def deletion_impact(G : nx.Graph, node_list : list,\
                      grouping_size : int, iterations : int) -> dict:
    '''Assess the impact of node deletions on the graph average distance.

    Parameters
    ----------
    - param: G : Networkx graph
        Graph to analyze
    - param: node_list : list
        List of nodes to delete from the network
    - param: grouping_size : list
        The size of the groupings among nodes in the list
    - param: iterations : int
        Number of iterations to perform for average distance
    - return: dict
        Dictionary with grouping node names tuples as keys and differential
        average distance as values.
    '''
    # ------- IMPLEMENT HERE THE BODY OF THE FUNCTION ------- #
    original_avg_distance = average_distance(G, iterations)

    results = {}
    groupings = list(combinations(node_list, grouping_size))

    for grouping in groupings:
        G_copy = G.copy()
        G_copy.remove_nodes_from(grouping)

        new_avg_distance = average_distance(G_copy, iterations)
        diff = new_avg_distance - original_avg_distance
        results[grouping] = diff

    return results
    # ----------h-h------ END OF FUNCTION --------------------- #



if __name__ == "__main__":
    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #

    print("Loading graph...")
    graph_file = "./output_graphs/Ecoli_operon_TRN.graphml"
    G = nx.read_graphml(graph_file)

    print(f"Graph has: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    degree_dist, top_tf_nodes = degree_distribution(G, degree_type="out", node_type="TF", bestN=5)
    print("\nTop 5 TFs by out-degree:")
    for node, degree in top_tf_nodes.items():
        print(f"  {node}: {degree}")

    UG = G.to_undirected()
    LCC = largest_CC_graph(UG)
    print(f"\nBiggest cc has: {LCC.number_of_nodes()} nodes, {LCC.number_of_edges()} edges")


    avg_dist = average_distance(LCC, iterations=500)
    print(f"\nAverage distance: {avg_dist:.4f}")


    tf_nodes = [n for n, attr in G.nodes(data=True) if attr.get('ntype') == 'TF']
    impact_dict = deletion_impact(UG, tf_nodes, grouping_size=1, iterations=300)

    print("\nImpact of removing TFs individually (top 5 by increase in avg dist):")
    for group, impact in sorted(impact_dict.items(), key=lambda x: -x[1])[:5]:
        print(f"  {group}: Δ = {impact:.4f}")
    # ------------------- END OF MAIN ------------------------ #
