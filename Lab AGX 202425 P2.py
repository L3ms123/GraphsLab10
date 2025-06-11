from collections import defaultdict, Counter
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from random import randint, sample
from scipy import stats
from itertools import combinations

# ------- IMPLEMENT HERE ANY AUXILIARY FUNCTIONS NEEDED ------- #
def plot_degree_distribution(degree_freq_list, title_prefix="Degree Distribution"):
    degrees = np.arange(len(degree_freq_list))
    freqs = np.array(degree_freq_list)
    
    # Linear scale plot
    plt.figure(figsize=(12,5))
    plt.subplot(1,2,1)
    plt.bar(degrees, freqs, color='skyblue')
    plt.title(f"{title_prefix} (Linear scale)")
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    
    # Log-log scale plot (avoid zero frequencies)
    nonzero = freqs > 0
    plt.subplot(1,2,2)
    plt.scatter(degrees[nonzero], freqs[nonzero], color='red')
    plt.xscale('log')
    plt.yscale('log')
    plt.title(f"{title_prefix} (Log-Log scale)")
    plt.xlabel("Degree (log scale)")
    plt.ylabel("Frequency (log scale)")
    plt.tight_layout()
    plt.show()
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
        Type of nodes on which we report degrees ('TF', 'TG', 'all')

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
    original_avg_distance = nx.average_shortest_path_length(G) # utilitzar el de netx

    results = {}
    groupings = list(combinations(node_list, grouping_size))

    for grouping in groupings:
        G_copy = G.copy()
        G_copy.remove_nodes_from(grouping)
        


        new_avg_distance = average_distance(G_copy, iterations=iterations)
        diff = new_avg_distance - original_avg_distance
        results[grouping] = diff

    return results
    # ----------h-h------ END OF FUNCTION --------------------- #


def simulate_attack(G, attack_type='random', removal_fraction=0.2):
    H = G.copy()
    n_removals = int(removal_fraction * len(H))
    sizes = [len(H) / len(G)]      
    if attack_type == 'targeted':
        nodes_by_degree = sorted(H.degree, key=lambda x: x[1], reverse=True)
        removal_nodes = [node for node, _ in nodes_by_degree[:n_removals]]
    else: #when its random
        removal_nodes = sample(list(H.nodes()), n_removals)
    
    batch_size = min(100, max(1, n_removals // 10))
    for i in range(0, n_removals, batch_size):
        batch = removal_nodes[i:i+batch_size]
        H.remove_nodes_from(batch)
        lcc = max(nx.connected_components(H), key=len, default=[])
        sizes.append(len(lcc) / len(G))
    
    return np.linspace(0, removal_fraction, len(sizes)), sizes

if __name__ == "__main__":
    # ------- IMPLEMENT HERE THE MAIN FOR THIS SESSION ------- #
    print("Loading graph...")
    graph_file = "./output_graphs/Ecoli_TRN.graphml"
    G = nx.read_graphml(graph_file)
    print(f"Graph has: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    degree_dist, top_nodes = degree_distribution(G, degree_type="out", node_type="all", bestN=10)


    UG = G.to_undirected()
    LCC = largest_CC_graph(UG)

    
    tf_nodes = [n for n, attr in LCC.nodes(data=True) if attr.get('ntype') == 'TF']
    impact_dict = deletion_impact(LCC, tf_nodes, grouping_size=1, iterations=4100)
    for group, impact in sorted(impact_dict.items(), key=lambda x: -x[1])[:30]:
        print(f"  {group}: Δ = {impact:.4f}")
    top = [g[0] for g, _ in sorted(impact_dict.items(), key=lambda x: -x[1])[:30]]
    print("w")
    impact2 = deletion_impact(LCC, top, grouping_size=2, iterations=4100)
    for group, impact in sorted(impact2.items(), key=lambda x: -x[1])[:5]:
        print(f"  {group}: Δ = {impact:.4f}")


    
    n = len(LCC)
    m = LCC.number_of_edges()
    
    # ER model
    p = 2*m/(n*(n-1)) 
    G_er = nx.fast_gnp_random_graph(n, p, seed=42)
    
    # BA model
    m_ba = max(1, m // n) 
    G_ba = nx.barabasi_albert_graph(n, m_ba, seed=42)
    
    # Test 
    print("\nTesting random failure robustness...")
    for name, graph in [('E. coli', LCC), ('ER', G_er), ('BA', G_ba)]:
        fractions, sizes = simulate_attack(graph, 'random', 0.5)
        final_size = sizes[-1] if sizes else 0
        print(f"  {name}: Remaining LCC = {final_size*100:.1f}% after 50% random removal")
    
    print("\nTesting targeted attack vulnerability...")
    for name, graph in [('E. coli', LCC), ('ER', G_er), ('BA', G_ba)]:
        fractions, sizes = simulate_attack(graph, 'targeted', 0.2)
        final_size = sizes[-1] if sizes else 0
        print(f"  {name}: Remaining LCC = {final_size*100:.1f}% after 20% hub removal")
    
    print("\nAverage Path Lengths:")
    for name, graph in [('E. coli', LCC), ('ER', G_er), ('BA', G_ba)]:
        avg_path = nx.average_shortest_path_length(largest_CC_graph(graph))
        
        print(f"  {name}: {avg_path:.3f}")
    # ------------------- END OF MAIN ------------------------ #
