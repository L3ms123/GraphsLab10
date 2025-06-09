# -*- coding: utf-8 -*-
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import csv
import networkx as nx
import sys
import os
import matplotlib.pyplot as plt

"""
>>> feature = record.features[0]
>>> feature.type
'gene'
>>> feature.location
[100:500](+)
>>> feature.qualifiers
{'gene': ['lacZ'], 'locus_tag': ['b0001'], 'product': ['beta-galactosidase']}
"""

# ------- IMPLEMENT HERE ANY AUXILIARY FUNCTIONS NEEDED ------- #

def intergenic_distance(gene1, gene2, strand):
    if strand == 1:  # Forward strand
        # Distance from end of gene1 to start of gene2
        return gene2.location.start - gene1.location.end
    else:  # Reverse strand (strand == -1)
        # Distance from end of gene2 to start of gene1 (since we're going backwards)
        return gene1.location.start - gene2.location.end
    

def analyze_graph(G, genome):
    TFs = [n for n, d in G.nodes(data=True) if d.get("ntype") == "TF"]
    TGs = [n for n, d in G.nodes(data=True) if d.get("ntype") == "TG"]

    num_TFs = len(TFs)
    num_TGs = len(set(TGs))  # avoid duplicates
    total_genes = sum(1 for f in genome.features if f.type == "gene")
    involved_genes = len(set(TFs + TGs))
    fraction_involved = involved_genes / total_genes

    return num_TFs, num_TGs, total_genes, involved_genes, fraction_involved

def max_degrees(G):
    in_deg_node = max(G.in_degree(), key=lambda x: x[1])
    out_deg_node = max(G.out_degree(), key=lambda x: x[1])
    return in_deg_node, out_deg_node

# --------------- END OF AUXILIARY FUNCTIONS ------------------ #

def feature_list(genome : SeqRecord, query : str) -> list:
    '''
    Extract CDS (Coding Sequence) features with specific feature description.
    - param genome : SeqRecord
        genome SeqRecord object to be analyzed.
    - param query : str
        feature descriptor.
    - return list
        list of tuples (locus_tag, protein_id) matching descriptor.
    '''
    records = []

    features = [f for f in genome.features if f.type == 'CDS']

    for idx, feature in enumerate(features):
        if 'product' in feature.qualifiers and query in feature.qualifiers['product']:
            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers['locus_tag'][0]
            else:
                locus_tag = ''
            if 'protein_id' in feature.qualifiers:
                protein_id = feature.qualifiers['protein_id'][0]
            else:
                protein_id = ''
            records.append((locus_tag, protein_id))
    return records




def gene_qualifier(query : str, query_field : str,
                   target_field : str, genome : SeqRecord) -> tuple:
    '''
    Obtain the specified qualifier identifier for a given gene qualifier
    - param: query: str
        gene name/locus_tag to map to corresponding locus_tag/name
    - param: query_field: str
        query type indicator (gene/locus_tag/protein_id/product)
    - param: target_field: str
        target type indicator (gene/locus_tag/protein_id/product)        
    - param: genome : SeqRecord
        genome SeqRecord object containing features
    - return: tuple
        int : feature index
        str : specified qualifier for gene (empty string if no match)
    '''
    features = [f for f in genome.features if f.type == 'gene']

    for idx, feature in enumerate(features):
        if query_field in feature.qualifiers and query in feature.qualifiers[query_field]:
            if target_field in feature.qualifiers:
                return idx, feature.qualifiers[target_field][0]
            else:
                return idx, ''
    return -1, ''


def operon(locus_tag: str, max_intergenic_dist: int, genome: SeqRecord) -> list:
    '''
    Obtain the locus_tag identifier for a given gene name
    - param: locus_tag : str locus_tag of lead operon gene
    - param: max_intergenic_dist : int maximum distance between consecutive, same strand genes
    - param: genome : SeqRecord genome SeqRecord object for operon inference
    - return: list list of locus_tags conforming operon (including query)
    '''
    
    # Find the starting gene by locus_tag
    start_gene = None
    for feature in genome.features:
        if feature.type == 'gene':
            if 'locus_tag' in feature.qualifiers:
                if feature.qualifiers['locus_tag'][0] == locus_tag:
                    start_gene = feature
                    break
    
    if start_gene is None:
        return []  
    
    # Get all genes sorted by position
    all_genes = []
    for feature in genome.features:
        if feature.type == 'gene':
            if 'locus_tag' in feature.qualifiers:
                all_genes.append(feature)
    
    # Sort by start position
    all_genes.sort(key=lambda x: x.location.start)
    
    # Find the index starting gene
    start_idx = None
    for i, gene in enumerate(all_genes):
        if gene.qualifiers['locus_tag'][0] == locus_tag:
            start_idx = i
            break
    
    if start_idx is None:
        return []
    
    operon_genes = [start_gene]
    start_strand = start_gene.location.strand
    
    # Extend operon in the forward direction 
    if start_strand == 1:  # Forward strand genes
        current_idx = start_idx
        while current_idx + 1 < len(all_genes):
            next_gene = all_genes[current_idx + 1]
            # Check if next gene is on the same strand
            if next_gene.location.strand != start_strand:
                break
            
            # Calculate intergenic distance
            current_gene = all_genes[current_idx]
            dist = intergenic_distance(current_gene, next_gene, start_strand)
            
            # Check if distance is within threshold
            if dist <= max_intergenic_dist:
                operon_genes.append(next_gene)
                current_idx += 1
            else:
                break
    
    else:  # Reverse strand genes
        current_idx = start_idx
        while current_idx - 1 >= 0:
            prev_gene = all_genes[current_idx - 1]
            # Check if previous gene is on the same strand
            if prev_gene.location.strand != start_strand:
                break
            
            # Calculate intergenic distance
            current_gene = all_genes[current_idx]
            dist = intergenic_distance(current_gene, prev_gene, start_strand)
            
            # Check if distance is within threshold
            if dist <= max_intergenic_dist:
                operon_genes.insert(0, prev_gene)  # Insert at beginning to maintain order
                current_idx -= 1
            else:
                break
    
    # Extend operon in the backward direction 
    if start_strand == 1:  # Forward strand genes
        current_idx = start_idx
        while current_idx - 1 >= 0:
            prev_gene = all_genes[current_idx - 1]
            # Check if previous gene is on the same strand
            if prev_gene.location.strand != start_strand:
                break
            
            # Calculate intergenic distance
            current_gene = all_genes[current_idx]
            dist = intergenic_distance(prev_gene, current_gene, start_strand)
            
            # Check if distance is within threshold
            if dist <= max_intergenic_dist:
                operon_genes.insert(0, prev_gene)  # Insert at beginning
                current_idx -= 1
            else:
                break
    
    else:  # Reverse strand genes 
        current_idx = start_idx
        while current_idx + 1 < len(all_genes):
            next_gene = all_genes[current_idx + 1]
            # Check if next gene is on the same strand
            if next_gene.location.strand != start_strand:
                break
            
            # Calculate intergenic distance
            current_gene = all_genes[current_idx]
            dist = intergenic_distance(next_gene, current_gene, start_strand)
            
            # Check if distance is within threshold
            if dist <= max_intergenic_dist:
                operon_genes.append(next_gene)
                current_idx += 1
            else:
                break
    
    # Extract locus_tags
    locus_tags = []
    for gene in operon_genes:
        if 'locus_tag' in gene.qualifiers:
            locus_tags.append(gene.qualifiers['locus_tag'][0])
    
    # Remove duplicates
    seen = set()
    unique_locus_tags = []
    for tag in locus_tags:
        if tag not in seen:
            seen.add(tag)
            unique_locus_tags.append(tag)
    
    return unique_locus_tags



def TF_RISet_parse(tf_riset_filename: str, tf_set_filename: str,
                   detect_operons: bool, max_intergenic_dist: int,
                   genome: SeqRecord) -> nx.DiGraph:
    """
    Parse TF-RISet file to obtain a TRN graph.
    The TFSet file will be used to extract information on the gene coding for
    each transcription factor [4)geneCodingForTF, 5)geneBnumberCodingForTF]
    The TF-RISet file will be used to extract genes regulated by each TF.
    We will create nodes using their locus_tag identifier, and save the gene
    name as a 'name' attribute. 
    If selected, operons will be predicted for each of these genes to determine 
    the entire set of genes regulated by the TF.
    
    - param: tf_riset_filename : str
        name/path of the TF_RISet file
    - param: tf_set_filename : str
        name/path of the TFSet file
    - param: detect_operons : bool
        whether we run operon detection
    - param: max_intergenic_dist : int
        maximum distance between consecutive, same strand genes
    - param: genome : SeqRecord
        genome SeqRecord object to extract information from
    - return: TF dictionary
    """    
    locus_to_name = {}
    name_to_locus = {}
    for feat in genome.features:
        locus = feat.qualifiers.get("locus_tag", [""])[0].strip()
        name = feat.qualifiers.get("gene", [""])[0].strip()
        if locus:
            locus_to_name[locus] = name
        if name:
            name_to_locus[name] = locus

    # Parse TFSet: guarda todo
    TF_dict = {}
    with open(tf_set_filename, "r") as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for cols in reader:
            if len(cols) < 5:
                continue
            tf_name = cols[1].strip()
            gene_name = cols[3].strip()
            locus_tag = cols[4].strip() or name_to_locus.get(gene_name, "")
            if not locus_tag and gene_name in name_to_locus:
                locus_tag = name_to_locus[gene_name]
            TF_dict[tf_name] = (gene_name, locus_tag)

    G = nx.DiGraph()
    skipped_tg = 0

    with open(tf_riset_filename, "r") as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for cols in reader:
            if len(cols) < 19:
                continue

            tf_name = cols[3].strip()
            if tf_name not in TF_dict:
                continue
            tf_gene_name, tf_id = TF_dict[tf_name]
            G.add_node(tf_id, ntype="TF", name=tf_gene_name)

            target_info = cols[18].strip()
            tg_gene_name = target_info.split(":")[-1].strip()

            tg_ids = []

            if '-' in tg_gene_name:
                parts = tg_gene_name.split('-')
                for part in parts:
                    tg_id = name_to_locus.get(part, "")
                    if tg_id:
                        tg_ids.append((tg_id, part))
            else:
                tg_id = name_to_locus.get(tg_gene_name, "")
                if tg_id:
                    tg_ids.append((tg_id, tg_gene_name))
                else:
                    # Intentar descomponer en subgenes  como"gadAX"
                    if len(tg_gene_name) > 4:
                        prefix = tg_gene_name[:3]
                        suffixes = tg_gene_name[3:]
                        for s in suffixes:
                            subgene = prefix + s
                            if subgene in name_to_locus:
                                tg_id = name_to_locus[subgene]
                                tg_ids.append((tg_id, subgene))
                                break

            if not tg_ids:
                skipped_tg += 1
                # print(f"[SKIPPED TG] No se resolvió tg_id para: '{tg_gene_name}'")
                # print(f" - Línea original: {cols}")
                # print(f" - TF: {tf_name}, TF_gene_name: {tf_gene_name}, TF_id: {tf_id}")
                continue

            for tg_id, tg_name in tg_ids: 
                if tg_id not in G:
                    G.add_node(tg_id, ntype="TG", name=locus_to_name.get(tg_id, tg_name))
                G.add_edge(tf_id, tg_id)


                if detect_operons and tg_id in locus_to_name:
                    try:
                        op_genes = operon(tg_id, max_intergenic_dist, genome)
                        for op_gene in op_genes:
                            if op_gene == tg_id:
                                continue
                            op_name = locus_to_name.get(op_gene, op_gene)
                            if op_gene not in G:
                                G.add_node(op_gene, ntype="TG", name=op_name)
                            G.add_edge(tf_id, op_gene)
                    except Exception as e:
                        print(f"Operon detection failed for {tg_id}: {e}")

    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges. Skipped {skipped_tg} TGs.")
    return G


if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

    gb_path = os.path.join(base_dir, "data", "U00096.3.gb")
    tf_riset_path = os.path.join(base_dir, "data", "TF-RISet.tsv")
    tf_set_path = os.path.join(base_dir, "data", "TFSet.tsv")
    EColi = SeqIO.read(gb_path, "genbank")
    Ecoli_TRN = TF_RISet_parse(tf_riset_path, \
                              tf_set_path, \
                              detect_operons=False, \
                              max_intergenic_dist=100, \
                              genome=EColi)
    
    Ecoli_operon_TRN = TF_RISet_parse(tf_riset_path, \
                              tf_set_path, \
                              detect_operons=True, \
                              max_intergenic_dist=100, \
                              genome=EColi)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "output_graphs")

    nx.write_graphml(Ecoli_TRN, os.path.join(output_dir, "Ecoli_TRN.graphml"))
    nx.write_graphml(Ecoli_operon_TRN, os.path.join(output_dir, "Ecoli_operon_TRN.graphml")) 

    nx.draw(Ecoli_operon_TRN, with_labels=True, node_size=300, node_color="lightblue")
    plt.show()

    tf1, tg1, total1, involved1, frac1 = analyze_graph(Ecoli_TRN, EColi)
    tf2, tg2, total2, involved2, frac2 = analyze_graph(Ecoli_operon_TRN, EColi)
    print('TRANSCRIPTION REGULATORY NETWORK ANALYSIS')
    print('--------------- Ecoli TRN Analysis --------------')
    print(f"Number of TFs: {tf1}, Number of TGs: {tg1}")
    print(f"Total genes: {total1}, Involved genes: {involved1}")
    print(f"Fraction of genes involved in TRN: {frac1:.2f}")

    print('------------- Ecoli operon TRN Analysis ----------')
    print(f"Number of TFs: {tf2}, Number of TGs: {tg2}")
    print(f"Total genes: {total2}, Involved genes: {involved2}")
    print(f"Fraction of genes involved in operon TRN: {frac2:.2f}")

    print('----------------- MAX DEGREES --------------------')
    in1, out1 = max_degrees(Ecoli_TRN)
    in2, out2 = max_degrees(Ecoli_operon_TRN)
    in_name1 = Ecoli_TRN.nodes[in1[0]].get("name", "unknown")
    out_name1 = Ecoli_TRN.nodes[out1[0]].get("name", "unknown")
    in_name2 = Ecoli_operon_TRN.nodes[in2[0]].get("name", "unknown")
    out_name2 = Ecoli_operon_TRN.nodes[out2[0]].get("name", "unknown")
    print(f"Ecoli TRN: Max in-degree node: {in_name1} ({in1[0]}) = {in1[1]}, Max out-degree node: {out_name2} ({out1[0]}) = {out1[1]}")
    print(f"Ecoli operon TRN: Max in-degree node: {in_name2} ({in2[0]}) = {in2[1]}, Max out-degree node: {out_name2} ({out2[0]}) = {out2[1]}")
    print('--------------------------------------------------')