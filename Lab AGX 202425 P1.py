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

    


def operon(locus_tag : str, max_intergenic_dist : int, genome : SeqRecord) -> list:
    '''
    Obtain the locus_tag identifier for a given gene name
    - param: locus_tag : str
        locus_tag of lead operon gene
    - param: max_intergenic_dist : int
        maximum distance between consecutive, same strand genes
    - param: genome : SeqRecord
        genome SeqRecord object for operon inference
    - return: list
        list of locus_tags conforming operon (including query)
    '''
    # Filter genes with locus_tag
    gene_features = [
        f for f in genome.features 
        if f.type == "gene" and "locus_tag" in f.qualifiers
    ]
    
    # Find the initial gene and its strand
    initial_gene = None
    for f in gene_features:
        if f.qualifiers["locus_tag"][0] == locus_tag:
            initial_gene = f
            break
    if not initial_gene:
        return []
    strand = initial_gene.location.strand

    # Sort all genes by start position
    gene_features.sort(key=lambda f: f.location.start)

    # Find the index of the initial gene after sorting
    try:
        idx_initial = gene_features.index(initial_gene)
    except ValueError:
        return []

    operon_genes = [locus_tag]

    # Detect upstream genes (reverse iteration)
    current_idx = idx_initial - 1
    while current_idx >= 0:
        prev_gene = gene_features[current_idx]
        if prev_gene.location.strand != strand:
            break
        # Calculate distance to previous gene
        curr_end = gene_features[current_idx + 1].location.start
        prev_end = prev_gene.location.end
        distance = abs(curr_end - prev_end)
        if distance <= max_intergenic_dist:
            operon_genes.insert(0, prev_gene.qualifiers["locus_tag"][0])
            current_idx -= 1
        else:
            break

    # Detect downstream genes (forward iteration)
    current_idx = idx_initial + 1
    while current_idx < len(gene_features):
        next_gene = gene_features[current_idx]
        if next_gene.location.strand != strand:
            break
        # Calculate distance to next gene
        curr_end = gene_features[current_idx - 1].location.end
        next_start = next_gene.location.start
        distance = abs(next_start - curr_end)
        if distance <= max_intergenic_dist:
            operon_genes.append(next_gene.qualifiers["locus_tag"][0])
            current_idx += 1
        else:
            break

    # Reverse order for reverse strand to maintain transcription direction
    if strand == -1:
        operon_genes = operon_genes[::-1]

    return operon_genes


   



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

    # Parse TFSet: no validación, guarda todo
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
                locus_tag = name_to_locus[tf_gene_name]
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
                    # Intentar descomponer en subgenes tipo "gadAX"
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
                print(f"[SKIPPED TG] No se resolvió tg_id para: '{tg_gene_name}'")
                print(f" - Línea original: {cols}")
                print(f" - TF: {tf_name}, TF_gene_name: {tf_gene_name}, TF_id: {tf_id}")
                continue

            for tg_id, tg_name in tg_ids:
                G.add_node(tg_id, ntype="TG", name=locus_to_name.get(tg_id, tg_name))
                G.add_edge(tf_id, tg_id)

                if detect_operons and tg_id in locus_to_name:
                    try:
                        op_genes = operon(tg_id, max_intergenic_dist, genome)
                        for op_gene in op_genes:
                            if op_gene == tg_id:
                                continue
                            op_name = locus_to_name.get(op_gene, op_gene)
                            G.add_node(op_gene, ntype="TG", name=op_name)
                            G.add_edge(tf_id, op_gene)
                    except Exception as e:
                        print(f"Operon detection failed for {tg_id}: {e}")

    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges. Skipped {skipped_tg} TGs.")
    return G


if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

    gb_path = os.path.join(base_dir, "data", "U00096.3.gb")
    tf_riset_path = os.path.join(base_dir, "data", "miniTF-RISet.tsv")
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
    
    
    nx.write_graphml(Ecoli_TRN, "Ecoli_TRN.graphml")
    nx.write_graphml(Ecoli_operon_TRN, "Ecoli_operon_TRN.graphml")  
    nx.draw(Ecoli_operon_TRN, with_labels=True, node_size=300, node_color="lightblue")
    plt.show()