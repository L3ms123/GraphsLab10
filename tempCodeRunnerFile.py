
# import csv
# import networkx as nx
# from Bio.SeqRecord import SeqRecord

# def TF_RISet_parse(tf_riset_filename: str, tf_set_filename: str,
#                    detect_operons: bool, max_intergenic_dist: int,
#                    genome: SeqRecord) -> nx.DiGraph:
#     """
#     Parse TF-RISet file to obtain a TRN graph.
#     """
#     # --- Index genome features for fast lookup ---
#     # Map: locus_tag -> gene feature, gene name -> locus_tag
#     locus_to_name = {}
#     name_to_locus = {}
#     for feat in genome.features:
#         if feat.type == "gene":
#             locus_tag = feat.qualifiers.get("locus_tag", [""])[0]
#             gene_name = feat.qualifiers.get("gene", [""])[0]
#             if locus_tag:
#                 locus_to_name[locus_tag] = gene_name
#             if gene_name:
#                 name_to_locus[gene_name] = locus_tag

#     # --- Parse TFSet file to map TF name to (gene_name, locus_tag) ---
#     TF_dict = {}
#     with open(tf_set_filename, "r") as f:
#         reader = csv.reader(f, delimiter='\t')
#         next(reader)  # skip header
#         for cols in reader:
#             if len(cols) < 5:
#                 continue
#             tf_name = cols[1].strip()      # tfGeneName
#             gene_name = cols[3].strip()    # geneCodingForTF
#             locus_tag = cols[4].strip()    # geneBnumberCodingForTF or locus_tag
#             # Try to resolve locus_tag if it's a gene name
#             if locus_tag not in locus_to_name and locus_tag in name_to_locus:
#                 locus_tag = name_to_locus[locus_tag]
#             TF_dict[tf_name] = (gene_name, locus_tag)

#     G = nx.DiGraph()

#     # --- Parse TF-RISet file ---
#     with open(tf_riset_filename, "r") as f:
#         reader = csv.reader(f, delimiter='\t')
#         next(reader)  # skip header
#         for cols in reader:
#             if len(cols) < 19:
#                 continue
#             tf_name = cols[3].strip()            # regulatorName
#             target_info = cols[18].strip()       # targetTuOrGene: format id:name

#             if tf_name not in TF_dict:
#                 continue

#             tf_gene_name, tf_locus_tag = TF_dict[tf_name]
#             if not tf_locus_tag:
#                 continue

#             # Add TF node
#             if tf_locus_tag not in G:
#                 G.add_node(tf_locus_tag, ntype="TF", name=tf_gene_name)

#             # Parse target gene
#             if ':' in target_info:
#                 tg_id, tg_gene_name = target_info.split(':', 1)
#             else:
#                 tg_id = tg_gene_name = target_info

#             # Try to resolve locus_tag for target gene
#             tg_locus_tag = ""
#             if tg_id in locus_to_name:
#                 tg_locus_tag = tg_id
#                 tg_gene_name = locus_to_name[tg_id]
#             elif tg_gene_name in name_to_locus:
#                 tg_locus_tag = name_to_locus[tg_gene_name]
#             else:
#                 continue  # Cannot resolve target locus_tag

#             # Add TG node
#             if tg_locus_tag not in G:
#                 G.add_node(tg_locus_tag, ntype="TG", name=tg_gene_name)
#             elif G.nodes[tg_locus_tag].get("ntype") != "TF":
#                 G.nodes[tg_locus_tag]["ntype"] = "TG"

#             # Add edge TF -> TG
#             G.add_edge(tf_locus_tag, tg_locus_tag)

#             # --- Add operon genes if requested ---
#             if detect_operons:
#                 try:
#                     op_genes = operon(tg_locus_tag, max_intergenic_dist, genome)
#                     for op_gene in op_genes:
#                         if op_gene == tg_locus_tag:
#                             continue  # Already added
#                         op_gene_name = locus_to_name.get(op_gene, op_gene)
#                         if op_gene not in G:
#                             G.add_node(op_gene, ntype="TG", name=op_gene_name)
#                         elif G.nodes[op_gene].get("ntype") != "TF":
#                             G.nodes[op_gene]["ntype"] = "TG"
#                         G.add_edge(tf_locus_tag, op_gene)
#                 except Exception:
#                     continue

#     return G