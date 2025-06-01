# REPORT

## PART 1
### 1. (0.5 points) Report on the order and size of both networks. Explain why these numbers may be different from the number of rows in the TF-RISet file. 
#### Order and Size:
- Without operons:
    - Nodes (order): 1969
    - Edges (size): 3834

- With operons:
    - Nodes (order): 3170
    - Edges (size): 7800

This differs from the numbers of rows in the TF-RISet file because some rows encode interactions with multiple genes or operons.


### 2. (1.5 points) Report the number of TFs and TGs in both networks, and the fraction of genes in the E. coli network that are captured as being involved in transcriptional regulatory interactions in the RegulonDB dataset. What does this tell us about the genomic organization of the E. coli transcriptional regulatory network? 
#### Ecoli TRN Analysis

Number of TFs: 241, Number of TGs: 1728
Total genes: 4651, Involved genes: 1969
Fraction of genes involved in TRN: 0.42


#### Ecoli operon TRN Analysis 
Number of TFs: 241, Number of TGs: 2929
Total genes: 4651, Involved genes: 3170
Fraction of genes involved in operon TRN: 0.68

With operon detection, genes involved in regulatori interactions rises to 68%, even though the number of TFs stays the same.
This reflects the regulation of E. coli’s genes often targets operons, not individual genes. A single TF can control multiple genes indirectly by regulating a promoter upstream of an operon, which transcribes multiple genes as one unit.


### 3. (1 points) Report on the nodes with maximum in- and out-degree (and their respective in- and out-degrees). Are they the same on both networks? Why?
#### MAX DEGREES
- Ecoli TRN: Max in-degree node: csgD (b1040) = 16, Max out-degree node: nac (b1988) = 521
- Ecoli operon TRN: Max in-degree node: csgD (b1040) = 16, Max out-degree node: nac (b1988) = 1039

The most regulated and most regulatory genes are the same in both networks, this can be because detection affects only target genes, but the number of transcription factor doesn't change. 
- With operons, nac appears to regulate more genes, so it remains the node with the highest out-degree.
- The most targeted gene, csgD, stays the same because operon expansion doesn’t add new TFs regulating it,  it only increases the number of targets each TF regulates.

## Part 3  

### 1. Compare the abundance of FFL and BF motifs between the gene-level and operon-level networks. Which network has more motifs of each type? What might explain the differences observed?
We compared the abundance of feed-forward loop FFL and bi-fan BF motifs in both the gene-level and operon-level transcriptional regulatory networks (TRNs) of E. coli.  
Motif counts:  
- Gene-level TRN:  
- FFL: 874  
- BF: 17,703  
Operon-level TRN:  
- FFL: 2,953  
- BF: 117,992  
Despite the operon-level network having fewer nodes (due to gene grouping), it displays significantly more motifs of both types. This can be explained by the aggregation of genes into operons, which increases the density of regulatory connections. When multiple TFs regulate a common operon, or when operons themselves are targets and intermediates in regulatory paths, motifs are more likely to form. The grouping effect effectively increases the number of overlapping or shared targets, especially enhancing the formation of bi-fan motifs.

### 2.Report on the distribution of motif cluster sizes for both motif types and both networks. What do you observe? Are motif clusters typically small or large? Does this align with the results you obtained when analyzing the overall properties of the E. coli transcriptional regulatory network in session 2? 

We analyzed the sizes of motif clusters (connected components formed by edges participating in motifs).  
FFL Clusters:  
- Gene-level: Two clusters of size 4 and 438 (avg: 221)  
- Operon-level: Six clusters ranging from 3 to 1101 nodes (avg: ~161)  
BF Clusters:  
- Gene-level: Clusters of size 4 and 775 (avg: ~390)  
- Operon-level: Clusters of size 4, 10, and 1764 (avg: ~593)  
These results show that motifs do not just form isolated modules but aggregate into large interconnected regions, particularly bi-fans. This is consistent with the overall TRN structure observed in earlier sessions, where we noted small-world properties and the existence of large hub-regulated subgraphs.

### 3. Are there any coincidences between any secondary clusters you obtain for FFL and BF in the E. coli (no operon) network? Look up the genes involved to see if you can explain why these genes are set apart from the rest of the network. 

We identified overlapping genes in secondary (non-main) clusters of both FFL and BF motifs in the gene-level TRN. One such set includes: b1569, b1570, b1575, b4135  
These genes are part of a cohesive subnetwork likely dedicated to a specific cellular response. The fact that they form overlapping FFL and BF motif clusters but are isolated from the main network cluster suggests functional modularity; graph-theoretic modularity that aligns with localized computation or control logic within the network.

### 4. Compare the motif statistics of the E. coli networks to those of a randomized network with the same degree sequence (e.g., using the configuration model). Are motifs still overrepresented? What does this tell us about the evolutionary design of the TRN? How do operons contribute to this? 

We compared the real TRNs to randomized counterparts (preserving in-/out-degree distributions).  

Gene-level TRN:  
- FFL: 874 (real) ; 555.2 ± 58.3 (random)  
- BF: 17,703 (real) ; 15,147.4 ± 765 (random)  
Operon-level TRN:  
- FFL: 2,953 (real) ; 1,588.2 ± 256.5 (random)  
- BF: 117,992 (real) ; 78,641.4 ± 4,548.3 (random) 

Motifs are significantly overrepresented in real TRNs compared to randomized ones. This suggests that motif structures are not a random outcome of the degree distribution, but a product of selection; reflecting computational design principles like input integration (BF) and dynamic filtering (FFL). The operon abstraction magnifies this further, indicating that gene grouping supports denser regulatory logic.

### 5. araB ('b0063') gene is the lead gene of the araBAD operon, which encodes three proteins involved in the import and degradation of arabinose (a sugar used by E. coli as a secondary energy source (the primary one being glucose)). As in the case of the lac operon, the araBAD operon is maximally expressed when there is there is arabinose, but no glucose, around. That is, its promoter implements a (NOT(glucose) AND arabinose) logic. Analyze the interactions of the araB gene in the no-operon E. coli network. What TFs interact to regulate araB? Look them up. Do any of them sense presence/absence of glucose/arabinose? If so, how are they connected (what connections and what type (check the original RegulonDB file to see if they are activation/repression))? Do they make up a feed-forward loop with araB?

the original RegulonDB file to see if they are activation/repression))? Do they make up a feed-forward loop with araB?
We explored the local network around gene b0063 (araB) in the gene-level TRN.  

Regulators of araB:  
b0064 (AraC), b3357 (CRP), b0683 (Fnr) 

These form two coherent FFL motifs:  
CRP → AraC → araB and CRP → araB  
CRP → Fnr → araB and CRP → araB  
These motifs represent a regulatory structure where a primary controller (CRP) integrates upstream signals and passes them through intermediate TFs. This structure allows fine-tuned regulation of araB based on environmental signals.

### 6. When TFs sense an environmental input, their reaction is very fast and in a matter of seconds they start activating their target genes in response. As they do so, transcription of target genes is irrevocably initiated, leading to translation and synthesis of the resulting proteins after 10-20 minutes. Protein synthesis has an obvious energetic cost for the cell. As they swim around, E. coli cells may encounter sudden, temporary shifts in glucose levels. Explain how a feed-forward loop regulating the araBAD operon might prevent E. coli from spending energy by synthesizing arabinose degradation genes when arabinose and glucose are around and the cell experiences a short, downward shift in the concentration of glucose. How does that compare to a situation in which both TFs independently regulate araBAD?
The FFL architecture involving CRP and AraC implements signal filtering logic. In this structure:  
CRP is activated in the absence of glucose  
AraC is activated by arabinose

Both are needed to activate araB. If glucose briefly disappears, CRP activates temporarily, but AraC activation may lag. This introduces a delay, preventing premature or unnecessary expression of arabinose-degrading genes. Compared to independent regulation (no FFL), this architecture avoids costly synthesis in response to transient changes, providing a persistence filter.

### 7. When TFs sense an environmental input, their reaction is very fast and in a matter of seconds they start activating their target genes in response. As they do so, transcription of target genes is irrevocably initiated, leading to translation and synthesis of the resulting proteins after 10-20 minutes. Protein synthesis has an obvious energetic cost for the cell. As they swim around, E. coli cells may encounter sudden, temporary shifts in glucose levels. Explain how a feed-forward loop regulating the araBAD operon might prevent E. coli from spending energy by synthesizing arabinose degradation genes when arabinose and glucose are around and the cell experiences a short, downward shift in the concentration of glucose. How does that compare to a situation in which both TFs independently regulate araBAD?

Top 3 TFs (by out-degree): 
b3357 (CRP): 5,844 motifs
b1988 (FNR): 4,987 motifs
b0889 (IHF): 4,122 motifs 

Shared motifs: 
b1988 & b0889: 1,596
b1988 & b3357: 1,035
b0889 & b3357: 496  
This extensive overlap indicates that these TFs frequently co-regulate targets in bi-fan structures. Rather than forming long chains of sequential regulation (deep networks), E. coli appears to implement regulatory logic in a shallow topology, relying on local combinatorial control. This is consistent with a modular, efficient design strategy in compact genomes, optimizing for fast and reliable decision-making at promoter regions.



## PART 4
### Visualization Design in Gephi:
- Layout algorithm: Yifan Hu for spatial separation of dense hubs and peripheral nodes.
- Node size:scale by out-degree → larger nodes = more regulatory influence.
- Node color:
    - Color by ntype: Red for TFs, blue for TGs and yellow for autoregulatry TFs (separate regulators from targets).
- Edge width: edges had not weight here so we just display them gray.
- Labeling: Display of labels for top 10 nodes by out-degree (locus tag).

### Network Characteristics Conveyed:
- Hierarchy: Layout and node size show central TFs.
- Modularity: Local clusters can reveal functionally co-regulated genes.
- Regulatory influence: Out-degree highlights global vs. local regulators.
- Vulnerability: Central nodes represent points of control.