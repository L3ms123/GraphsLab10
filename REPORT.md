# REPORT

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
<<<<<<< HEAD
=======
#### MAX DEGREES
- Ecoli TRN: Max in-degree node: csgD (b1040) = 16, Max out-degree node: nac (b1988) = 521
- Ecoli operon TRN: Max in-degree node: csgD (b1040) = 16, Max out-degree node: nac (b1988) = 1039

The most regulated and most regulatory genes are the same in both networks, this can be because detection affects only target genes, but the number of transcription factor doesn't change. 
- With operons, nac appears to regulate more genes, so it remains the node with the highest out-degree.
- The most targeted gene, csgD, stays the same because operon expansion doesn’t add new TFs regulating it,  it only increases the number of targets each TF regulates.

>>>>>>> 4e43c1a69722101ab68508271960fba0e4ea7b76
