# XIOS_RNA_fingerprint
```
This project contains the source code for the paper "Accurate Classification of RNA Structures Using 
Topological Fingerprints"(https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0164726), 
a graph-based RNA structure matching and classification project. 
```
# Abstract
```
While RNAs are well known to possess complex structures, functionally similar RNAs often have little 
sequence similarity. While the exact size and spacing of base-paired regions vary, functionally similar 
RNAs have pronounced similarity in the arrangement, or topology, of base-paired stems. Furthermore, 
predicted RNA structures often lack pseudoknots (a crucial aspect of biological activity), and are only 
partially correct, or incomplete. A topological approach addresses all of these difficulties. In this 
work we describe each RNA structure as a graph that can be converted to a topological spectrum (RNA 
fingerprint). The set of subgraphs in an RNA structure, its RNA fingerprint, can be compared with the 
fingerprints of other RNA structures to identify and correctly classify functionally related RNAs. 
Topologically similar RNAs can be identified even when a large fraction, up to 30%, of the stems are 
omitted, indicating that highly accurate structures are not necessary. We investigate the performance of 
the RNA fingerprint approach on a set of eight highly curated RNA families, with diverse sizes and 
functions, containing pseudoknots, and with little sequence similarity–an especially difficult test set. 
In spite of the difficult test set, the RNA fingerprint approach is very successful (ROC AUC > 0.95). 
Due to the inclusion of pseudoknots, the RNA fingerprint approach both covers a wider range of possible 
structures than methods based only on secondary structure, and its tolerance for incomplete structures 
suggests that it can be applied even to predicted structures.

```

# unzip the input RNA graph files
```
tar -xvzf inputs.tar.gz
```

# generate RNA fingerprint
```
perl fingerprint_shred.pl -i ./inputs/curated_graphs/xios_graph/ -w ./outputs/curated_graphs/fingerprint/
perl fingerprint_shred.pl -i ./inputs/expanded_graphs/xios_graph/ -w ./outpus/expanded_graphs/fingerprint/
perl fingerprint_shred.pl -i ./inputs/incomplete_graphs/remove_00_percent_stems/xios_graph/ -w ./outputs/incomplete_graphs/remove_00_percent_stems/fingerprint/
perl fingerprint_shred.pl -i ./inputs/incomplete_graphs/remove_10_percent_stems/xios_graph/ -w ./outputs/incomplete_graphs/remove_10_percent_stems/fingerprint/
perl fingerprint_shred.pl -i ./inputs/incomplete_graphs/remove_30_percent_stems/xios_graph/ -w ./outputs/incomplete_graphs/remove_30_percent_stems/fingerprint/
perl fingerprint_shred.pl -i ./inputs/incomplete_graphs/remove_50_percent_stems/xios_graph/ -w ./outputs/incomplete_graphs/remove_50_percent_stems/fingerprint/
perl fingerprint_shred.pl -i ./inputs/incomplete_graphs/remove_70_percent_stems/xios_graph/ -w ./outputs/incomplete_graphs/remove_70_percent_stems/fingerprint/
```

# calculate RNA fingerprint similarity 
```
perl fingerprint_distance.pl -d ./inputs/curated_graphs/fingerprint/ -w ./outputs/curated_graphs/clustering/
perl fingerprint_distance.pl -d ./inputs/expanded_graphs/fingerprint/ -w ./outputs/expanded_graphs/clustering/
perl fingerprint_distance.pl -d ./inputs/incomplete_graphs/remove_00_percent_stems/fingerprint/ -w ./outputs/incomplete_graphs/remove_00_percent_stems/clustering/
perl fingerprint_distance.pl -d ./inputs/incomplete_graphs/remove_10_percent_stems/fingerprint/ -w ./outputs/incomplete_graphs/remove_10_percent_stems/clustering/
perl fingerprint_distance.pl -d ./inputs/incomplete_graphs/remove_30_percent_stems/fingerprint/ -w ./outputs/incomplete_graphs/remove_30_percent_stems/clustering/
perl fingerprint_distance.pl -d ./inputs/incomplete_graphs/remove_50_percent_stems/fingerprint/ -w ./outputs/incomplete_graphs/remove_50_percent_stems/clustering/
perl fingerprint_distance.pl -d ./inputs/incomplete_graphs/remove_70_percent_stems/fingerprint/ -w ./outputs/incomplete_graphs/remove_70_percent_stems/clustering/
```

# curated graphs
```
A set of curated RNA structures have been collected from the literature and a variety of biological 
databases and is extended in this work. This set of known structures has been carefully selected to 
contain pseudoknots, to cover a broad range of lengths, and to have been the subject of extensive 
expert curation by the biological community. This curated set includes 206 structures of transfer RNA, 
Ribonuclease P RNA, transfer-messenger RNA, group I and group II self-splicing introns, and 5S, 16S 
and 23S ribosomal RNA. The structures in this curated set have been reviewed to ensure they reflect 
expert opinion on the correct structure, and to ensure that the reported structures are as accurate 
as possible given existing experimental data such as X-ray crystallography and covariance analysis. 
The curated structures have been screened to ensure that all structures are full-length, and no pair 
of structures has greater than 50% sequence identity. Multiple families of the curated structures 
contain pseudoknots. 
```

# expanded_graphs
```
We selected a set of 177 RNA graphs containing up to 25 vertices from the curated data set, and create 
an expanded set by embedding subgraphs, randomly selected according to probability of occurrence, from 
the decoy database into these RNA graphs until each RNA graph contained 30 vertices. 
As a control, we also created a decoy set of 45 graphs with 30 vertices, by random embedding of subgraphs 
from the decoy database only, i.e., graphs with no information from real biological structures except the 
frequency of occurrence of subgraphs in the known structures. 
Both the expanded and the decoy graph sets are completely free of size effects since they all have exactly 
the same number of stems. 
```

# incomplete_graphs
```
To test the effects of graph incompleteness on the extended-fingerprint Jaccard Similarity function, 
incomplete RNA graphs were generated by randomly removing a percentage (10%, 30%, 50%, and 70%, 
respectively) of the vertices (stems) in the curated structures. 
Note: There are 1 RNA graph with 50% vertice removal and 19 RNA graphs with 70% vertice removal that 
became fully unconnected graphs (none of the vertices connect), which were not used for fingerprint 
calculation. 
```
