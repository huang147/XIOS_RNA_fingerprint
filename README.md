# XIOS_RNA_fingerprint
Source code for XIOS RNA fingerprint generation and comparison

# unzip the input RNA graph files
tar -xvzf expanded_graphs.tar.gz  

tar -xvzf incomplete_graphs.tar.gz

tar -xvzf curated_graphs.tar.gz  

# generate RNA fingerprint
perl fingerprint_shred.pl -i ./inputs/xios_graph/ -w ./inputs/xios_graph/

# calculate RNA fingerprint similarity 
perl fingerprint_distance.pl -d ./inputs/fingerprint/

