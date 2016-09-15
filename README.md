# XIOS_RNA_fingerprint
Source code for XIOS RNA fingerprint generation and comparison

# generate RNA fingerprint
perl fingerprint_shred.pl -i ./inputs/xios_graph/ -w ./inputs/xios_graph/

# calculate RNA fingerprint similarity 
fingerprint_distance.pl -d ./inputs/fingerprint/

