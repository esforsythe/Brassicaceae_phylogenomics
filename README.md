# Brassicaceae_phylogenomics

The repository contains scripts used to infer gene trees from genes spanning nuclear and organellar genomes. Also included are scripts for performing various downstream analyses of gene tree topologies. 

Scripts are written in perl, shell, and R.

List of scripts included:

Clean_CDSandAAseqs.R
-An R script that cleans CDS files to contain only the longest isoform for each gene. This script also translates CDS to AA and output an AA fasta file to ensure that CDS and AA files are compatible.

batch_raxml.pl
-A perl wrapper script for running for inferring gene trees from multiple alignments 

Topology_analysis_nuc.R
-Script for assessing the topology of single-copy nuclear gene trees

Topology_analysis_ret_dups.R
-Script for splitting retained duplicate trees and assessing the topology of subtrees

Topology_analysis_CP.R
-Script for assessing the topology of CP trees

Topology_analysis_MT.R
-Script for assessing the topology of MT trees

Print_exact_tops.R
-Prints gene tree topologies with branchlengths and support values removed to make them more human readable        

D_stat_BCvsAB.R
-Script used to calculate the D-statistic comparing the frequencies of A(BC) and C(AB) for nuclear genes

D_stat_ACvsAB.R 
-Script used to calculate the D-statistic comparing the frequencies of B(AC) and C(AB) for nuclear genes

Node_depths.R
-Rate smooth gene trees and calculate node depths

IG_directionality.R
-Uses dS distances to perform a test of IG directionality (i.e. determine donor and recipient lineages)

GO_categories.R
-Extracts GO category data and calculated proportion of A(BC) and B(AC) with each GO category        
       
Prot_Interact_Network.R 
-Script used to perform a network analysis using protein-protein interaction data

Chromosome_maps.R
-Script to generate chromosome map figures

Linkage_blocks.R
-Script used to detect linear clustering of genes with the same topology along the C. rubella chromosome map.


Git notes:
git add --all
git status
git commit -m "notes about the changes made"
git push