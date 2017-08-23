# Brassicaceae_phylogenomics

The repository contains scripts used to infer gene trees from genes spanning nuclear and organellar genomes. Also included are scripts for performing various downstream analyses of gene tree topologies. 

Scripts are written in perl, shell, and R.

Git notes:

git add --all

git status

git commit -m "notes about the changes made"

git push


List of scripts included:
Clean_CDSandAAseqs.R
-An R script that cleans CDS files to contain only the longest isoform for each gene. This script also translates CDS to AA and output an AA fasta file to ensure that CDS and AA files are compatible.