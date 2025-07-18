All scripts for de novo emergence analysis are comprised in the `/scripts` folder. They are ordered in the following order:
1. Genome annotation
2. De novo identification
3. Plotting of the de novo results
4. Integrity analysis and identification of "good" de novo candidates (n = 64)
5. Clustering of the 64 de novo
6. Plotting of the integrity analysis results
7. De novo characterisation, using sequence and structure descriptors
8. De novo characterisation visualisation
9. Non-coding origin identification
10. Result generalisation, where we compare the 64 de novo to the 109 discarded candidates.


The corresponding results are in the `/results` folder.

In `scripts/`, you can also find the subfolder `my_functions` that contains useful functions that are used several times across scripts. The script `blast_CDS.py` is also useful to blast any sequence in a fasta file against any other fasta file locally.


## All that has been done since writing the report (`report.pdf`)
- We added the species (based on the ANI) to some plots (`3_plot_dense_results/denovo_trg_116_with_species.png` and plots in `6_plot_integrity_results/clusters`)
- We conducted a conservation analysis for all 64 de novo, to see how they are conserved accress genomes compared with non-coding ORFs (`11_conservation_analysis`)


## Ideas for the future
- Structure prediction
- In order to control whether the de novo genes are genes and not annotation errors, we can look at gene expression data, ie. RNAseq data
- In the same spirit, we can calculate the DNDS for de novo genes for which we have a multiple alignment (see `results/11_conservation_analysis/homologs/`)
- Re-evaluate whether 70% is a good coverage threshold for the integrity analysis (probably using the results of the conservation analysis `11_conservation_analysis`)
