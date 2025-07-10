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

In `results/`, you can also find the subfolder `my_functions` that contains useful functions that are used several times across scripts. The script `blast_CDS.py` is also useful to blast any sequence in a fasta file against any other fasta file locally.


## All that has been done since writing the report (`report.pdf`)
- We added the species (based on the ANI) to some plots (`3_plot_dense_results/denovo_trg_116_with_species.png` and plots in `6_plot_integrity_results/clusters`)


## Ideas for the future
