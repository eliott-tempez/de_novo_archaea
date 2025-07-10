This folder contains the necessary scripts to extract de novo genes, TRGs, CDSs and iORFs from all genomes, calculate all descriptors for all of these sequences, and statistically analyse the differences between the different kinds of sequences.

- extract_iorfs extracts all iorfs for all genomes
- compare_sequences.py and run_sequence_comparison.sh extract all sequences and apply all descriptors to it. It outputs the file sequence_features_good_candidates_all.csv
- calculate_pvals_seq_comparison.py and run_sequence_pvals.sh calculate the p-values for each sequence pair. It also permits the binning of the sequences by GC content.

So that they can be included in this repo, the corresponding results have been compressed.