import pandas as pd
import collections

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANKS_LETTERS = ["k", "p", "c", "o", "f", "g", "s"]
RANK = "superkingdom"

# Filenames
INPUT = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/reblast_output_8.tsv"
TAXIDS = "/home/eliott.tempez/Documents/archaea_data/re_blast_1000/ncbi_lineages_2025-02-17.csv"

RANK_LETTER = RANKS_LETTERS[RANKS.index(RANK)]
OUTPUT = f"/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/{RANK}_level.txt"

rank_count = collections.defaultdict(lambda: collections.defaultdict(int))
chunksize = 10 ** 5

# **1. Load Taxonomy Data Once**
print("Loading taxonomy data...")
taxid_to_rank = {}
with open(TAXIDS, "r") as f:
    for line in f:
        line_split = line.strip().split(",")
        taxid_to_rank[line_split[0]] = line_split[1]

# **2. Process BLAST Output in Chunks**
print("Processing BLAST output...")
for i, chunk in enumerate(pd.read_csv(INPUT, sep="\t", header=None, chunksize=chunksize)):
    print(f"Processing chunk {i + 1}...")

    chunk.columns = ["query", "subject", "taxids", "evalue", "qcov"]

    # **3. Group by Query and Count TaxIDs Efficiently**
    taxids_count = {}
    for query, sub_df in chunk.groupby("query"):
        taxid_list = [str(t).split(";")[0] for t in sub_df["taxids"]]
        taxids_count[query] = collections.Counter(taxid_list)

    # **4. Get Unique TaxIDs**
    unique_taxids = set(t for taxids in taxids_count.values() for t in taxids.keys())

    # **5. Retrieve Taxonomic Ranks in Bulk**
    ranks = {taxid: taxid_to_rank.get(taxid, "Unknown") for taxid in unique_taxids}

    # **6. Aggregate Counts at Superkingdom Level**
    for query, taxids in taxids_count.items():
        for taxid, count in taxids.items():
            rank = ranks[taxid]
            if rank == "":
                rank = "Unknown"
            rank_count[query][rank] += count


# **7. Save the Results**
print("Saving results...")
with open(OUTPUT, "w") as f:
    for query, ranks in rank_count.items():
        print(f"{query}\t" + ";".join(f"{rank}:{count}" for rank, count in ranks.items()))
        #f.write(f"{query}\t" + ";".join(f"{rank}:{count}" for rank, count in ranks.items()) + "\n")

print("Processing complete.")
