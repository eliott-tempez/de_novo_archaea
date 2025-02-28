import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANKS_LETTERS = ["k", "p", "c", "o", "f", "g", "s"]

# Read input file
INPUT = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/diamond_output_re_rank_1.tsv"
df = pd.read_csv(INPUT, sep="\t", header=None)
df.columns = ["query", "subject", "taxids", "evalue", "qcov"]

# Get the chosen rank
rank = "superkingdom"
RANK_LETTER = RANKS_LETTERS[RANKS.index(rank)]
OUTPUT = f"/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/{rank}_level.txt"

# For each unique query
taxids_count = {}
for query in df["query"].unique():
    # get the taxids
    taxids_lst = []
    taxids = df[df["query"] == query]["taxids"].tolist()
    for taxid in taxids:
        taxids_lst.append(str(taxid).split(";")[0])
    # Count each taxid
    taxids_count[query] = {t: taxids_lst.count(t) for t in taxids_lst}

# get the list of unique taxids
unique_taxids = set([t for taxids in taxids_count.values() for t in taxids.keys()])

# Retrieve the chosen rank
def get_rank(taxid, rank_letter):
    cmd = f"echo {taxid} | taxonkit reformat -I 1 -f {{{rank_letter}}}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    result_rank = result.stdout.split("\t")[1].strip()
    return taxid, result_rank

# Parallelize the class retrieval
classes = {}
with ThreadPoolExecutor() as executor:
    results = executor.map(get_rank, unique_taxids, [RANK_LETTER] * len(unique_taxids))
    for taxid, result_class in results:
        classes[taxid] = result_class

# For each query, get the class count
class_count = {}
for query, taxids in taxids_count.items():
    class_count[query] = {}
    for taxid in taxids:
        c = classes[taxid]
        if c in class_count[query]:
            class_count[query][c] += taxids[taxid]
        else:
            class_count[query][c] = taxids[taxid]

# Write the result to a file
with open(OUTPUT, "w") as f:
    for query, classes in class_count.items():
        for c, count in classes.items():
            f.write(f"{query}\t{c}\t{count}\n")