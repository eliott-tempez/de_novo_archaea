import pandas as pd
import subprocess

INPUT = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/diamond_output_re_rank_1.tsv"
OUTPUT = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/classes_level.txt"

df = pd.read_csv(INPUT, sep="\t", header=None)
df.columns = ["query", "subject", "taxids", "evalue", "qcov"]

# For each unique query
taxids_count = {}
for query in df["query"].unique():
    # get the taxids
    taxids_lst = []
    taxids = df[df["query"] == query]["taxids"].tolist()
    for taxid in taxids:
        taxids_lst += str(taxid).split(";")
    # Count each taxid
    taxids_count[query] = {t: taxids_lst.count(t) for t in taxids_lst}

# get the list of unique taxids
unique_taxids = set([t for taxids in taxids_count.values() for t in taxids.keys()])
# For each taxid, get the class
classes = {}
for taxid in unique_taxids:
    cmd = f"echo {taxid} | taxonkit reformat -I 1 -f {{c}}"
    # execute and get the output
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    result_class = result.stdout.split("\t")[1].strip()
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
        f.write(f"{query}\n")
        for c, count in classes.items():
            f.write(f"{c}\t{count}\n")
        f.write("\n")







