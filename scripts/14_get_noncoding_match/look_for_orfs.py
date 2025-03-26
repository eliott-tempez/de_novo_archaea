import pandas as pd
import os
from Bio import SeqIO
from orffinder import orffinder
from Bio.SeqRecord import SeqRecord


DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"



def get_sequence_from_loci(genome, contig, start, end):
    fa_file = os.path.join(DATA_DIR, "fasta_renamed", genome + ".fa")
    for record in SeqIO.parse(fa_file, "fasta"):
        if str(record.name) == str(contig):
            seq = record.seq[start:end]
            name = f"{record.name}_{start}_{end}"
            return SeqRecord(seq, id=name, description="")
    return None



# Read the data
data = pd.read_csv("/home/eliott.tempez/Documents/M2_Stage_I2BC/status_several_origins.tsv", sep="\t")
# get list [outgroup, noncoding_match_start, noncoding_match_end] for each row
nc_zones = []
for row in data.iterrows():
    row = row[1]
    outgroup = row["outgroup"]
    contig = row["noncoding_match_contig"]
    noncoding_match_start = row["noncoding_match_start"]
    noncoding_match_end = row["noncoding_match_end"]
    strand = row["noncoding_match_strand"]
    nc_zones.append((outgroup, contig, noncoding_match_start, noncoding_match_end, strand))
nc_zones = set(nc_zones)


# Get the sequence for each match
for outgroup, contig, start, end, strand in nc_zones:
    seq = get_sequence_from_loci(outgroup, contig, start, end)
    if strand == "-":
        seq = seq.reverse_complement()
    orfs = orffinder.getORFs(seq, minimum_length=30)
    # Keep only ORFs on the right strand
    orfs = [dict for dict in orfs if dict["sense"] == "+"]
    max_len = 0
    if orfs:
        for dict in orfs:
            if dict["length"] > max_len:
                max_len = dict["length"]
                max_start = dict["start"]
                max_end = dict["end"]
    print(f"De novo with match in {outgroup} ({contig}, {start}-{end}):")
    print(f"Match of length {len(seq)} has {len(orfs)} ORFs with max length {max_len} ({max_start}-{max_end}):")
    print("-"*100)
    print(f"{' '*round((max_start/len(seq))*100-1)}{'-'*round(((max_end-max_start)/len(seq))*100)}")
    print("\n\n")
    






