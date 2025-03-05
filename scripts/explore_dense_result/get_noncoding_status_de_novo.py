import os
import re
import glob
import pandas as pd
from Bio import SeqIO



DENSE_DIR = "/home/eliott.tempez/Documents/archaea_data/dense/"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"








if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    # Initialise data dtf
    col = ["genome", "de_novo_name", "ancestor_genome", "intergenic", "f0", "f1", "f2", "eval"]
    data = pd.DataFrame(columns=col)



    ### Get the info on all de novo genes ###
    denovo_dict_gS = {}
    # Get the names of all de novo genes
    for genome in genomes:
        denovo_dict_gS[genome] = {}
        # Read dense result
        denovo_file = os.path.join(DENSE_DIR, genome, "denovogenes.tsv")
        denovo_df = pd.read_csv(denovo_file, sep="\t", header=0)
        if denovo_df.empty:
            continue
        for denovo_gene in denovo_df["gene"]:
            denovo_dict_gS[genome][denovo_gene] = {}
    
    # Get the last matches in synteny for all de novo genes
    for genome in denovo_dict_gS:
        trg_file = os.path.join(DENSE_DIR, genome, "TRG_match_matrix.tsv")
        trg_df = pd.read_csv(trg_file, sep="\t", header=0)
        for denovo in denovo_dict_gS[genome]:
            # Get the row corresponding to the de novo gene (we must have one)
            try:
                matches = trg_df[trg_df["CDS"] == f"{denovo}_elongated"]
                if matches.empty:
                    raise ValueError(f"No corresponding row found for {denovo}_elongated in {genome}")
            except ValueError as e:
                print(e)
                continue
            # Iterate over the row in reverse until we hit a cds
            i = len(matches.columns)
            cell = ""
            while cell != "CDS":
                i -= 1
                cell = matches.iloc[0, i]
                if "gS" in cell:
                    denovo_dict_gS[genome][denovo][matches.columns[i]] = []
                    
    # Get the loci for each noncoding match of each de novo gene
    for genome in denovo_dict_gS:
        for denovo in denovo_dict_gS[genome]:
            for ancester_sp in denovo_dict_gS[genome][denovo]:
                tblastn_result_file = f"{DENSE_DIR}{genome}/blast_out/TRG_multielongated_tblastn_{ancester_sp}_genome.out"
                # Get the tblastn result corresponding to the de novo gene
                tblastn_df = pd.read_csv(tblastn_result_file, sep="\t", header=None, comment="#")
                tblastn_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlength", "qcov", "sframe"]
                tblastn_df = tblastn_df[tblastn_df["qseqid"] == f"{denovo}_elongated"]
                # Keep only qcov > 50
                tblastn_df = tblastn_df[tblastn_df["qcov"] > 50]
                # Keep only the best evalue
                tblastn_df = tblastn_df.sort_values("evalue")
                contig = tblastn_df.iloc[0]["sseqid"]
                start_loci = tblastn_df.iloc[0]["sstart"]
                end_loci = tblastn_df.iloc[0]["send"]
                denovo_dict_gS[genome][denovo][ancester_sp] = [contig, start_loci, end_loci]

# Print the info on all de novo genes
print("genome\tdenovo\tancester_sp\tloci")
for genome in denovo_dict_gS:
    for denovo in denovo_dict_gS[genome]:
        for ancester_sp in denovo_dict_gS[genome][denovo]:
            print(f"{genome}\t{denovo}\t{ancester_sp}\t{denovo_dict_gS[genome][denovo][ancester_sp]}")
    print("\n")


altframes = {}
### Get the corresponding area in the ancestor genome ###
for genome in denovo_dict_gS:
    altframes[genome] = {}
    for denovo in denovo_dict_gS[genome]:
        altframes[genome][denovo] = {}
        for ancester_sp in denovo_dict_gS[genome][denovo]:
            altframes[genome][denovo][ancester_sp] = {}
            # Get all CDSs for the ancestor genome of the de novo gene
            ancester_sp_genes = []
            contig, start_match, end_match = denovo_dict_gS[genome][denovo][ancester_sp]
            gbk_pattern = os.path.join(DATA_DIR, "annotation", ancester_sp + ".gb*")
            gbk_files = [f for f in glob.glob(gbk_pattern) if not f.endswith(".fai")]
            gbk_file = gbk_files[0]
            for record in SeqIO.parse(gbk_file, "genbank"):
                if record.id == contig:
                    for feature in record.features:
                        if feature.type == "CDS":
                            gene_strand = feature.location.strand
                            gene_start = feature.location.start
                            gene_end = feature.location.end
                            ancester_sp_genes.append([gene_strand, (gene_start, gene_end)])
            # Get each cds in the ancestor genome for which our de novo gene is in
            de_novo_is_on_plus = start_match < end_match
            all_loci = list(range(start_loci, end_loci)).sort() if de_novo_is_on_plus else list(range(end_loci, start_loci)).sort()
            loci_in_gene = []
            altframes = {}
            for ancester_cds in ancester_sp_genes:
                if de_novo_is_on_plus and ancester_cds[0] == 1:
                    if any(gene_start <= i <= gene_end for i in all_loci):                        
                        # Get the loci that are in the cds
                        loci_in_gene += [i for i in all_loci if gene_start <= i <= gene_end]
                        # Get the frame
                        frame = (min(loci_in_gene) - gene_start) % 3
                        if frame not in altframes[genome][denovo][ancester_sp]:
                            altframes[genome][denovo][ancester_sp][frame] = len([i for i in all_loci if gene_start <= i <= gene_end])
                        else:
                            altframes[genome][denovo][ancester_sp][frame] += len([i for i in all_loci if gene_start <= i <= gene_end])
                elif not de_novo_is_on_plus and ancester_cds[0] == -1:
                    if any(gene_start <= i <= gene_end for i in all_loci):
                        # Get the loci that are in the cds
                        loci_in_gene += [i for i in all_loci if gene_start <= i <= gene_end]
                        # Get the frame
                        frame = (gene_start - max(loci_in_gene)) % 3
                        if frame not in altframes[genome][denovo][ancester_sp]:
                            altframes[genome][denovo][ancester_sp][frame] = len([i for i in all_loci if gene_start <= i <= gene_end])
                        else:
                            altframes[genome][denovo][ancester_sp][frame] += len([i for i in all_loci if gene_start <= i <= gene_end])
            # Get the loci that are in the intergenic
            loci_in_intergenic = [i for i in all_loci if i not in loci_in_gene]
            altframes[genome][denovo][ancester_sp]["intergenic"] = len(loci_in_intergenic)

# Print the results
for genome in altframes:
    for denovo in altframes[genome]:
        for ancester_sp in altframes[genome][denovo]:
            print(f"{genome}\t{denovo}\t{ancester_sp}\t{altframes[genome][denovo][ancester_sp]}")
        print("\n")
                    
                    
            