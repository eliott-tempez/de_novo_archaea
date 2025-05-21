def read_TRG_table(file_path):
    all_TRGs = []
    TRGs = {}
    with open(file_path, "r") as f_in:
        f_in.readline()
        for line in f_in:
            line_split = line.split()
            TRG = line_split[1]
            yeast_origin = line_split[2]
            gene_origin = line_split[3]
            time_div = line_split[5]
            all_TRGs.append(TRG)
            TRGs[TRG] = {"yeast_origin": yeast_origin, "gene_origin": gene_origin, "time_div": time_div}
    return TRGs, all_TRGs

initial_TRGs, all_initial_TRGs = read_TRG_table("/home/eliott.tempez/Documents/M2_Stage_I2BC/data/initial_TRG_table.tsv")
mikatae_TRGs, all_mikatae_TRGs = read_TRG_table("/home/eliott.tempez/Documents/M2_Stage_I2BC/data/mikatae_TRG_table.tsv")


# Check the number of diverging TRGs
# TRGs in common
common_TRGs = set(initial_TRGs.keys()).intersection(set(mikatae_TRGs.keys()))
print(f"Number of TRGs in common: {len(common_TRGs)}")
# TRGs in initial but not in mikatae
initial_only_TRGs = set(initial_TRGs.keys()).difference(set(mikatae_TRGs.keys()))
print(f"Number of TRGs in initial but not in mikatae: {len(initial_only_TRGs)}")
# TRGs in mikatae but not in initial
mikatae_only_TRGs = set(mikatae_TRGs.keys()).difference(set(initial_TRGs.keys()))
print(f"Number of TRGs in mikatae but not in initial: {len(mikatae_only_TRGs)}")


