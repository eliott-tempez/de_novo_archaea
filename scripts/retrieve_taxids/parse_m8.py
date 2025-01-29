import os

M8_FOLDER = "/datas/ELIOTT/archaea_data/retrieve_taxid/diamond/"

taxid_dict = {}
m8_files = os.listdir(M8_FOLDER)
for file in m8_files:
    organism = ".".join(file.split(".")[:-1])
    taxid_dict[organism] = {}
    with open(M8_FOLDER + file, "r") as f_in:
        for line in f_in:
            taxids = [int(t) for t in (line.split("\t")[2]).split(";")]
            for taxid in taxids:
                if taxid in taxid_dict[organism]:
                    taxid_dict[organism][taxid] += 1
                else:
                    taxid_dict[organism][taxid] = 1

    # sort by number of occurences
    taxid_dict[organism] = dict(sorted(taxid_dict[organism].items(), key=lambda item: item[1], reverse=True))

# print results by alphabetical order
for organism in sorted(taxid_dict.keys()):
    print(f"{organism}: {taxid_dict[organism]}\n")