import os

FASTA_FOLDER = "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/data/complete_122/fasta_renamed/"
GFF_folder = "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/data/complete_122/gff3/"
CDS_folder = "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/data/dense/archive_tree_v1/GCA_000195935@Pyrococcus_abyssi_GE5/CDS/"
EXTENDED_TRG_folder = "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/data/complete_122/extended_trg/"
TMP_DIR = os.path.join(os.getcwd(), "tmp/")


def get_query_from_user():
    """Get the CDS and the database from the user"""
    print("\n#### BLAST CDS ####")
    cds_str = input("<Name of CDS> <0:unextended 1:extended>: ")
    cds = extract_cds(cds_str)
    db_str = input("\n<DB name> <0:CDS 1:whole genome>: ")
    db, blast_type = extract_db(db_str)
    return cds, db, blast_type


def extract_cds(cds_str):
    """Get the faa file corresponding to the query"""
    query, extended = cds_str.split()
    if extended == "0":
        cds_file = get_exact_file(query, CDS_folder, ".faa")
    elif extended == "1":
        cds_file = get_exact_file(query, EXTENDED_TRG_folder, ".faa")
    else:
        print("Error: extended should be 0 or 1")
        exit
    should_extract_cds = input("\nDo you want to extract a CDS from the genome? (y/n): ")
    if should_extract_cds == "y":
        cds = extract_cds_from_genome(cds_file)
    elif should_extract_cds == "n":
        cds = cds_file
    else:
        print("Error: please enter y or n")
        exit
    return cds


def extract_cds_from_genome(cds_file):
    """Extract a fasta file with only the wanted CDS"""
    # Create temp dir if it doesn't exist
    if not os.path.exists(TMP_DIR):
        os.makedirs(TMP_DIR)
    # Set working directory
    os.chdir(TMP_DIR)
    # Get the genome name
    genome_name = input("Enter the genome name: ")
    # Create temp file with the wanted CDS
    os.system(f"grep -A 1 {genome_name} {cds_file} > {TMP_DIR}cds.faa")
    return f"{TMP_DIR}cds.faa"


def extract_db(db_str):
    """Get the faa file corresponding to the query"""
    query, extended = db_str.split()
    if extended == "0":
        db = get_exact_file(query, CDS_folder, ".faa")
        blast_type = "prot"
    elif extended == "1":
        db = get_exact_file(query, FASTA_FOLDER, ".fa")
        blast_type = "nucl"
    else:
        print("Error: extended should be 0 or 1")
        exit
    return db, blast_type


def get_exact_file(query, folder, extension):
    """Get the exact file corresponding to the query"""
    # folder content in alphabetical order
    folder_content  = sorted(os.listdir(folder))
    if query in folder_content:
        return folder + query
    else:
        # look for all the files containing the query
        matching_files = [f for f in folder_content if query in f and f.endswith(extension)]
        if len(matching_files) == 0:
            print("Error: no file found for the query")
            exit
        elif len(matching_files) == 1:
            print(f"File found: {matching_files[0]}")
            return folder + matching_files[0]
        else:
            print("Multiple files found for the query; please choose one:")
            for i, f in enumerate(matching_files):
                print(f"{i}: {f}")
            choice = input("choice: ")
            return folder + matching_files[int(choice)]


def make_db(db, blast_type):
    """Make the blast database"""
    # Create temp dir if it doesn't exist
    if not os.path.exists(TMP_DIR):
        os.makedirs(TMP_DIR)
    # Set working directory
    os.chdir(TMP_DIR)
    # Make blast database
    os.system(f"makeblastdb -in {db} -dbtype {blast_type} -out {TMP_DIR}db > /dev/null")


def run_blast(cds, db, blast_type):
    blast_cmd = "blastp" if blast_type == "prot" else "tblastn"
    # Run blast
    os.system(f"{blast_cmd} -query {cds} -db {TMP_DIR}db > {TMP_DIR}blast_result.txt")
    os.system(f"less {TMP_DIR}blast_result.txt")



if __name__ == "__main__":
    # get files from user
    cds, db, blast_type = get_query_from_user()
    # make blast database
    make_db(db, blast_type)
    # run blast and print results
    print("\nRunning blast...")
    run_blast(cds, db, blast_type)
    # remove temp dir
    os.system(f"rm -r {TMP_DIR}")
    print("\n#### END ####")