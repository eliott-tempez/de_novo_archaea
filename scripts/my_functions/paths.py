"""
This script will import all hardcoded paths from the yaml file
"""

import yaml

with open("/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/my_functions/filepaths.yaml", "r") as file:
    config = yaml.safe_load(file)

# Define paths
GENOMES_LIST = config["local_paths"]["genomes_list"]
GENERA_DIR = config["local_paths"]["genera_dir"]
DENSE_DIR = config["local_paths"]["dense_dir"]
GFF_DIR = config["local_paths"]["gff_dir"]
FA_DIR = config["local_paths"]["fa_dir"]
CDS_DIR = config["local_paths"]["cds_dir"]
GBK_DIR = config["local_paths"]["gbk_dir"]