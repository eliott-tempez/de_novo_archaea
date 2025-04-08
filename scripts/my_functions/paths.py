"""
This script will import all hardcoded paths from the yaml file
"""
CLUSTER = False

import yaml
import os
import sys

# Get file directory
file_dir = os.path.dirname(os.path.abspath(__file__))
with open(file_dir + "/filepaths.yaml", "r") as file:
    config = yaml.safe_load(file)

# Define paths
local_or_cluster = "local_paths" if not CLUSTER else "cluster_paths"
GENOMES_LIST = config[local_or_cluster]["genomes_list"]
GENERA_DIR = config[local_or_cluster]["genera_dir"]
DENSE_DIR = config[local_or_cluster]["dense_dir"]
GFF_DIR = config[local_or_cluster]["gff_dir"]
FA_DIR = config[local_or_cluster]["fa_dir"]
CDS_DIR = config[local_or_cluster]["cds_dir"]
GBK_DIR = config[local_or_cluster]["gbk_dir"]