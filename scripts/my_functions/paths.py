"""
This script will import all hardcoded paths from the yaml file
"""
CLUSTER = False
HOME = False


import yaml
import os
import sys

# Get file directory
file_dir = os.path.dirname(os.path.abspath(__file__))
with open(file_dir + "/filepaths.yaml", "r") as file:
    config = yaml.safe_load(file)

# Define paths
origin = "local_paths" if not CLUSTER else "cluster_paths"
origin = "home_paths" if HOME else origin
GENOMES_LIST = config[origin]["genomes_list"]
GENERA_DIR = config[origin]["genera_dir"]
DENSE_DIR = config[origin]["dense_dir"]
GFF_DIR = config[origin]["gff_dir"]
FA_DIR = config[origin]["fa_dir"]
CDS_DIR = config[origin]["cds_dir"]
GBK_DIR = config[origin]["gbk_dir"]

