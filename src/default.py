import os

abspath = os.path.dirname(os.path.abspath(__file__))

# Database defaults
DATABASE_NAME = "chembl_35"
CHEMBL_USR = "chembl_user"
CHEMBL_PWD = "aaa"

# Path defaults
CONFIGPATH = os.path.join(abspath, "..", "config")
TMPDIR = os.path.join(abspath, "..", "tmp")
PATHOGENSPATH = os.path.join(CONFIGPATH, "pathogens.csv")

# Parameters for binarization
MIN_SIZE_ASSAY_TASK = 1000 # Top assays with at least this data size will get a specific task
MIN_SIZE_PROTEIN_TASK = 250 # Top proteins with at least this data size will get a specific task
MIN_COUNT_POSITIVE_CASES = 30
TOP_ASSAYS = 3
TOP_PROTEINS = 3
DATASET_SIZE_LIMIT = 1e6
ST_TYPES = ["MIC", 'IZ', "IC50", "Inhibition", "Activity"]
SPLIT_METHOD = 'random'
