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
PCHEMBL_CUTOFFS = [5, 6, 7, 8, 9]
PERCENTAGE_ACTIVITY_CUTOFFS = [50, 75, 90]
PERCENTILES = [1, 5, 10, 25, 50]
MIN_SIZE_ASSAY_TASK = 500
MIN_SIZE_ASSAY_SUBTASK = 99
MIN_SIZE_ANY_TASK = 99
MAX_NUM_INDEPENDENT_ASSAYS = 5
MAX_NUM_ASSAY_SUBTASKS = 3
MIN_POSITIVES = 10
MAX_NUM_INDEPENDENT_UNITS = 5
DATASET_SIZE_LIMIT = 1e6