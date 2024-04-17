# chembl-binary-tasks
A repository to generate bioassay datasets from ChEMBL ready for downstream AI/ML modelling. This repository is based on the [antimicrobial-ml-tasks](https://github.com/ersilia-os/antimicrobial-ml-tasks) and [chembl-ml-tools](https://github.com/chembl-ml-tools), now Archived.

## Installation

To install the package in a conda environment, please run:
```
conda create -n chemblml python=3.10
conda activate chemblml
pip install git+https://github.com/ersilia-os/chembl-binary-tasks.git
```

## Requirements

These tools require access to a postgres database server containing the ChEMBL database. You may install ChEMBL in your own computer 
by following these instructions: [How to install ChEMBL](docs/install_chembl.md)

You can use the following code to check that the package is working. This test assumes that there is a DB user called `chembl_user` with permissions to read the database.

Before running, make sure that the postgres service with the ChEMBL database is up.

```
import pandas as pd
from chemblmltools import chembl_activity_target

df1 = chembl_activity_target(
        db_user='chembl_user',
        db_password='aaa',
        organism_contains='enterobacter',
        max_heavy_atoms=100)

print(df1.head(5))
```

# Create datasets

**1.** Make sure that the PostgreSQL server containing the ChEMBL database is running. In case of doubt, review the requirements section.

By default, the programs assume that PostgreSQL is running in the local computer, and that the database user `chembl_user` with
password `aaa` has read access to the tables of ChEMBL. This can be changed in program `scripts/defaults.py`.

**2.** Revise `src/default.py`. This file contains several default settings, including the path where data will be stored, or the minimum number of assays to consider. Modify according to your needs.

**3.** Edit the file `config/pathogens.csv` to select the pathogens for which we need models.

This file has two columns:

- **pathogen_code**: Choose a short code to identify the pathogen, alphanumeric only, **without spaces**. Example: "efaecium".
- **search_text**: A search string, *case insensitive*, to search for the pathogen name in the `organism` field 
in the ChEMBL database. Example: "Enterococcus Faecium".

**4.** Run the script `pathogens.py`
```
cd src
python pathogens.py
```

This will create folders under the specified `DATAPATH` folder with all the available data for the selected pathogens. In this example, we will simply refer to it as `/data` folder by convention.

**5.** Create the datasets for each individual pathogen. It includes a binary classication of all assays pulled together as well as independent assays. Run the script `main.py`, passing the pathogen code as argument. 
```
cd src
python main.py <pathogen_code>
```
Master files:
There are three master files in the /pathogen_name folder:
- pathogen_original.csv: the original file pulled from ChEMBL
- pathogen_processed.csv: the original file including the columns final_unit, transformer, final_value. The final_value column will contain the end result in the standardised unit as defined in `config/ucum.csv`
- pathogen_binary.csv: the processed file including the cut-offs for each assay (and unit type). Assay - unit combinations not selected in `config/st_type_summary_manual.csv` are not included here. Two columns are created: activity_lc and activity_hc, corresponding to the binary activity for the row if using the Low cut-off or the High cut-off. If there is a comment from the author (Active, Non Active) this determines the activity both in LC and HC. This file will be processed to create the following:

Outputs:
- Two file containing all the molecules and its binary classification (regardless of assay or target (whole cell organism or protein)), both using a high confidence threshold (pathogen_all_hc) and low confidence threshold (pathogen_all_lc) for the binarization of activities. At this stage, if molecules are duplicated, the values will be averaged and if they are > 0.5, the molecule will be considered active (1), else inactive (0)
- Two files containing all the molecules and its binary classification for whole cell assays, both using a high confidence threshold (pathogen_org_all_hc) and low confidence threshold (pathogen_org_all_lc) for the binarization of activities
- Two files containing all the molecules and its binary classification for protein assays, both using a high confidence threshold (pathogen_prot_all_hc) and low confidence threshold (pathogen_prot_all_lc) for the binarization of activities
- Files containing the top assays as determined by the thresholds in default.py (for example, an IC50 assay with over 250 molecules on it). Those are identified by pathogen_org_hc_top_{}, pathogen_org_lc_top_{}, pathogen_prot_hc_top_{}, pathogen_prot_hc_top_{}. The assay id and target protein can be found in the summary file.
- Files containing all results for selected assays (specified in ST_TYPEs in default.py). Currently those include MIC, IC50, IZ, Activity, Inhibition. They all relate to whole cell assays. The files produced are named pathogen_sttype_hc.csv and pathogen_sttype_lc.csv

A `pathogen_summary.csv` file is created containing a summary of the processing for each pathogen, and by running the following code snippet a full summary for all pathogens is created:

```
cd src
python summary.py
```

### Parameters
The following parameters are specified in `default.py`and can be modified according to the user needs:
MIN_SIZE_ASSAY_TASK = 1000 #Top assays with at least this data size will get a specific task
MIN_SIZE_PROTEIN_TASK = 250 #Top proteins with at least this data size will get a specific task
MIN_COUNT_POSITIVE_CASES = 30 #Minimum number of positive hits per assay
TOP_ASSAYS = 3 #Max number of selected organism assays
TOP_PROTEINS = 3 #Max number of selected protein assays
DATASET_SIZE_LIMIT = 1e6 # Limit the largest dataset
ST_TYPES = ["MIC", 'IZ', "IC50", "Inhibition", "Activity"] #Organism Bioassays that are merged together
SPLIT_METHOD = 'random' #split mode for ZairaChem (see below)

Paths to several files, including the `/data` folder, can be specified as well

# Train models with ZairaChem
Optionally, we offer the possibility of preparing the data directly for model training using ZairaChem. The datasets will automatically be split into train/test (80/20) split to perform model validation. To learn more about ZairaChem and install it, please see its [own repository](https://github.com/ersilia-os/zaira-chem).

Please run:
```
conda activate zairachem
cd src
python split_datasets.py <pathogen>
```
This will create the following files and folders in the `/data/<pathogen>` folder:
- input: will contain, once the split is performed:
  - input.csv: full input data
  - train.csv: input data for training
  - test.csv: input data for test
  - input_rejected.csv: cases that ZairaChem has rejected (typically because the molecule's SMILES is not valid)
  
- model: Contains the model definition, in the format used by ZairaChem
- test: Predictions for the test data and assessment reports of the model
- log: The log files resulting from the split, test and predict runs of ZairaChem

Once you are ready to train the models, you can run directly the files `"split_<pathogen>.sh` and `fit_predict_<pathogen>.sh` (be aware of the memory and time necessary to build ZairaChem models) or copy line by line the commands from the .sh files to fit and predict one model at a time:
```
conda activate zairachem
cd src
bash split_<pathogen>.sh
bash fit_predict_<pathogen>.sh
```


