import os
import shutil
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Wrap up tasks and clean output folder')
parser.add_argument('--pathogen_code', type=str, help='Pathogen code')
parser.add_argument('--output_dir', type=str, help='Output folder')
parser.add_argument('--flush', action='store_true', help='Flush output folder')
args = parser.parse_args()

pathogen_code = args.pathogen_code
data_dir = args.output_dir

ds = pd.read_csv(os.path.join(data_dir, pathogen_code, "05_selected_tasks.csv"))

if os.path.exists(os.path.join(data_dir, pathogen_code, "tasks")):
    shutil.rmtree(os.path.join(data_dir, pathogen_code, "tasks"))
os.makedirs(os.path.join(data_dir, pathogen_code, "tasks"))

for task in ds["task"].tolist():
    shutil.copy(os.path.join(data_dir, pathogen_code, "02_raw_tasks", task+".csv"), os.path.join(data_dir, pathogen_code, "tasks", task+".csv"))

ds.to_csv(os.path.join(data_dir, pathogen_code, "tasks_summary.csv"), index=False)

if args.flush:
    shutil.rmtree(os.path.join(data_dir, pathogen_code, "02_raw_tasks"))
    os.remove(os.path.join(data_dir, pathogen_code, "05_selected_tasks.csv"))
    os.remove(os.path.join(data_dir, pathogen_code, "01_abaumannii_cleaned.csv"))
    os.remove(os.path.join(data_dir, pathogen_code, "04_modelability.csv"))
    # os.remove(os.path.join(data_dir, pathogen_code, "00_abaumannii_original.csv"))