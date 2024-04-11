import sys
import os
import pandas as pd

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from utils import UnitStandardiser, RawCleaner, Binarizer
from generate_datasets import OrgDatasets, ProteinDatasets, AllDatasets, BioassayDatasets
from default import DATAPATH, ST_TYPES

pathogen = sys.argv[1]

#Standardise units
df = pd.read_csv(os.path.join(DATAPATH, pathogen, "{}_original.csv".format(pathogen)), low_memory=False)
initial_len = len(df)
rc = RawCleaner()
df, no_activity_info, smi_not_processed = rc.run(df)

us = UnitStandardiser()
df = us.run(df)
df.to_csv(os.path.join(DATAPATH, pathogen, "{}_processed.csv".format(pathogen)), index=False)


bin = Binarizer()
df, assay_entries_discarded, total_inconsistent = bin.run(df)
df.to_csv(os.path.join(DATAPATH, pathogen, "{}_binary.csv".format(pathogen)), index=False)

unique_smis = df["compound_chembl_id"].nunique()
pos_lc = len(df[df["activity_lc"]==1])
pos_hc = len(df[df["activity_hc"]==1])

summary_dict = {"Dataset Original total molecules":initial_len,
                "Dataset Original: rows without sufficient activity info": no_activity_info,
                "Dataset Original: rows with unprocessable SMILES": smi_not_processed,
                "Dataset Original: rows with activity not in config file": assay_entries_discarded,
                "Dataset Original: rows with inconsistent activity info": total_inconsistent,
                "Dataset Original Final total molecules": initial_len - (no_activity_info+smi_not_processed+assay_entries_discarded+total_inconsistent),
                "Dataset Original Unique molecules:": unique_smis,
                "Dataset Original, low cutoff positive molecules:": pos_lc,
                "Dataset Original, high cutoff positive molecules:": pos_hc
                }

td = AllDatasets(pathogen)
print("All Data")
lc, hc, info_datasets = td.run_any(df)
try:
    lc.to_csv(os.path.join(DATAPATH, pathogen, "{}_all_lc.csv".format(pathogen)), index=False)
except:
    print("No Low Cut data for {} assay types".format(pathogen))
try:
    hc.to_csv(os.path.join(DATAPATH, pathogen, "{}_all_hc.csv".format(pathogen)), index=False)
except:
    print("No High Cut data for {} assay types".format(pathogen))

summary_dict.update(info_datasets)

od = OrgDatasets(pathogen)
print("All ORG Data")
lc, hc, info_datasets = od.run_any_org(df)
try:
    lc.to_csv(os.path.join(DATAPATH, pathogen, "{}_org_all_lc.csv".format(pathogen)), index=False)
except:
    print("No Low Cut data for {} org assay types".format(pathogen))
try:
    hc.to_csv(os.path.join(DATAPATH, pathogen, "{}_org_all_hc.csv".format(pathogen)), index=False)
except:
    print("No High Cut data for {} org assay types".format(pathogen))
summary_dict.update(info_datasets)

org_data = od.run_top_org(df)
print("Top ORG Data")
if len(org_data) == 0:
    info_datasets["Top ORG Assays with minimum data"] = 0
    summary_dict.update(info_datasets)
for k,v in org_data.items():
    print(k)
    info_datasets = v[2]
    summary_dict.update(info_datasets)
    try:
        v[0].to_csv(os.path.join(DATAPATH, pathogen, "{}_org_lc_{}.csv".format(pathogen, k)), index=False)
    except:
        print("No Low Cut data for {} {} org assay".format(pathogen, k))
    try:           
        v[1].to_csv(os.path.join(DATAPATH, pathogen, "{}_org_hc_{}.csv".format(pathogen, k)), index=False)
    except:
        print("No High Cut data for {} {} org assay".format(pathogen, k))

protd = ProteinDatasets(pathogen)
print("All PROT Data")
lc, hc, info_datasets = protd.run_any_prot(df)
try:
    lc.to_csv(os.path.join(DATAPATH, pathogen, "{}_prot_all_lc.csv".format(pathogen)), index=False)
except:
    print("No Low Cut data for {} prot assay types".format(pathogen))
try:
    hc.to_csv(os.path.join(DATAPATH, pathogen, "{}_prot_all_hc.csv".format(pathogen)), index=False)
except:
    print("No High Cut data for {} prot assay types".format(pathogen))
summary_dict.update(info_datasets)

prot_data = protd.run_top_prot(df)
print("top PROT Data")
if len(prot_data) == 0:
    info_datasets["Top Protein Assays with minimum data"] = 0
    summary_dict.update(info_datasets)
for k,v in prot_data.items():
    print(k)
    info_datasets = v[2]
    summary_dict.update(info_datasets)
    try:
        v[0].to_csv(os.path.join(DATAPATH, pathogen, "{}_prot_lc_{}.csv".format(pathogen, k)), index=False)
    except:
        print("No Low Cut data for {} {} org assay".format(pathogen, k))
    try:           
        v[1].to_csv(os.path.join(DATAPATH, pathogen, "{}_prot_hc_{}.csv".format(pathogen, k)), index=False)
    except:
        print("No High Cut data for {} {} org assay".format(pathogen, k))

bd = BioassayDatasets(pathogen)
print("Bioassays")
for st_type in ST_TYPES:
    print(st_type)
    lc, hc, info_datasets = bd.run_merged_bioassays(df, st_type)
    summary_dict.update(info_datasets)
    try:
        lc.to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_lc.csv".format(pathogen, st_type)), index=False)
    except:
        print("No Low Cut data for {}, {} assay types".format(pathogen, st_type))
    try:
        hc.to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_hc.csv".format(pathogen, st_type)), index=False)
    except:
        print("No High Cut data for {}, {}, assay types".format(pathogen, st_type))

summary = pd.DataFrame(summary_dict.items(), columns=["Concept", "Cases"])
summary.to_csv(os.path.join(DATAPATH, pathogen, f"{pathogen}_summary.csv"), index=False)