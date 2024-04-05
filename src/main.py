import sys
import os
import pandas as pd

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from utils import UnitStandardiser, RawCleaner, Binarizer
from generate_datasets import OrgDatasets, ProteinDatasets, TypeDatasets, BioassayDatasets
from default import DATAPATH, ST_TYPES
from split_datasets import Splitter


pathogen = sys.argv[1]

#Standardise units
df = pd.read_csv(os.path.join(DATAPATH, pathogen, "{}_original.csv".format(pathogen)), low_memory=False)
rc = RawCleaner()
df = rc.run(df)

us = UnitStandardiser()
df = us.standardise(df)
df.to_csv(os.path.join(DATAPATH, pathogen, "{}_processed.csv".format(pathogen)), index=False)

bin = Binarizer()
df = bin.run(df)
df.to_csv(os.path.join(DATAPATH, pathogen, "{}_binary.csv".format(pathogen)), index=False)

td = TypeDatasets(pathogen)
lc, hc = td.run_any(df)
try:
    lc.to_csv(os.path.join(DATAPATH, pathogen, "{}_anytype_lc.csv".format(pathogen)), index=False)
except:
    print("No Low Cut data for {} assay types".format(pathogen))
try:
    hc.to_csv(os.path.join(DATAPATH, pathogen, "{}_anytype_hc.csv".format(pathogen)), index=False)
except:
    print("No High Cut data for {} assay types".format(pathogen))

od = OrgDatasets(pathogen)
lc, hc = od.run_any_org(df)
try:
    lc.to_csv(os.path.join(DATAPATH, pathogen, "{}_anyorg_lc.csv".format(pathogen)), index=False)
except:
    print("No Low Cut data for {} org assay types".format(pathogen))
try:
    hc.to_csv(os.path.join(DATAPATH, pathogen, "{}_anyorg_hc.csv".format(pathogen)), index=False)
except:
    print("No High Cut data for {} org assay types".format(pathogen))

org_data = od.run_top_org(df)
for k,v in org_data.items():
    try:
        v[0].to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_org_lc.csv".format(pathogen, k)), index=False)
    except:
        print("No Low Cut data for {} {} org assay".format(pathogen, k))
    try:           
        v[1].to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_org_hc.csv".format(pathogen, k)), index=False)
    except:
        print("No High Cut data for {} {} org assay".format(pathogen, k))

pd = ProteinDatasets(pathogen)
lc, hc = pd.run_any_prot(df)
try:
    lc.to_csv(os.path.join(DATAPATH, pathogen, "{}_anyprot_lc.csv".format(pathogen)), index=False)
except:
    print("No Low Cut data for {} prot assay types".format(pathogen))
try:
    hc.to_csv(os.path.join(DATAPATH, pathogen, "{}_anyprot_hc.csv".format(pathogen)), index=False)
except:
    print("No High Cut data for {} prot assay types".format(pathogen))

prot_data = pd.run_top_prot(df)
for k,v in prot_data.items():
    try:
        v[0].to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_prot_lc.csv".format(pathogen, k)), index=False)
    except:
        print("No Low Cut data for {} {} org assay".format(pathogen, k))
    try:           
        v[1].to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_prot_hc.csv".format(pathogen, k)), index=False)
    except:
        print("No High Cut data for {} {} org assay".format(pathogen, k))

bd = BioassayDatasets(pathogen)
for st_type in ST_TYPES:
    lc, hc = bd.run_merged_bioassays(df, st_type)
    try:
        lc.to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_lc.csv".format(pathogen, st_type)), index=False)
    except:
        print("No Low Cut data for {}, {} assay types".format(pathogen, st_type))
    try:
        hc.to_csv(os.path.join(DATAPATH, pathogen, "{}_{}_hc.csv".format(pathogen, st_type)), index=False)
    except:
        print("No High Cut data for {}, {}, assay types".format(pathogen, st_type))

"""
s = Splitter(pathogen)
s.create_summary()
s.create_directoy_structure()
s.create_input_files()
s.train_test_split()
"""