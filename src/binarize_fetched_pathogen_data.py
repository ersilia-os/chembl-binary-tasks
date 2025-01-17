import pandas as pd
from tqdm import tqdm
import os
import collections
import numpy as np
import argparse
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root))
from default import *

parser = argparse.ArgumentParser(description='Binarize fetched pathogen data')
parser.add_argument('--pathogen_code', type=str, help='Pathogen code')
parser.add_argument('--output_dir', type=str, help='Data directory')

args = parser.parse_args()
data_dir = args.output_dir
pathogen_code = args.pathogen_code

# Loading the data
df = pd.read_csv(os.path.join(data_dir, pathogen_code, "01_{0}_cleaned.csv".format(pathogen_code)))
print("Considering only organism target types")
print("Before: {0}".format(df.shape))
df = df[df["target_type"] == "ORGANISM"]
print("After: {0}".format(df.shape))
df.drop(columns=["target_type"], inplace=True)
print("Considering only functional assay types")
print("Before: {0}".format(df.shape))
print(df.value_counts("assay_type"))
df = df[df["assay_type"] == "F"]
print("After: {0}".format(df.shape))
df.drop(columns=["assay_type"], inplace=True)

tasks_dir = os.path.join(data_dir, pathogen_code, "02_raw_tasks")
if not os.path.exists(tasks_dir):
    os.makedirs(tasks_dir)

def pchembl_binarizer(df, prefix):
    df = df[df["pchembl_value"].notnull()]
    df = df[df["pchembl_relation"].notnull()]
    data = {}
    for pchembl_cutoff in PCHEMBL_CUTOFFS:
        da = df[df["pchembl_value"] >= pchembl_cutoff]
        da = da[da["pchembl_relation"] != "<"]
        if da.shape[0] < MIN_POSITIVES:
            print("Not enough positives for pchembl cutoff {0}, {1}".format(pchembl_cutoff, da.shape[0]))
            continue
        di = df[df["pchembl_value"] < pchembl_cutoff]
        di = di[di["pchembl_relation"] != ">"]
        actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
        inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
        data["{0}_pchembl_value_{1}".format(prefix, pchembl_cutoff)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                             "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                             "pchembl_value_{0}".format(pchembl_cutoff): [1] * len(actives) + [0] * len(inactives)})
    for percentile in PERCENTILES:
        print("Percentile: {0}".format(percentile))
        df = df[df["pchembl_relation"] != ">"]
        df.loc[df["pchembl_relation"] == "<", "pchembl_value"] = 0
        N = df.shape[0]
        n = int(N * percentile / 100)
        if n < MIN_POSITIVES:
            continue
        df = df.sort_values("pchembl_value", ascending=False)
        da = df.head(n)
        di = df.tail(N - n)
        actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
        inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
        data["{0}_pchembl_percentile_{1}".format(prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                             "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                             "pchembl_percentile_{}".format(percentile): [1] * len(actives) + [0] * len(inactives)})
    print("Collected {0} datasets".format(len(data)))
    return data


def percentage_activity_binarizer(df, prefix):
    df = df[df["standard_value"].notnull()]
    df = df[df["standard_relation"].notnull()]
    df = df[df["standard_units"] == "%"]
    df = df[df["direction_flag"] == 1]
    data = {}
    for percentage_activity_cutoff in PERCENTAGE_ACTIVITY_CUTOFFS:
        da = df[df["standard_value"] >= percentage_activity_cutoff]
        da = da[da["standard_relation"] != "<"]
        if da.shape[0] < MIN_POSITIVES:
            continue
        di = df[df["standard_value"] < percentage_activity_cutoff]
        di = di[di["standard_relation"] != ">"]
        actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
        inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
        data["{0}_percentage_activity_{1}".format(prefix, percentage_activity_cutoff)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                             "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                             "percentage_activity_{0}".format(percentage_activity_cutoff): [1] * len(actives) + [0] * len(inactives)})
    for percentile in PERCENTILES:
        df = df[df["standard_relation"] != ">"]
        df.loc[df["standard_relation"] == "<", "standard_value"] = 0
        N = df.shape[0]
        n = int(N * percentile / 100)
        if n < MIN_POSITIVES:
            continue
        df = df.sort_values("standard_value", ascending=False)
        da = df.head(n)
        di = df.tail(N - n)
        actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
        inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
        data["{0}_percentage_activity_percentile_{1}".format(prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                             "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                             "percentage_activity_percentile_{0}".format(percentile): [1] * len(actives) + [0] * len(inactives)})
    print("Collected {0} datasets".format(len(data)))
    return data


def active_inactive_binarizer(df, prefix):
    df = df[df["activity_flag"] != 0]
    data = {}
    da = df[df["activity_flag"] == 1]
    di = df[df["activity_flag"] == -1]
    actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
    inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
    data["{0}_labeled_active".format(prefix)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                            "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                            "labeled_active": [1] * len(actives) + [0] * len(inactives)})
    print("Collected {0} datasets".format(len(data)))
    return data


def others_binarizer(df, prefix):
    def split_by_percentiles(df_, direction, inner_prefix):
        data = {}
        for percentile in PERCENTILES:
            if direction == 1:
                df_ = df_[df_["standard_relation"] != ">"]
                df_.loc[df_["standard_relation"] == "<", "standard_value"] = 0
                df_ = df_.sort_values("standard_value", ascending=False)
                N = df_.shape[0]
                n = int(N * percentile / 100)
                if n < MIN_POSITIVES:
                    continue
                da = df_.head(n)
                di = df_.tail(N - n)
                actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
                inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
                data["{0}_{1}_percentile_{2}".format(prefix, inner_prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                            "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                            "{0}_percentile_{1}".format(inner_prefix, percentile): [1] * len(actives) + [0] * len(inactives)})
            elif direction == -1:
                df = df_[df_["standard_relation"] != "<"]
                df.loc[df["standard_relation"] == ">", "standard_value"] = 10**15
                df = df.sort_values("standard_value", ascending=True)
                N = df.shape[0]
                n = int(N * percentile / 100)
                if n < MIN_POSITIVES:
                    continue
                da = df.head(n)
                di = df.tail(N - n)
                actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
                inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
                data["{0}_{1}_percentile_{2}".format(prefix, inner_prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                            "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                            "{0}_percentile_{1}".format(inner_prefix, percentile): [1] * len(actives) + [0] * len(inactives)})
            else:
                continue
            return data

    data = {}
    df = df[df["standard_value"].notnull()]
    df = df[df["standard_relation"].notnull()]
    df = df[df["standard_units"].notnull()]
    df = df[df["pchembl_value"].isnull()]
    df = df[df["pchembl_relation"].isnull()]
    df = df[df["standard_units"] != "%"]
    type_units = []
    for t,u in df[["standard_type", "standard_units"]].values:
        type_units += [(t,u)]
    type_units = list(set(type_units))
    for t, u in type_units:
        df_ = df[(df["standard_type"] == t) & (df["standard_units"] == u)]
        if df_.shape[0] < MIN_SIZE_ANY_TASK:
            continue
        direction = list(set(df_["direction_flag"].tolist()))
        if len(direction) != 1:
            continue
        direction = direction[0]
        direction_confidence = list(set(df_["direction_confidence"].tolist()))
        if len(direction_confidence) != 1:
            continue
        direction_confidence = direction_confidence[0]
        if direction_confidence == 0 or direction == 0:
            for k,v in split_by_percentiles(df_, 1, "{0}_{1}_uncertain".format(t.lower(), u.lower())):
                data[k] = v
            for k,v in split_by_percentiles(df_, -1, "{0}_{1}_uncertain".format(t.lower(), u.lower())):
                data[k] = v
        else:
            for k,v in split_by_percentiles(df_, direction, "{0}_{1}".format(t.lower(), u.lower())):
                data[k] = v
    print("Collected {0} datasets".format(len(data)))
    return data

def append_data(all_datasets, data):
    for k,v in data.items():
        if k not in all_datasets:
            all_datasets[k] = v
        else:
            raise Exception("Dataset {0} already exists".format(k))
    return all_datasets

# Select up to 5 assays with at least 100 molecules

def create_datasets_by_top_assays(df, all_datasets, priority):
    assay_ids = [x for x in df.value_counts("assay_id").index]
    counts = [x for x in df.value_counts("assay_id")]
    sel_assay_ids = []
    for aid, count in zip(assay_ids, counts):
        if count >= MIN_SIZE_ASSAY_TASK:
            sel_assay_ids.append(aid)
        else:
            break
    sel_assay_ids = sel_assay_ids[:MAX_NUM_INDEPENDENT_ASSAYS]
    for aid in sel_assay_ids:
        print("Assay ID: {0}".format(aid))
        dt = df[df["assay_id"] == aid]
        activity_types = [x for x in dt.value_counts("standard_type").index]
        counts = [x for x in dt.value_counts("standard_type")]
        sel_activity_types = []
        for at, count in zip(activity_types, counts):
            if count >= MIN_SIZE_ASSAY_SUBTASK:
                sel_activity_types.append(at)
        sel_activity_types = sel_activity_types[:MAX_NUM_ASSAY_SUBTASKS]
        for activity_type in activity_types:
            prefix = "{0}_assay_{1}_{2}".format(priority, aid, activity_type)
            dtt = dt[dt["standard_type"] == activity_type]
            for has_pchembl in [True, False]:
                if has_pchembl:
                    dttp = dtt[dtt["pchembl_value"].notnull()]
                    if dttp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                        continue
                    data = pchembl_binarizer(dttp, prefix=prefix)
                    all_datasets = append_data(all_datasets, data)
                else:
                    dttp = dtt[dtt["pchembl_value"].isnull()]
                    if dttp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                        continue
                    for has_percentage_activity in [True, False]:
                        if has_percentage_activity:
                            dttpp = dttp[dttp["standard_units"] == "%"]
                            if dttpp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                                continue
                            data = percentage_activity_binarizer(dttpp, prefix=prefix)
                            all_datasets = append_data(all_datasets, data)
                        else:
                            dttpp = dttp[dttp["standard_value"].isnull()]
                            if dttpp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                                continue
                            data = others_binarizer(dttpp, prefix=prefix)
                            all_datasets = append_data(all_datasets, data)
    return all_datasets

# Global active / inactive binarizer

def create_datasets_by_active_inactive(df, all_datasets, priority):
    prefix = "{0}_all".format(priority)
    data = active_inactive_binarizer(df, prefix=prefix)
    all_datasets = append_data(all_datasets, data)
    return all_datasets

def create_datasets_by_major_types(df, all_datasets, priority):
    counter = collections.defaultdict(int)
    for v in df[["target_id", "standard_type", "standard_units"]].values:
        counter[(v[0], v[1], v[2])] += 1
    selected_units = sorted(counter.items(), key=lambda x: x[1], reverse=True)[:MAX_NUM_INDEPENDENT_ASSAYS]
    for r in selected_units:
        r = r[0]
        target_id = r[0]
        standard_type = r[1]
        standard_units = r[2]
        dt = df[(df["target_id"] == target_id) & (df["standard_type"] == standard_type) & (df["standard_units"] == standard_units)]
        prefix = "{0}_target_{1}_{2}_{3}".format(priority, target_id, standard_type.lower(), standard_units.lower())
        for has_pchembl in [True, False]:
            if has_pchembl:
                dtp = dt[dt["pchembl_value"].notnull()]
                if dtp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                    continue
                data = pchembl_binarizer(dtp, prefix=prefix)
                all_datasets = append_data(all_datasets, data)
            else:
                dtp = dt[dt["pchembl_value"].isnull()]
                if dtp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                    continue
                for has_percentage_activity in [True, False]:
                    if has_percentage_activity:
                        dtpp = dtp[dtp["standard_units"] == "%"]
                        if dtpp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                            continue
                        data = percentage_activity_binarizer(dtpp, prefix=prefix)
                        all_datasets = append_data(all_datasets, data)
                    else:
                        dtpp = dtp[dtp["standard_value"].isnull()]
                        if dtpp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                            continue
                        data = others_binarizer(dtpp, prefix=prefix)
                        all_datasets = append_data(all_datasets, data)
    return all_datasets 

def create_datasets_by_all_pchembl(df, all_datasets, priority):
    print("Considering only pchembl values")
    dtp = df[df["pchembl_value"].notnull()]
    if dtp.shape[0] < MIN_SIZE_ANY_TASK:
        return all_datasets
    prefix = "{0}_all".format(priority)
    data = pchembl_binarizer(dtp, prefix=prefix)
    all_datasets = append_data(all_datasets, data)
    return all_datasets

def create_datasets_by_all_percentage(df, all_datasets, priority):
    print("Considering only percentage activity")
    dtp = df[df["standard_units"] == "%"]
    dtp = dtp[dtp["standard_value"].notnull()]
    dtp = dtp[dtp["standard_relation"].notnull()]
    dtp = dtp[dtp["direction_flag"] == 1]
    if dtp.shape[0] < MIN_SIZE_ANY_TASK:
        return all_datasets
    prefix = "{0}_all".format(priority)
    data = percentage_activity_binarizer(dtp, prefix=prefix)
    all_datasets = append_data(all_datasets, data)
    return all_datasets

def create_datasets_by_grouping_percentiles(df, all_datasets, priority):
    print("Grouping percentiles")
    all_units = []
    data_actives = collections.defaultdict(list)
    data_inactives = collections.defaultdict(list)
    for v in df[["target_id", "standard_type", "standard_units"]].values:
        all_units.append((v[0], v[1], v[2]))
    all_units = list(set(all_units))
    for r in tqdm(all_units):
        dp = df[(df["target_id"] == r[0]) & (df["standard_type"] == r[1]) & (df["standard_units"] == r[2])]
        for direction in [1, -1]:
            if direction == 1:
                dp = dp[dp["standard_relation"] != ">"]
                dp.loc[dp["standard_relation"] == "<", "standard_value"] = 0
                dp = dp.sort_values("standard_value", ascending=False)
            elif direction == -1:
                dp = dp[dp["standard_relation"] != "<"]
                dp.loc[dp["standard_relation"] == ">", "standard_value"] = 10**15
                dp = dp.sort_values("standard_value", ascending=True)
            else:
                continue
            N = dp.shape[0]
            for percentile in PERCENTILES:
                n = int(N * percentile / 100)
                if n == 0:
                    continue
                da = dp.head(n)
                di = dp.tail(N - n)
                actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
                inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
                data_actives[percentile] += actives
                data_inactives[percentile] += inactives
    prefix = "{0}_grouped_percentiles".format(priority)
    data = {}
    for percentile in PERCENTILES:
        data["{0}_{1}".format(prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in data_actives[percentile]] + [x[0] for x in data_inactives[percentile]],
                                         "smiles": [x[1] for x in data_actives[percentile]] + [x[1] for x in data_inactives[percentile]],
                                         "percentile_{0}".format(percentile): [1] * len(data_actives[percentile]) + [0] * len(data_inactives[percentile])})
    all_datasets = append_data(all_datasets, data)
    return all_datasets

all_datasets = {}
all_datasets = create_datasets_by_top_assays(df, all_datasets, priority=1)
all_datasets = create_datasets_by_major_types(df, all_datasets, priority=2)
all_datasets = create_datasets_by_all_pchembl(df, all_datasets, priority=3)
all_datasets = create_datasets_by_all_percentage(df, all_datasets, priority=4)
all_datasets = create_datasets_by_grouping_percentiles(df, all_datasets, priority=5)
all_datasets = create_datasets_by_active_inactive(df, all_datasets, priority=6)

def disambiguate_data(df):
    ik2smi = {}
    ik2act = collections.defaultdict(list)
    columns = list(df.columns)
    assert len(columns) == 3, "Expected 3 columns"
    for k,v in df[[columns[0], columns[1]]].values:
        ik2smi[k] = v
    for k,v in df[[columns[0], columns[2]]].values:
        ik2act[k] += [v]
    ik2act = {k: int(np.max(v)) for k,v in ik2act.items()}
    R = []
    for k,v in ik2act.items():
        R += [[k, ik2smi[k], v]]
    return pd.DataFrame(R, columns=columns)

all_datasets = {k: disambiguate_data(v) for k,v in all_datasets.items()}

all_datasets_filtered = {}
for k,v in all_datasets.items():
    if v.shape[0] < MIN_SIZE_ANY_TASK:
        continue
    columns = list(v.columns)
    assert len(columns) == 3, "Expected 3 columns"
    n = v[columns[2]].sum()
    if n < MIN_POSITIVES:
        continue
    file_name = os.path.join(tasks_dir, "{0}.csv".format(k))
    print("Saving data in {0}".format(file_name))
    v.to_csv(file_name, index=False)