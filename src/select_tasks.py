import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import collections
import argparse

parser = argparse.ArgumentParser("Select 25 tasks")
parser.add_argument("--pathogen_code", type=str)
parser.add_argument("--output_dir", type=str)
args = parser.parse_args()

data_dir = args.output_dir
pathogen_code = args.pathogen_code
tasks_dir = os.path.join(data_dir, pathogen_code, "02_raw_tasks")

dm = pd.read_csv(os.path.join(data_dir, pathogen_code, "04_modelability.csv"))
ratios = []
priorities = []
for r in dm.iterrows():
    r = r[1]
    ratios += [r["num_pos_samples"] / r["num_samples"]]
    priorities += [int(r["task"][0])]
dm["ratio"] = ratios
dm["priority"] = priorities
modelable_tasks = set(dm[(dm["auroc_avg"] >= 0.7) & (dm["ratio"] <= 0.5)]["task"])

positive_sets = {}
for task in os.listdir(tasks_dir):
    fname = task[:-4]
    if fname not in modelable_tasks:
        continue
    df = pd.read_csv(os.path.join(tasks_dir, task))
    columns = list(df.columns)
    c = columns[-1]
    inchikeys = df[df[c] == 1]["inchikey"].tolist()
    positive_sets[fname] = set(inchikeys)


def positive_overlaps(positive_sets):
    tasks = sorted(list(positive_sets.keys()))
    R = []
    for task1 in tqdm(tasks):
        for task2 in tasks:
            if task1 >= task2:
                continue
            o = len(positive_sets[task1].intersection(positive_sets[task2]))
            n1 = len(positive_sets[task1])
            n2 = len(positive_sets[task2])
            oi = o / min(n1, n2)
            t = o / (len(positive_sets[task1].union(positive_sets[task2])))
            p1 = int(task1[0])
            p2 = int(task2[0])
            dm_ = dm[dm["task"] == task1]
            for v in dm_.values:
                auroc1 = v[1]
                n_total1 = v[3]
                break
            dm_ = dm[dm["task"] == task2]
            for v in dm_.values:
                auroc2 = v[1]
                n_total2 = v[3]
                break
            r = [task1, task2, n1, n2, o, oi, t, p1, p2, auroc1, auroc2, n_total1, n_total2]
            R += [r]
    return pd.DataFrame(R, columns=["task1", "task2", "n1", "n2", "overlap", "overlap_index", "jaccard_index", "priority1", "priority2", "auroc1", "auroc2", "n_total1", "n_total2"])

dp = positive_overlaps(positive_sets)
dp = dp.sort_values("jaccard_index", ascending=False)

to_remove = set()
for r in dp[dp["jaccard_index"] > 0.8].iterrows():
    r = r[1]
    if r["priority1"] < r["priority2"]:
        to_remove.add(r["task2"])
    elif r["priority1"] > r["priority2"]:
        to_remove.add(r["task1"])
    else:
        if r["n_total1"] > r["n_total2"]*1.25:
            to_remove.add(r["task2"])
        elif r["n_total2"] > r["n_total1"]*1.25:
            to_remove.add(r["task1"])
        else:
            if r["auroc1"] < r["auroc2"]:
                to_remove.add(r["task2"])
            else:
                to_remove.add(r["task1"])

to_remove = list(to_remove)

valid_tasks = modelable_tasks - set(to_remove)
valid_tasks = sorted(list(valid_tasks))

lb = np.percentile(dm["num_pos_samples"], 10)
ub = np.percentile(dm["num_pos_samples"], 90)

dm = dm[dm["task"].isin(valid_tasks)]

dm = dm[dm["num_pos_samples"] >= lb]
dm = dm[dm["num_pos_samples"] <= ub]

to_remove = []

percentile_names = collections.defaultdict(list)
for fname in dm["task"].tolist():
    if "_percentile_" in fname:
        agg_name = fname.split("_percentile_")[0]
        value = int(fname.split("_percentile_")[1])
        percentile_names[agg_name] += [(fname, value)]

for k,v in percentile_names.items():
    if len(v) > 1:
        for x in v:
            if x[1] == 50:
                to_remove += [x[0]]

to_remove = set(to_remove)
dm = dm[~dm["task"].isin(to_remove)]

to_remove = []

percentage_activity_names = collections.defaultdict(list)
for fname in dm["task"].tolist():
    if "_percentile_" in fname:
        continue
    if "_percentage_activity_" in fname:
        agg_name = fname.split("_percentage_activity_")[0]
        value = int(fname.split("_percentage_activity_")[1])
        percentage_activity_names[agg_name] += [(fname, value)]

for k,v in percentage_activity_names.items():
    if len(v) > 1:
        for x in v:
            if x[1] == 50:
                to_remove += [x[0]]

to_remove = set(to_remove)
dm = dm[~dm["task"].isin(to_remove)]

dm = dm.sort_values(by = ["priority", "auroc_avg"], ascending=[True, False]).head(25)
dm.to_csv(os.path.join(data_dir, pathogen_code, "05_selected_tasks.csv"), index=False)

print("Removing fingerprints and other unnecessary files")
os.remove(os.path.join(data_dir, pathogen_code, "03_fingerprints.npy"))
os.remove(os.path.join(data_dir, pathogen_code, "03_fingerprints_inchikeys.txt"))