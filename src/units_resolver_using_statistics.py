import os
import sys
import collections
from tqdm import tqdm
import pandas as pd
import joblib
import numpy as np
from scipy.stats import mannwhitneyu

root = os.path.dirname(os.path.abspath(__file__))

sys.path.append(root)
from default import CONFIGPATH, TMPDIR


def annotate_activities_with_activity_flag():
    print("Loading standard text table...")
    df = pd.read_csv(os.path.join(CONFIGPATH, "standard_text.csv"))
    std_text = {}
    for v in df.values:
        if v[1] == 0:
            continue
        std_text[v[0]] = v[1]
    print("Loading activity comments table...")
    df = pd.read_csv(os.path.join(CONFIGPATH, "activity_comments.csv"), low_memory=False)
    act_comm = {}
    for v in df[["activity_comment", "activity_classified"]].values:
        if not v[0]:
            continue
        if str(v[0]).lower() == "nan" or str(v[0]).lower() == "na":
            continue
        if str(v[1]) == "0":
            continue
        if str(v[1]).lower() == "nan" or str(v[1]).lower() == "na" or str(v[1]).lower() == "to resolve":
            continue
        act_comm[v[0]] = int(v[1])
    print("Loading activities table...")
    df = pd.read_csv(os.path.join(TMPDIR, "activities.csv"), low_memory=False)
    print(df.columns)
    flags = []
    for v in tqdm(df[["activity_id", "activity_comment", "standard_text_value"]].values):
        if v[2] in std_text:
            flags += [std_text[v[2]]]
        elif v[1] in act_comm:
            flags += [act_comm[v[1]]]
        else:
            flags += [0]
    print("Adding activity flag to activities table...")
    df["activity_flag"] = flags
    df.to_csv(os.path.join(TMPDIR, "activities.csv"), index=False)


def create_annotated_actives_and_inactives_dictionaries():
    print("Loading activities table...")
    df = pd.read_csv(os.path.join(TMPDIR, "activities.csv"), low_memory=False)
    actives = collections.defaultdict(list)
    inactives = collections.defaultdict(list)
    for v in tqdm(df[["assay_id", "standard_type", "standard_units", "standard_value", "activity_flag"]].values):
        k = (v[0], v[1], v[2])
        if v[4] == 1:
            actives[k] += [v[3]]
        elif v[4] == -1:
            inactives[k] += [v[3]]
        else:
            continue
    return actives, inactives


def assign_directions_with_statistical_test(actives, inactives):
    units_directions = collections.defaultdict(list)
    all_keys = set(actives.keys()).union(set(inactives.keys()))
    for key in tqdm(all_keys):
        if key not in actives:
            continue
        if key not in inactives:
            continue
        act = actives[key]
        inact = inactives[key]
        if len(act) < 3 or len(inact) < 3:
            continue
        _, p_value = mannwhitneyu(act, inact)
        act_med = np.median(act)
        inact_med = np.median(inact)
        if p_value < 0.05:
            if act_med > inact_med:
                units_directions[key[1:]] += [1]
            else:
                units_directions[key[1:]] += [-1]
    actives_union = collections.defaultdict(list)
    inactive_union = collections.defaultdict(list)
    for k, v in actives.items():
        actives_union[k[1:]] += v
    for k, v in inactives.items():
        inactive_union[k[1:]] += v
    all_keys_union = set(actives_union.keys()).union(set(inactive_union.keys()))
    for key in tqdm(all_keys_union):
        if key not in actives_union:
            continue
        if key not in inactive_union:
            continue
        act = actives_union[key]
        inact = inactive_union[key]
        if len(act) < 3 or len(inact) < 3:
            continue
        _, p_value = mannwhitneyu(act, inact)
        act_med = np.median(act)
        inact_med = np.median(inact)
        if p_value < 0.05:
            if act_med > inact_med:
                units_directions[key] += [1]
            else:
                units_directions[key] += [-1]
    directions = {}
    for k, v in units_directions.items():
        if len(v) < 3:
            continue
        v = list(set(v))
        if len(v) == 1:
            directions[k] = v[0]
        else:
            continue
    print(directions)
    return directions


def standard_units_resolver_with_statistical_test(file_path, directions, rewrite=False):
    df = pd.read_csv(file_path, dtype=str)
    columns = list(df.columns)
    if "activity_direction" not in columns:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    if rewrite:
        df["activity_direction"] = ""
        df.to_csv(file_path, index=False)
    current_directions = df["activity_direction"].tolist()
    for i, v in enumerate(df[["standard_type", "standard_units"]].values):
        v = tuple(v)
        if v in directions:
            current_directions[i] = str(directions[v])
    df["activity_direction"] = current_directions
    df.to_csv(file_path, index=False)


if __name__ == "__main__":
    annotate_activities_with_activity_flag()
    act, inact = create_annotated_actives_and_inactives_dictionaries()
    joblib.dump(act, os.path.join(TMPDIR, "actives.pkl"))
    joblib.dump(inact, os.path.join(TMPDIR, "inactives.pkl"))
    act = joblib.load(os.path.join(TMPDIR, "actives.pkl"))
    inact = joblib.load(os.path.join(TMPDIR, "inactives.pkl"))
    directions = assign_directions_with_statistical_test(act, inact)
    file_path = os.path.join(CONFIGPATH, "activity_std_units.csv")
    standard_units_resolver_with_statistical_test(file_path, directions, rewrite=False)
    os.remove(os.path.join(TMPDIR, "actives.pkl"))
    os.remove(os.path.join(TMPDIR, "inactives.pkl"))
