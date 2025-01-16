import os
import sys
import pandas as pd
from tqdm import tqdm

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)

from default import CONFIGPATH, TMPDIR

def main():
    print("Loading activities table...")
    df = pd.read_csv(os.path.join(TMPDIR, "activities.csv"), low_memory=False)
    print(df.columns)
    print("Getting activity comments...")
    dc = pd.read_csv(os.path.join(CONFIGPATH, "activity_comments.csv"))
    activity_comments = {}
    for v in dc[["activity_comment", "activity_classified"]].values:
        activity_comments[v[0]] = v[1]
    print("Getting standard text...")
    ds = pd.read_csv(os.path.join(CONFIGPATH, "standard_text.csv"))
    standard_text = {}
    for v in ds.values:
        if v[1] == 0:
            continue
        standard_text[v[0]] = v[1]
    print("Getting directions from lookup table...")
    dl = pd.read_csv(os.path.join(CONFIGPATH, "activity_stds_lookup.csv"))
    directions_lookup = {}
    for v in dl[["standard_type", "standard_units", "activity_direction"]].values:
        directions_lookup[(v[0], v[1])] = v[2]
    print("Getting directions from standard units table...")
    du = pd.read_csv(os.path.join(CONFIGPATH, "activity_std_units.csv"))
    directions = {}
    for v in du[["standard_type", "standard_units", "activity_direction"]].values:
        directions[(v[0], v[1])] = v[2]
    R = []
    columns = ["activity_id",
               "standard_relation",
               "standard_value",
               "standard_type",
               "standard_units",
               "pchembl_value",
               "activity_comment",
               "standard_text_value"]
    R = []
    for v in tqdm(df[columns].values):
        activity_id = v[0]
        standard_relation = v[1]
        standard_value = v[2]
        standard_type = v[3]
        standard_units = v[4]
        pchembl_value = v[5]
        activity_comment = v[-2]
        standard_text_value = v[-1]
        if standard_text_value in standard_text:
            activity_flag = standard_text[standard_text_value]
        elif activity_comment in activity_comments:
            activity_flag = activity_comments[activity_comment]
        else:
            activity_flag = 0
        if (standard_type, standard_units) in directions_lookup:
            direction_flag = directions_lookup[(standard_type, standard_units)]
            direction_confidence = 1
        elif (standard_type, standard_units) in directions:
            direction_flag = directions[(standard_type, standard_units)]
            direction_confidence = 0
        else:
            direction_flag = 0
            direction_confidence = 0
        r = [
            activity_id,
            pchembl_value,
            standard_relation,
            standard_value,
            standard_type,
            standard_units,
            activity_flag,
            direction_flag,
            direction_confidence,
        ]
        R += [r]
    print("Creating activities table...")
    df = pd.DataFrame(R, columns=["activity_id",
                                  "pchembl_value",
                                  "standard_relation",
                                  "standard_value",
                                  "standard_type",
                                  "standard_units",
                                  "activity_flag",
                                  "direction_flag",
                                  "direction_confidence"])
    df = df[df["standard_value"].notnull()]
    df = df[df["standard_type"].notnull()]
    df = df[df["standard_units"].notnull()]
    df.loc[df["standard_relation"].isnull(), "standard_relation"] = "="
    df.loc[df["standard_relation"] == "<=", "standard_relation"] = "<"
    df.loc[df["standard_relation"] == ">=", "standard_relation"] = ">"
    df.loc[df["standard_relation"] == "~", "standard_relation"] = "="
    df.loc[df["standard_relation"] == "<<", "standard_relation"] = "<"
    df.loc[df["standard_relation"] == ">>", "standard_relation"] = ">"
    activity_flags = [int(x) for x in df["activity_flag"].tolist()]
    df["activity_flag"] = activity_flags
    df.to_csv(os.path.join(CONFIGPATH, "all_activities.csv"), index=False)

    standard_units_file = os.path.join(CONFIGPATH, "standard_units_conversions.csv")
    if os.path.exists(standard_units_file):
        su_cache = {}
        du = pd.read_csv(standard_units_file)
        for v in du[["standard_units", "count", "conversion_formula", "new_unit"]].values:
            if v[1]:
                su_cache[v[0]] = (v[2], v[3])
    else:
        su_cache = {}
    R = []
    for k, v in df.value_counts("standard_units").items():
        if k in su_cache:
            R += [[k, v, su_cache[k][0], su_cache[k][1]]]
        else:
            R += [[k, v, "", ""]]
    pd.DataFrame(R, columns=["standard_units", "count", "conversion_formula", "new_unit"]).to_csv(os.path.join(CONFIGPATH, "standard_units_conversions.csv"), index=False)

if __name__ == "__main__":
    main()