import os
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors

root = os.path.dirname(os.path.abspath(__file__))

import sys
sys.path.append(root)
from default import CONFIGPATH


def main(pathogen_code, data_dir):
    
    print(f"Cleaning data for pathogen {pathogen_code}")
    df = pd.read_csv(f"{data_dir}/{pathogen_code}/00_{pathogen_code}_original.csv", low_memory=False)
    initial_len = len(df)
    print(f"Initial length: {initial_len}")
    
    print("Reading all activity data from the config of this repo")
    da = pd.read_csv(os.path.join(CONFIGPATH, "all_activities.csv"))
    activities = {}
    for v in da.values:
        activities[v[0]] = (v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8])
    
    print("Reading conversions")
    dc = pd.read_csv("../config/standard_units_conversions.csv")
    dc = dc[dc["conversion_formula"].notnull()]
    dc = dc[dc["new_unit"].notnull()]
    conversions = {}
    for v in dc[["standard_units", "conversion_formula", "new_unit"]].values:
        conversions[v[0]] = (v[1], v[2])
    
    print("Cleaning data...")
    R = []
    for _, v in tqdm(df.iterrows()):
        activity_id = v["activity_id"]
        if activity_id  not in activities:
            continue
        pchembl_value = activities[activity_id][0]
        standard_relation = activities[activity_id][1]
        standard_value = activities[activity_id][2]
        standard_type = activities[activity_id][3]
        standard_unit = activities[activity_id][4]
        activity_flag = activities[activity_id][5]
        direction_flag = activities[activity_id][6]
        direction_confidence = activities[activity_id][7]
        if standard_unit in conversions:
            conversion_formula, new_unit = conversions[standard_unit]
            data = {"x": standard_value}
            if "w" in conversion_formula:
                w = Descriptors.ExactMolWt(Chem.MolFromSmiles(v["canonical_smiles"]))
                data["w"] = w
            standard_value = eval(conversion_formula, {}, data)
            standard_unit = new_unit
        smiles = v["parent_canonical_smiles"]
        inchikey = v["parent_inchikey"]
        assay_id = v["assay_id"]
        assay_type = v["assay_type"]
        target_type = v["target_type"]
        target_id = v["target_chembl_id"]
        assay_chembl_id = v["assay_chembl_id"]
        data = {
            "activity_id": activity_id,
            "assay_id": assay_chembl_id,
            "assay_type": assay_type,
            "target_type": target_type,
            "target_id": target_id,
            "inchikey": inchikey,
            "smiles": smiles,
            "standard_relation": standard_relation,
            "standard_value": standard_value,
            "standard_units": standard_unit,
            "standard_type": standard_type,
            "pchembl_value": pchembl_value,
            "activity_flag": activity_flag,
            "direction_flag": direction_flag,
            "direction_confidence": direction_confidence,
        }
        R += [data]
    df = pd.DataFrame(R)

    print("Working on pChEMBL values...")
    pchembl_values = []
    pchembl_relations = []
    for _, v in df.iterrows():
        pchembl_value = v["pchembl_value"]
        standard_value = v["standard_value"]
        standard_relation = v["standard_relation"]
        standard_units = v["standard_units"]
        if standard_units is None or str(standard_units).lower() == "nan":
            print("Standard units is nan")
            pchembl_values += [None]
            pchembl_relations += [None]
            continue
        if str(pchembl_value).lower() != "nan":
            pchembl_values += [pchembl_value]
            if standard_units.endswith("M"):
                if standard_relation == "<":
                    pchembl_relations += [">"]
                elif standard_relation == ">":
                    pchembl_relations += ["<"]
                else:
                    pchembl_relations += ["="]
            else:
                if standard_relation == "=":
                    pchembl_relations += ["="]
                else:
                    pchembl_relations += [None]
        else:
            if standard_units == "uM":
                value = standard_value * 1e-6
                pchembl_value = np.clip(-np.log10(value), 1, 9)
                pchembl_values += [pchembl_value]
                if standard_relation == "<":
                    pchembl_relations += [">"]
                elif standard_relation == ">":
                    pchembl_relations += ["<"]
                else:
                    pchembl_relations += ["="]
            else:
                pchembl_values += [None]
                pchembl_relations += [None]

    df["pchembl_value"] = pchembl_values
    df["pchembl_relation"] = pchembl_relations

    columns = [
        "activity_id",
        "assay_id",
        "assay_type",
        "target_type",
        "target_id",
        "inchikey",
        "smiles",
        "standard_relation",
        "standard_value",
        "standard_units",
        "standard_type",
        "pchembl_relation",
        "pchembl_value",
        "activity_flag",
        "direction_flag",
        "direction_confidence"
    ]
    df = df[columns]
    
    print("Final length: ", len(df))
    print("Saving cleaned data to {0}".format(f"{data_dir}/{pathogen_code}/{pathogen_code}_cleaned.csv"))
    df.to_csv(f"{data_dir}/{pathogen_code}/01_{pathogen_code}_cleaned.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Clean pathogen data from ChEMBL database.')
    parser.add_argument('-d', '--output_dir', type=str, required=True, help='Output directory')
    parser.add_argument('-p', '--pathogen_code', type=str, required=False, help='Pathogen code to search for')
    args = parser.parse_args()
    main(args.pathogen_code, args.output_dir)