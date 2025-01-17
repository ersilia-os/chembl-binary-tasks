from rdkit import Chem
from rdkit.Chem import AllChem
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import random
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import argparse

random_seed = 54

np.random.seed(random_seed)
random.seed(random_seed)

parser = argparse.ArgumentParser("Estimate modelability of datasets")
parser.add_argument("--pathogen_code", type=str)
parser.add_argument("--output_dir", type=str)
args = parser.parse_args()

pathogen_code = args.pathogen_code
data_dir = args.output_dir
tasks_dir = os.path.join(data_dir, pathogen_code, "02_raw_tasks")

def get_binary_fingerprints_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024), dtype=int)
    return fp

ik_smi_pairs = []
for f in os.listdir(tasks_dir):
    print("Reading task", f)
    df = pd.read_csv(os.path.join(tasks_dir, f))
    for i, row in df.iterrows():
        ik_smi_pairs.append((row['inchikey'], row['smiles']))
ik_smi_pairs = list(set(ik_smi_pairs))

X = np.zeros((len(ik_smi_pairs), 1024), dtype=int)
keys = []
for i, (ik, smi) in tqdm(enumerate(ik_smi_pairs)):
    X[i] = get_binary_fingerprints_from_smiles(smi)
    keys.append(ik)
np.save(os.path.join(data_dir, pathogen_code, "03_fingerprints.npy"), X)
with open(os.path.join(data_dir, pathogen_code, "03_fingerprints_inchikeys.txt"), "w") as f:
    for k in keys:
        f.write(k + "\n")

def load_fingerprints(data_dir):
    X = np.load(os.path.join(data_dir, pathogen_code, "03_fingerprints.npy"))
    with open(os.path.join(data_dir, pathogen_code, "03_fingerprints_inchikeys.txt"), "r") as f:
        keys = f.read().splitlines()
    return X, keys

X, inchikeys = load_fingerprints(data_dir)

def modelability(df, X, inchikeys):
    inchikeys_ = list(df['inchikey'])
    columns = list(df.columns)
    assert len(columns) == 3, "The dataframe must have 3 columns"
    y = np.array(df[columns[-1]], dtype=int)
    indices = {}
    for i, ik in enumerate(inchikeys):
        indices[ik] = i
    idxs = [indices[ik] for ik in inchikeys_]
    X = X[idxs]
    print("Ready to model dataset with {0} samples".format(X.shape[0]))
    skf = StratifiedKFold(n_splits=5, shuffle=True)
    aurocs = []
    for train, test in tqdm(skf.split(X, y)):
        clf = RandomForestClassifier(n_estimators=100)
        print("Fitting model")
        clf.fit(X[train], y[train])
        aurocs += [roc_auc_score(y[test], clf.predict_proba(X[test])[:, 1])]
        print("AUROC", aurocs[-1])
    results = {"auroc_avg": np.mean(aurocs),
               "auroc_std": np.std(aurocs),
               "num_samples": X.shape[0],
               "num_pos_samples": np.sum(y),}
    return results

R = []
for l in os.listdir(tasks_dir):
    print("Modeling task", l)
    df = pd.read_csv(os.path.join(tasks_dir, l))
    results = modelability(df, X, inchikeys)
    fname = l[:-4]
    R += [(fname, results["auroc_avg"], results["auroc_std"], results["num_samples"], results["num_pos_samples"])]

pd.DataFrame(R, columns=["task", "auroc_avg", "auroc_std", "num_samples", "num_pos_samples"]).to_csv(os.path.join(data_dir, pathogen_code, "04_modelability.csv"), index=False)

