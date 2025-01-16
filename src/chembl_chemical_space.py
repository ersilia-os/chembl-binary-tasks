import os
import csv
import sys
from rdkit import Chem
from tqdm import tqdm

root = os.path.dirname(os.path.abspath(__file__))

sys.path.append(root)

from default import CONFIGPATH

TMPDIR = os.path.join(root, "..", "tmp")

def get_num_heavy_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return mol.GetNumHeavyAtoms()

def get_inchikey(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)

def main():
    R = []
    with open(os.path.join(TMPDIR, "chembl_35_chemreps.txt"), "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for r in tqdm(reader):
            chembl_id = r[0]
            smiles = r[1]
            num_heavy_atoms = get_num_heavy_atoms(smiles)
            inchikey = get_inchikey(smiles)
            R += [[chembl_id, smiles, num_heavy_atoms, inchikey]]
    with open(os.path.join(CONFIGPATH, "all_molecules.csv"), "w") as f:
        writer = csv.writer(f)
        writer.writerow(["chembl_id", "smiles", "num_heavy_atoms", "inchikey"])
        for r in R:
            writer.writerow(r)


if __name__ == "__main__":
    main()