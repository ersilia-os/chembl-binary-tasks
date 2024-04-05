import os
import sys
import pandas as pd
import numpy as np
from chembl_molecule_sampler import ChemblMoleculeSampler

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from default import MIN_SIZE_ASSAY_TASK, MIN_SIZE_PROTEIN_TASK, MIN_COUNT_POSITIVE_CASES
from default import TOP_ASSAYS, TOP_PROTEINS, TOP_TYPES
from default import CHEMBL_PWD, CHEMBL_USR, DATAPATH

class TypeDatasets():
    def __init__(self, pathogen):
        self.pathogen = pathogen

    def top_type_selector(self, df):
        type_counts = df.standard_type.value_counts()
        list_top_types = list(type_counts.index)[:TOP_TYPES]
        print('Selected top types:', list_top_types)
        return list_top_types
    
    def _sample_negatives(self, num_molecules, list_positive_molecules):
        sampler = ChemblMoleculeSampler(DATAPATH, CHEMBL_USR, CHEMBL_PWD)
        neg_df = sampler.negative_sample(num_molecules, list_positive_molecules)
        return neg_df
    
    def create_datasets(self, df):
        target_var = ["activity_lc", "activity_hc"]
        datasets = {}
        for tv in target_var:
            df_ = df.rename(columns={tv:'activity'})[['compound_chembl_id', 'canonical_smiles', 'activity']]
            print(df_.shape)
            df_ = df_.groupby(['compound_chembl_id', 'canonical_smiles'])['activity'].mean().reset_index()
            df_['activity'] = df_['activity'].apply(lambda x: 1 if x>=0.5 else 0)
            count_pos = len(df_[df_.activity==1])
            count_neg = len(df_[df_.activity==0])        
            print(count_pos, count_neg)
            # If there are fewer negatives than positives, fill up negatives with a random sample
            if count_neg < count_pos:
                neg_df = self._sample_negatives(
                        num_molecules=count_pos-count_neg,
                        list_positive_molecules=list(df_[df_.activity==1].compound_chembl_id))
                neg_df.rename(columns={'chembl_id':'compound_chembl_id'}, inplace=True)
                neg_df['activity'] = 0
                df_ = pd.concat([df_, neg_df], ignore_index=True)
            # Provide dataset only if minimum number of positive cases achieved.
            if count_pos >= MIN_COUNT_POSITIVE_CASES:
                df_.rename(columns={'canonical_smiles':'smiles'}, inplace=True)
                datasets[tv] = df_
            else:
                print("here")
                datasets[tv] = np.nan
        return datasets.values()

    def run_any(self, df):
        lc, hc = self.create_datasets(df)
        return lc, hc

class OrgDatasets(TypeDatasets):
    def __init__(self, pathogen):
        self.pathogen = pathogen

    def top_assay_selector(self, df):
        assay_counts = df.assay_id.value_counts()[:TOP_ASSAYS]
        print(assay_counts)
        assay_counts = assay_counts[assay_counts >= MIN_SIZE_ASSAY_TASK]
        list_top_assays = list(assay_counts.index)
        print('Selected top assays:', list_top_assays)
        return list_top_assays

    def run_any_org(self, df):
        df = df[df.target_type == 'ORGANISM']
        lc, hc = self.create_datasets(df)
        return lc, hc
    
    def run_top_org(self, df):
        df = df[df.target_type == 'ORGANISM']
        print("HERE", len(df))
        top_assays = self.top_assay_selector(df)
        assays_data = {}
        for a in top_assays:
            df_ = df[df.assay_id==a]
            lc, hc = self.create_datasets(df_)
            assays_data[a] = lc, hc
        return assays_data
        

class ProteinDatasets(TypeDatasets):
    def __init__(self, pathogen):
        self.pathogen = pathogen

    def top_protein_selector(self, df):
        protein_counts = df.target_pref_name.value_counts()[:TOP_PROTEINS]
        protein_counts = protein_counts[protein_counts >= MIN_SIZE_PROTEIN_TASK]
        list_top_proteins = list(protein_counts.index)
        print('Selected top proteins:', list_top_proteins)
        return list_top_proteins

    def run_any_prot(self, df):
        df = df[df.target_type.str.contains('PROTEIN')]
        lc, hc = self.create_datasets(df)
        return lc, hc
    
    def run_top_prot(self, df):
        df = df[df.target_type.str.contains('PROTEIN')]
        top_prots = self.top_protein_selector(df)
        prots_data = {}
        for p in top_prots:
            df_ = df[df.target_pref_name==p]
            lc, hc = self.create_datasets(df_)
            prots_data[p] = lc, hc
        return prots_data

class BioassayDatasets(TypeDatasets):
    def __init__(self, pathogen):
        self.pathogen = pathogen

    def run_merged_bioassays(self, df, st_type):
        df = df[df.target_type == 'ORGANISM']
        df = df[df.standard_type == st_type]
        lc, hc = self.create_datasets(df)
        return lc, hc
    
