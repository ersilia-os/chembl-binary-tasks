import os
import sys
import pandas as pd
import numpy as np
from chembl_molecule_sampler import ChemblMoleculeSampler

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from default import MIN_SIZE_ASSAY_TASK, MIN_SIZE_PROTEIN_TASK, MIN_COUNT_POSITIVE_CASES
from default import TOP_ASSAYS, TOP_PROTEINS
from default import CHEMBL_PWD, CHEMBL_USR, DATAPATH

class AllDatasets():
    def __init__(self, pathogen):
        self.pathogen = pathogen

    def _sample_negatives(self, num_molecules, list_positive_molecules):
        sampler = ChemblMoleculeSampler(DATAPATH, CHEMBL_USR, CHEMBL_PWD)
        neg_df = sampler.negative_sample(num_molecules, list_positive_molecules)
        return neg_df
    
    def create_datasets(self, df, dataset_name):
        target_cut = ["activity_lc", "activity_hc"]
        summary_info = {"activity_lc": f'{dataset_name}, low cutoff',
                        "activity_hc": f'{dataset_name}, high cutoff',}
        datasets = {}
        info_datasets = {}
        for tc in target_cut:
            df_ = df.rename(columns={tc:'activity'})[['compound_chembl_id', 'canonical_smiles', 'activity']]
            total_mols = len(df_)
            df_ = df_.groupby(['compound_chembl_id', 'canonical_smiles'])['activity'].mean().reset_index() #two instances of same molecule will be averaged
            df_['activity'] = df_['activity'].apply(lambda x: 1 if x>=0.5 else 0)
            non_dupl = len(df_)
            count_pos = len(df_[df_.activity==1])
            count_neg = len(df_[df_.activity==0])        
            print("Total Positives", count_pos, "Total Negatives",count_neg)
            # If there are fewer negatives than positives, fill up negatives with a random sample
            neg_df = None
            if count_neg < count_pos:
                neg_df = self._sample_negatives(
                        num_molecules=count_pos-count_neg,
                        list_positive_molecules=list(df_[df_.activity==1].compound_chembl_id))
                neg_df.rename(columns={'chembl_id':'compound_chembl_id'}, inplace=True)
                neg_df['activity'] = 0
                print("Added Random Negatives: ", len(neg_df))
                df_ = pd.concat([df_, neg_df], ignore_index=True)
            # Provide dataset only if minimum number of positive cases achieved.
            if count_pos >= MIN_COUNT_POSITIVE_CASES:
                df_.rename(columns={'canonical_smiles':'smiles'}, inplace=True)
                datasets[tc] = df_
                info_datasets[f'Dataset {summary_info[tc]} duplicated molecules:'] = total_mols-non_dupl
                info_datasets[f'Dataset {summary_info[tc]} positive molecules:'] = count_pos
                info_datasets[f'Dataset {summary_info[tc]} negative molecules:'] = count_neg
                info_datasets[f'Dataset {summary_info[tc]} negatives from ChEMBL added:'] = 0
                if neg_df is not None:
                    info_datasets[f'Dataset {summary_info[tc]} negatives from ChEMBL added:'] = len(neg_df)
                info_datasets[f'Dataset {summary_info[tc]} total data size']= len(df_)
            else:
                datasets[tc] = np.nan
                print(f"Not Enough Positives in Dataset, only {count_pos}")
                info_datasets[f'Dataset {summary_info[tc]} total data size']= len(df_)
                info_datasets[f"Dataset {summary_info[tc]} has too few positives:"]= count_pos
        return datasets["activity_lc"], datasets["activity_hc"], info_datasets

    def run_any(self, df):
        lc, hc, info_datasets = self.create_datasets(df, "All")
        return lc, hc, info_datasets

class OrgDatasets(AllDatasets):
    def __init__(self, pathogen):
        self.pathogen = pathogen

    def run_any_org(self, df):
        df = df[df.target_type == 'ORGANISM']
        lc, hc, info_datasets = self.create_datasets(df, "Org All")
        return lc, hc, info_datasets

    def top_assay_selector(self, df):
        assay_counts = df.assay_id.value_counts()[:TOP_ASSAYS]
        assay_counts = assay_counts[assay_counts >= MIN_SIZE_ASSAY_TASK]
        list_top_assays = list(assay_counts.index)
        print('Selected top assays:', list_top_assays)
        return list_top_assays
    
    def run_top_org(self, df):
        df = df[df.target_type == 'ORGANISM']
        top_assays = self.top_assay_selector(df)
        assays_data = {}
        if len(top_assays) != 0:
            for i, a in enumerate(top_assays):
                df_ = df[df.assay_id==a]
                lc, hc, info_datasets = self.create_datasets(df_,f'Org Top_{i}, (assay id:{a})' )
                assays_data[f'top_{i}'] = lc, hc, info_datasets
        return assays_data
        
class ProteinDatasets(AllDatasets):
    def __init__(self, pathogen):
        self.pathogen = pathogen

    def top_protein_selector(self, df):
        protein_counts = df.assay_id.value_counts()[:TOP_PROTEINS]
        protein_counts = protein_counts[protein_counts >= MIN_SIZE_PROTEIN_TASK]
        list_top_proteins = list(protein_counts.index)
        print('Selected top proteins:', list_top_proteins)
        return list_top_proteins

    def run_any_prot(self, df):
        df = df[df.target_type.str.contains('PROTEIN')]
        lc, hc, info_datasets = self.create_datasets(df, "Prot All")
        return lc, hc, info_datasets
    
    def run_top_prot(self, df):
        df = df[df.target_type.str.contains('PROTEIN')]
        top_prots = self.top_protein_selector(df)
        prots_data = {}
        if len(top_prots) != 0:
            for i,p in enumerate(top_prots):
                t = df[df["assay_id"]==p]["target_pref_name"].unique()[0]
                df_ = df[df.assay_id==p]
                lc, hc, info_datasets = self.create_datasets(df_, f'Prot Top_{i}, (assay id: {p}), (target: {t})')
                prots_data[f'top_{i}'] = lc, hc, info_datasets
        return prots_data

class BioassayDatasets(AllDatasets):
    def __init__(self, pathogen):
        self.pathogen = pathogen

    def run_merged_bioassays(self, df, st_type):
        df = df[df.target_type == 'ORGANISM']
        df = df[df.standard_type == st_type]
        lc, hc, info_datasets = self.create_datasets(df, f'Bioassay {st_type} merged')
        return lc, hc, info_datasets
    
