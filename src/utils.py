import os
import sys
from tqdm import tqdm
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from default import UNITS_MASTER, CONFIGPATH


class RawCleaner():
    def __init__(self):
        self.ucum = UNITS_MASTER
    
    def get_standard_text(self, df):
        """The variable standard_text contains information about the type of result obtained in the experiment.
        We will use this information to assign a value to standard text:
        The column standard_text_active will be:
           - 1 if active.
           - 0 if inactive.
           - Empty otherwise.
        """
        sc = pd.read_csv(os.path.join(CONFIGPATH, "standard_text_classified.csv"))
        standard_text_dict = {}
        for k, v in sc.values:
            for x in k.split(";"):
                standard_text_dict[x] = v
        standard_text_active = []
        for c in tqdm(df['standard_text'].tolist()):
            c = str(c)
            if c is None or c == 'nan':
                standard_text_active.append(np.nan)
            else:
                if c not in standard_text_dict:
                    raise Exception(f"Standard text {c} not found in the standard text file. Please check the standard_text_classified.csv file and update it accordingly.")
                else:
                    val = standard_text_dict[c]
                    if val == "":
                        standard_text_active.append(np.nan)
                    else:
                        standard_text_active.append(int(val))
        df['standard_text_active'] = standard_text_active
        return df
    
    def get_comment(self, df):
        """The variable comment_active may sometimes include an indication that the result 
        of the experiment is "Active" or "Not Active".
        We will use this information in the calculation of the target variables.
        The column comment_active will be:
            - 1 if activity_comment = 'active', etc.
            - 0 if activity_comment = 'not active', 'inactive, etc.
            - Empty otherwise
        """
        ac = pd.read_csv(os.path.join(CONFIGPATH, "activity_comments_classified.csv"))
        comment_dict = {}
        for k, v in ac.values:
            for x in k.split(";"):
                comment_dict[x] = v
        comment_active = []
        for c in tqdm(df['activity_comment'].tolist()):
            c = str(c)
            if c is None or c == 'nan':
                comment_active.append(np.nan)
            else:
                if c not in comment_dict:
                    raise Exception(f"Comment {c} not found in the activity comments file. Please check the activity_comments_classified.csv file and update it accordingly.")
                else:
                    val = comment_dict[c]
                    if val == "":
                        comment_active.append(np.nan)
                    else:
                        comment_active.append(int(val))
        df['comment_active'] = comment_active
        return df

    def eliminate_rows(self, df):
        pref = len(df)
        df = df[~df["canonical_smiles"].isna()]
        df = df[~df["parent_canonical_smiles"].isna()]
        print("removing rows with empty smiles: {}".format(pref-len(df)))
        df['standard_units'] = df['standard_units'].fillna('N/A')    
        # Remove rows where both standard_value and comment_active are null 
        # (we can't know if they are active or not)
        rows_before_filter = len(df)
        df = df[~((df.standard_value.isnull()) & (df.comment_active.isnull()) )]
        rows_after_filter = len(df)
        no_activity_info = rows_before_filter-rows_after_filter
        print('\nRemoving rows where standard_value is null and activity_comment does not inform on activity.')
        print(f'Removed {no_activity_info } rows. Remaining {rows_after_filter} rows.')
        return df, no_activity_info 
    
    def drop_unwanted_cols(self, df):
        cols_to_drop = [
            'doc_id',
            'activity_id',
            'assay_bao_format',
            'target_tax_id',
            'assay_description']
        cols_int = list(set(cols_to_drop).intersection(df.columns))
        df = df.drop(columns=cols_int, errors='ignore')
        return df

    def _mol_weight(self, smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return Descriptors.MolWt(mol)

    def add_mol_weight(self, df, is_parent=False):
        if not is_parent:
            print("Molecular weight calculation for base SMILES")
            df["molecular_weight"] = df["canonical_smiles"].apply(self._mol_weight)
            #if molweight cannot be calculated, drop row
            smi_not_processed = len(df[df["molecular_weight"].isna()])
            print(len(df), smi_not_processed)
            df = df[~df["molecular_weight"].isna()]
            print(len(df))
            return df, smi_not_processed
        else:
            print("Molecular weight calculation for base SMILES for parent")
            df["parent_molecular_weight"] = df["parent_canonical_smiles"].apply(self._mol_weight)
            #if molweight cannot be calculated, drop row
            smi_not_processed = len(df[df["parent_molecular_weight"].isna()])
            print(len(df), smi_not_processed)
            df = df[~df["parent_molecular_weight"].isna()]
            print(len(df))
            return df, smi_not_processed
    
    def activity_cleanup(self, df):
        for i,r in df.iterrows():
            if r['standard_type'] == 'log10cfu':
                df.at[i, 'standard_type'] = 'Activity'
                df.at[i, 'standard_units'] = 'log10CFU'
            elif r['standard_type'] == 'log10CFU/ml':
                df.at[i, 'standard_type'] = 'Activity'
                df.at[i, 'standard_units'] = 'log10CFU/ml'
            elif r['standard_type'] == '-logMIC':
                df.at[i, 'standard_type'] = 'MIC'
                df.at[i, 'standard_units'] = '-logMIC'
            elif r["standard_type"] == "INHIBITION":
                df.at[i, "standard_type"] = "Inhibition"
            else:
                continue
        return df
    
    def run(self, df):
        print("Getting standard text data")
        df = self.get_standard_text(df)
        print("Getting comment data")
        df = self.get_comment(df)
        print("Eliminating no info rows")
        df, no_activity_info = self.eliminate_rows(df)
        print("Dropping unwanted columns")
        df = self.drop_unwanted_cols(df)
        df, smi_not_processed_0 = self.add_mol_weight(df, is_parent=False)
        df, smi_not_processed_1 = self.add_mol_weight(df, is_parent=True)
        smi_not_processed = list(set(smi_not_processed_0 + smi_not_processed_1))
        df = self.activity_cleanup(df)
        return df, no_activity_info, smi_not_processed

class UnitStandardiser():
    #the units file has been processed in UCUM. Selected units will be converted to standard:

    def __init__(self):
        self.unit_col = "standard_units"
        self.value_col = "standard_value"
        self.mw_col = "molecular_weight"
        self.ucum = UNITS_MASTER

    def add_ucum_units(self, df):
        unit_cols = ["val_unit", "final_unit", "transformer"]
        for c in unit_cols:
            units_to_col_units = dict(zip(self.ucum['units'], self.ucum[c]))
            df[c] = [units_to_col_units.get(unit, unit) for unit in df['standard_units']]
        return df
    
    def transformer(self, df):
        unique_transformers = df["transformer"].dropna().unique()  # Drop NaN values and get unique transformers
        unique_transformers = unique_transformers[unique_transformers != 'N/A'] 
        transformer_mapping = {}
        for t in unique_transformers:
            try:
                if "molecular_weight" in t:
                    transformer_expr = t.replace('standard_value', 'row["standard_value"]').replace('molecular_weight', 'row["molecular_weight"]')
                    function_expr = f"lambda row: {transformer_expr}"
                    transformer_mapping[t] = eval(function_expr)
                else:
                    transformer_expr = t.replace('standard_value', 'row["standard_value"]')
                    function_expr = f"lambda row: {transformer_expr}"
                    transformer_mapping[t] = eval(function_expr)
            except:
                print(f"Value {t} not in the UCUM file")
        # Apply the appropriate function to each row based on the value in the 'transformer' column
        def apply_transformer(row):
            transformer_func = transformer_mapping.get(row['transformer'])
            if transformer_func is not None:
                result = transformer_func(row)
                if result is not None:
                    return f"{result:.3f}"  # Format the result to have at least three digits
                else:
                    return None
            else:
                if f"{row['transformer']}" != 'N/A':
                    print(f"Warning: Transformer '{row['transformer']}' not found in mapping. Skipping...")
                return None
        df['final_value'] = df.apply(apply_transformer, axis=1)
        return df

    def run(self, df):
        df = self.add_ucum_units(df)
        df = self.transformer(df)
        return df

class Binarizer():
    def __init__(self):
        pass
    
    def _load_config_file(self):
        df = pd.read_csv(os.path.join(CONFIGPATH, "cutoffs_manual.csv"))
        df['final_unit'] = df['final_unit'].fillna('N/A')  
        df = df[df["use"]==1]
        print('\nstandard_type_config shape:', df.shape)
        return df
    
    def _warning_for_missing_config_entries(self, df):
        """Give a warning if missing important (type,units) combination in config table"""
        # Top 10 combinations of (type, units) that appear at least 100 times
        top_values_type_unit = df[['standard_type', 'final_unit']].value_counts()[0:10]
        top_values_type_unit = top_values_type_unit[top_values_type_unit>100]
        # Detailed data for those (type, units)
        df_top = df.merge(top_values_type_unit.to_frame(), left_on=['standard_type', 'final_unit'], right_index=True)
        # Combinations not present in config table ()
        missing_top_type_units = df_top[df_top.is_in_cutoff_table=='left_only']\
                [['standard_type', 'final_unit']].value_counts()
        if len(missing_top_type_units) > 0:
            print('\n--- WARNING - The following combinations of (standard_type, final_unit) are'
                ' often present in the data but do not exist in the configuration table. Please'
                ' consider updating the configuration table:')
            print(missing_top_type_units)
    
    def add_cutoffs(self, df):
        configfile = self._load_config_file()     
        configfile = configfile[["standard_type", "final_unit","active_direction", 'low_cut', 'high_cut']]
        df1 = df.merge(configfile, how='left',
              on=['standard_type', 'final_unit'],
              indicator='is_in_cutoff_table')
        self._warning_for_missing_config_entries(df1)
        rows_before_filter = len(df1)
        df1 = df1[df1["is_in_cutoff_table"]=="both"]
        rows_after_filter = len(df1)
        assay_entries_discarded = rows_before_filter-rows_after_filter
        print('\nRemoving rows where the combination of standard_type and final_unit is not in the cutoff table.')
        print(f'Removed {assay_entries_discarded} rows. Remaining {rows_after_filter} rows.')
        return df1, assay_entries_discarded

    def _calculate_active(self, row, cut):
        if not pd.isna(row.comment_active): #priority to author comment of Active or Inactive
            return row.comment_active
        elif row.active_direction == 1:  # Higher value is more active
            if row.final_value >= cut:
                return 1
            else:
                return 0
        elif row.active_direction == -1:  # Lower value is more active
            if row.final_value <= cut:

                return 1
            else:
                return 0
    
    def calculate_active(self,df):
        df["final_value"]=df["final_value"].astype(float)    
        df['activity_lc'] = df.apply(
                lambda row: self._calculate_active(row, row.low_cut), 
                axis=1).astype('float')
        df['activity_hc'] = df.apply(
                lambda row: self._calculate_active(row, row.high_cut), 
                axis=1).astype('float')
        return df
    
    def remove_inconsistencies(self, df):
        # Remove cases where standard_direction is not consistent with our activity definition
    # Background:
    # If the value resulting of an experiment is beyond the range that can be measured, 
    # instead of reporting the value, it will be reported as ">x" or "<x".

    # The variable standard_direction contains "=" if the precise value is reported. It will contain "<", "<=", ">" or ">=" if the reported value is a lower or upper bound. For example, if the "real" value is 500 but only up to 100 can be measured, then standard_value=100 and standard_direction=">".

    # Taking this into account, it makes sense that, for results that we label as ACTIVE:
    # - If active_direction=1
    #     * standard_relation may be '>' (it indicates a "large" value)
    # - If active_direction=-1
    #     * standard_relation may be '<' (it indicates a "small" value)

    # For results that we label as NOT ACTIVE:
    # - If active_direction=1
    #     * standard_relation may be '<' (it indicates a "small" value)
    # - If active_direction=-1
    #     * standard_relation may be '>' (it indicates a "large" value)
        total_inconsistent = 0
        rows_to_drop = df[(df.comment_active.isnull()) &
        (df.activity_hc==0) & 
        (df.active_direction==1) & 
        (df.standard_relation.isin(['>','>=']))
        ].index
        df.drop(rows_to_drop, inplace=True)
        print(f'Removed {len(rows_to_drop)} cases with active direction +, relation ">", but labeled as not active')
        total_inconsistent = total_inconsistent+len(rows_to_drop)

        rows_to_drop = df[(df.comment_active.isnull()) &
        (df.activity_hc==0) & 
        (df.active_direction==-1) & 
        (df.standard_relation.isin(['<','<=']))
        ].index
        df.drop(rows_to_drop, inplace=True)
        print(f'Removed {len(rows_to_drop)} cases with active direction -, relation "<", but labeled as not active')
        total_inconsistent = total_inconsistent+len(rows_to_drop)

        rows_to_drop = df[(df.comment_active.isnull()) &
        (df.activity_hc==1) & 
        (df.active_direction==1) & 
        (df.standard_relation.isin(['<','<=']))
        ].index
        df.drop(rows_to_drop, inplace=True)
        print(f'Removed {len(rows_to_drop)} cases with active direction +, relation "<", but labeled as active')
        total_inconsistent = total_inconsistent+len(rows_to_drop)

        rows_to_drop = df[(df.comment_active.isnull()) &
        (df.activity_hc==1) & 
        (df.active_direction==-1) & 
        (df.standard_relation.isin(['>','>=']))
        ].index
        df.drop(rows_to_drop, inplace=True)
        print(f'Removed {len(rows_to_drop)} cases with active direction -, relation ">", but labeled as active')
        total_inconsistent = total_inconsistent+len(rows_to_drop)

        print('Cases remaining after filter:', len(df))
        
        return df, total_inconsistent


    def run(self, df):
        df, assay_entries_discarded = self.add_cutoffs(df)
        df = self.calculate_active(df)
        df, total_inconsistent = self.remove_inconsistencies(df)
        return df, assay_entries_discarded, total_inconsistent