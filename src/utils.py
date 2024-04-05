import os
import sys
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
    
    def get_comment(self, df):
        """The variable comment_active may sometimes include an indication that the result 
        of the experiment is "Active" or "Not Active".
        We will use this information in the calculation of the target variables.
        The variable comment_active will be:
            - True if activity_comment = 'active'
            - False if activity_comment = 'not active'
            - NA otherwise
        """
        df['comment_active'] = df['activity_comment'].str.upper()\
                .apply(lambda x: 
                        1 if x=='ACTIVE' else
                        0 if x=='NOT ACTIVE'
                        else np.nan
                    ).astype('float')
        return df

    def eliminate_rows(self, df):
        pref = len(df)
        df = df[~df["canonical_smiles"].isna()]
        print("removing rows with empty smiles: {}".format(pref-len(df)))
        df['standard_units'] = df['standard_units'].fillna('N/A')    
        # Remove rows where both standard_value and comment_active are null 
        # (we can't know if they are active or not)
        rows_before_filter = len(df)
        df = df[~((df.standard_value.isnull()) & (df.comment_active.isnull()) )]
        rows_after_filter = len(df)
        print('\nRemoving rows where standard_value is null and activity_comment does not inform on activity.')
        print(f'Removed {rows_before_filter-rows_after_filter} rows. Remaining {rows_after_filter} rows.')
        return df
    
    def drop_unwanted_cols(self, df):
        cols_to_drop = ['doc_id', 'activity_id', 'assay_type',
       'assay_confidence_score', 'assay_bao_format', 'pchembl_value', 'activity_comment',
       'target_tax_id', 'protein_accession_class', 'year',
       'pubmed_id', 'count_activity_rows', 'doc_id_all',
       'assay_id_all', 'activity_id_all', 'assay_description']
        cols_int = list(set(cols_to_drop).intersection(df.columns))
        df = df.drop(columns=cols_int, errors='ignore')
        return df

    def _mol_weight(self, smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return Descriptors.MolWt(mol)

    def add_mol_weight(self, df):
        df["molecular_weight"] = df["canonical_smiles"].apply(self._mol_weight)
        #if molweight cannot be calculated, drop row
        df = df[~df["molecular_weight"].isna()]
        return df
    
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
    
    def add_ucum_units(self, df):
        units_to_val_units = dict(zip(self.ucum['units'], self.ucum['val_unit']))
        #df['val_units'] = [units_to_val_units[unit] for unit in df['standard_units']]
        df['val_units'] = [units_to_val_units.get(unit, unit) for unit in df['standard_units']]
        return df
    
    def run(self, df):
        df = self.get_comment(df)
        df = self.eliminate_rows(df)
        df = self.drop_unwanted_cols(df)
        df = self.add_mol_weight(df)
        df = self.activity_cleanup(df)
        df = self.add_ucum_units(df)
        return df

class UnitStandardiser():
    #the units file has been processed in UCUM. Selected units will be converted to standard:
    #Molar: umol
    #weight/volume: ug.ml-1 to uM
    #weight/weight: ug.mg-1
    #molar/weight: umol/mg

    def __init__(self):
        self.unit_col = "standard_units"
        self.value_col = "standard_value"
        self.mw_col = "molecular_weight"
        self.ucum = UNITS_MASTER
    
    def _parse_function(self, s):
        if 'standard_value' not in s:
            return None
        if "molecular_weight" in s:
            p = "lambda x,y: "
        else:
            p = "lambda x: "
        s = s.replace("molecular_weight", "y")
        s = s.replace("standard_value", "x")
        s = p + s
        return eval(s)

    def _umol_converter(self):
        converter_str = {}
        converter_frm = {}
        for i,r in self.ucum.iterrows():
            if r["final_unit"] == "umol":
                converter_str[r["units"]] = r["transformer"]
        for k,v in converter_str.items():
            f = self._parse_function(v)
            converter_frm[k]=f
        return converter_frm

    def _ugmg_converter(self):
        converter_str = {}
        converter_frm = {}
        for i,r in self.ucum.iterrows():
            if r["final_unit"] == "ug.mg-1":
                converter_str[r["units"]] = r["transformer"]
        for k,v in converter_str.items():
            f = self._parse_function(v)
            converter_frm[k]=f
        return converter_frm
    
    def _umolmg_converter(self):
        converter_str = {}
        converter_frm = {}
        for i,r in self.ucum.iterrows():
            if r["final_unit"] == "umol.mg-1":
                converter_str[r["units"]] = r["transformer"]
        for k,v in converter_str.items():
            f = self._parse_function(v)
            converter_frm[k]=f
        return converter_frm

    def standardise(self,df):
        final_units = []
        final_value = []
        umol_converter = self._umol_converter()
        ugmg_converter = self._ugmg_converter()
        umolmg_converter = self._umolmg_converter()

        for i,r in df.iterrows():
            if r[self.unit_col] in umol_converter.keys():
                final_units += ["umol"]
                if r["val_units"] in ["umol", "nmol", "pmol", "mmol", "mol"]:
                    final_value += [umol_converter[r[self.unit_col]](r[self.value_col])] 
                else:
                    final_value += [umol_converter[r[self.unit_col]](r[self.value_col], r[self.mw_col])]
            elif r[self.unit_col] in ugmg_converter.keys():
                final_units += ["ug.mg-1"]
                final_value += [ugmg_converter[r[self.unit_col]](r[self.value_col])]
            elif r[self.unit_col] in umolmg_converter.keys():
                final_units += ["umol.mg-1"]
                final_value += [umolmg_converter[r[self.unit_col]](r[self.value_col])]
            else:
                final_units += [r[self.unit_col]]
                final_value += [r[self.value_col]]
        df["final_units"] = final_units
        df["final_value"] = final_value
        return df

class Binarizer():
    def __init__(self):
        pass
    
    def _load_config_file(self):
        df = pd.read_csv(os.path.join(CONFIGPATH, "cutoff_config.csv"),
                         usecols = ["standard_type", "final_units", "active_direction", 'low_cut', 'high_cut'],
                         keep_default_na = False, na_values=['']
                         )
        print('\nstandard_type_config shape:', df.shape)
        return df
    
    def _warning_for_missing_config_entries(self, df):
        """Give a warning if missing important (type,units) combination in config table"""
        # Top 10 combinations of (type, units) that appear at least 100 times
        top_values_type_unit = df[['standard_type', 'final_units']].value_counts()[0:10]
        top_values_type_unit = top_values_type_unit[top_values_type_unit>100]
        # Detailed data for those (type, units)
        df_top = df.merge(top_values_type_unit.to_frame(), left_on=['standard_type', 'final_units'], right_index=True)
        # Combinations not present in config table ()
        missing_top_type_units = df_top[df_top.is_in_config_table=='left_only']\
                [['standard_type', 'final_units']].value_counts()
        if len(missing_top_type_units) > 0:
            print('\n--- WARNING - The following combinations of (standard_type, standard_unit) are'
                ' often present in the data but do not exist in the configuration table. Please'
                ' consider updating the configuration table:')
            print(missing_top_type_units)

    def add_cutoffs(self, df):
        configfile = self._load_config_file()     
        df1 = df.merge(configfile, how='left',
              on=['standard_type', 'final_units'],
              indicator='is_in_config_table')
        self._warning_for_missing_config_entries(df1)
        print(df1.columns)
        return df1
    
    def eliminate_no_direction(self, df):
        rows_before_filter = len(df)
        df = df[df.active_direction.notnull()]
        rows_after_filter = len(df)
        print('\nRemoving rows where the combination of standard_type and standard_units is not in the config table.')
        print(f'Removed {rows_before_filter-rows_after_filter} rows. Remaining {rows_after_filter} rows.')
        return df

    def _calculate_active(self, row, cut):
        if np.isnan(row.standard_value):
            return row.comment_active
        elif row.active_direction == 1:  # Higher value is more active
            if row.standard_value >= cut:
                return 1
            else:
                return 0
        elif row.active_direction == -1:  # Lower value is more active
            if row.standard_value <= cut:
                return 1
            else:
                return 0
    
    def calculate_active(self,df):
        print(df.columns)     
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

        rows_to_drop = df[(df.comment_active.isnull()) &
        (df.activity_hc==0) & 
        (df.active_direction==1) & 
        (df.standard_relation.isin(['>','>=']))
        ].index
        df.drop(rows_to_drop, inplace=True)
        print(f'Removed {len(rows_to_drop)} cases with active direction +, relation ">", but labeled as not active')

        rows_to_drop = df[(df.comment_active.isnull()) &
        (df.activity_hc==0) & 
        (df.active_direction==-1) & 
        (df.standard_relation.isin(['<','<=']))
        ].index
        df.drop(rows_to_drop, inplace=True)
        print(f'Removed {len(rows_to_drop)} cases with active direction -, relation "<", but labeled as not active')

        rows_to_drop = df[(df.comment_active.isnull()) &
        (df.activity_hc==1) & 
        (df.active_direction==1) & 
        (df.standard_relation.isin(['<','<=']))
        ].index
        df.drop(rows_to_drop, inplace=True)
        print(f'Removed {len(rows_to_drop)} cases with active direction +, relation "<", but labeled as active')

        rows_to_drop = df[(df.comment_active.isnull()) &
        (df.activity_hc==1) & 
        (df.active_direction==-1) & 
        (df.standard_relation.isin(['>','>=']))
        ].index
        df.drop(rows_to_drop, inplace=True)
        print(f'Removed {len(rows_to_drop)} cases with active direction -, relation ">", but labeled as active')

        print('Cases remaining after filter:', len(df))
        return df
    
    def run(self, df):
        df = self.add_cutoffs(df)
        print(df.head())
        df = self.eliminate_no_direction(df)
        print(df.head())
        df = self.calculate_active(df)
        df = self.remove_inconsistencies(df)
        return df