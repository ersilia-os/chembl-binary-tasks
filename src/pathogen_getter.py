import os
import sys
import pandas as pd
import psycopg2
import pandas.io.sql as sqlio
from rdkit import Chem

abspath = os.path.abspath(__file__)
sys.path.append(abspath)
from default import CHEMBL_PWD, CHEMBL_USR, PATHOGENSPATH, DATABASE_NAME
import argparse


class PathogenGetter():
    """
    This class obtains the full data related to the specified pathogens and creates the directories to save it
    """
    def __init__(self, patho_search, strict_match):
        self.db_user = CHEMBL_USR
        self.db_pwd = CHEMBL_PWD
        self.patho_search = patho_search
        self.strict_match = strict_match
        self._get_taxids()

    def _get_taxids(self):
        # Connect to chembl database
        conn = psycopg2.connect(database=DATABASE_NAME, user=self.db_user, password=self.db_pwd, host='localhost', port=5432)
        cursor = conn.cursor()
        # Run query
        if self.strict_match:
            cursor.execute(f"SELECT DISTINCT tax_id FROM target_dictionary WHERE organism ilike '{self.patho_search}'")
        else:
            cursor.execute(f"SELECT DISTINCT tax_id FROM target_dictionary WHERE organism ilike '%{self.patho_search}%'")
        taxids = cursor.fetchall()
        print(f'Found {len(taxids)} taxids for {self.patho_search}')
        self.taxids = [taxid[0] for taxid in taxids]

    def _get_org_str(self):
        if self.strict_match:
            print(f'Using strict match for pathogen search: {self.patho_search}')
            return f'{self.patho_search}'
        else:
            print(f'Using non-strict match for pathogen search: {self.patho_search}')
            return f'%{self.patho_search}%'
        
    def _get_taxid_str(self):
        return ','.join([str(taxid) for taxid in self.taxids])

    def chembl_activity_target(self, max_heavy_atoms=None,db_name=DATABASE_NAME, db_host='localhost', db_port=5432):
        # Connect to chembl database
        conn = psycopg2.connect(database=db_name, user=self.db_user, password=self.db_pwd, host=db_host, port=db_port)
        cursor = conn.cursor()
        # Run query first step
        query_string = """
    CREATE TEMPORARY TABLE tmp_activity_target_1 AS
    WITH
    target_details AS (
        SELECT
        td.tid,
        string_agg(cnt_s.accession || ' ' || pc.pref_name, ';') AS protein_accession_class
        FROM target_dictionary td
        LEFT JOIN target_components tc ON td.tid = tc.tid
        LEFT JOIN component_sequences cnt_s ON tc.component_id = cnt_s.component_id
        LEFT JOIN component_class cnt_c ON cnt_s.component_id = cnt_c.component_id
        LEFT JOIN protein_classification pc ON cnt_c.protein_class_id = pc.protein_class_id
        WHERE td.tax_id IN ({0})
        GROUP BY td.tid
    )
    SELECT
    docs.doc_id,
    a.assay_id,
    a.assay_type,
    a.confidence_score AS assay_confidence_score,
    a.bao_format AS assay_bao_format,
    act.activity_id,
    md.chembl_id AS compound_chembl_id,
    md.molregno,
    cnd_s.canonical_smiles,
    act.standard_type,
    act.standard_value,
    act.standard_units,
    act.standard_relation,
    act.pchembl_value,
    act.activity_comment,
    td.chembl_id AS target_chembl_id,
    td.target_type,
    td.organism AS target_organism,
    td.pref_name AS target_pref_name,
    td.tax_id AS target_tax_id,
    prot_d.protein_accession_class,
    docs.year,
    docs.pubmed_id,
    a.chembl_id as assay_chembl_id,
    a.description AS assay_description
    FROM target_dictionary td
    JOIN assays a ON td.tid = a.tid
    JOIN docs ON a.doc_id = docs.doc_id
    JOIN activities act ON a.assay_id = act.assay_id
    JOIN molecule_dictionary md ON act.molregno = md.molregno
    JOIN compound_structures cnd_s ON md.molregno = cnd_s.molregno
    LEFT JOIN target_details prot_d ON td.tid = prot_d.tid
    WHERE
       act.standard_value IS NOT NULL AND
       td.tax_id IN ({0})
    """.format(self._get_taxid_str())
        print(query_string)
        cursor.execute(query_string)

        cursor.execute("SELECT * FROM tmp_activity_target_1")
        
        # Run query second step (consolidate duplicates)
        cursor.execute(
    """
    CREATE TEMPORARY TABLE tmp_activity_target AS
    SELECT
    min(doc_id) as doc_id,
    min(assay_id) as assay_id,
    min(activity_id) as activity_id,
    assay_type,
    max(assay_confidence_score) as assay_confidence_score,
    assay_bao_format,
    compound_chembl_id,
    molregno,
    max(canonical_smiles) AS canonical_smiles,
    standard_type,
    standard_value,
    standard_units,
    standard_relation,
    pchembl_value,
    string_agg(cast(activity_comment as varchar), ';') AS activity_comment,
    target_chembl_id,
    max(target_type) AS target_type,
    max(target_organism) AS target_organism,
    max(target_pref_name) AS target_pref_name,
    max(target_tax_id) AS target_tax_id,
    max(protein_accession_class) AS protein_accession_class,
    max(year) AS year,
    max(pubmed_id) AS pubmed_id,
    max(assay_chembl_id) AS assay_chembl_id,
    count(*) as count_activity_rows,
    string_agg(cast(doc_id as varchar), ';') AS doc_id_all,
    string_agg(cast(assay_id as varchar), ';') AS assay_id_all,
    string_agg(cast(activity_id as varchar), ';') AS activity_id_all,
    max(assay_description) AS assay_description
    FROM tmp_activity_target_1
    GROUP BY
    target_chembl_id,
    compound_chembl_id,
    molregno,
    assay_type,
    assay_bao_format,
    standard_type,
    standard_value,
    standard_units,
    standard_relation,
    pchembl_value
    """
        )

        cursor.execute(
    """
    CREATE TEMPORARY TABLE tmp_activity_target_parent_1 AS
    SELECT 
    tat.*,
    mh.active_molregno AS parent_molregno
    FROM 
        tmp_activity_target tat
    LEFT JOIN 
        molecule_hierarchy mh
    ON 
        tat.molregno = mh.molregno;
    """
        )

        cursor.execute(
    """
    CREATE TEMPORARY TABLE tmp_activity_target_parent AS
    SELECT
    tatp1.*,
    md.canonical_smiles AS parent_canonical_smiles
    FROM
        tmp_activity_target_parent_1 tatp1
    LEFT JOIN
        compound_structures md
    ON
        tatp1.parent_molregno = md.molregno;
    """
        )

        # Count rows
        cursor.execute("SELECT count(*) FROM tmp_activity_target_parent")
        print(f'{cursor.fetchone()[0]} rows extracted from ChEMBL:')
        
        # NOTE: Selected rows are limited in case there are too many
        # Pending to give a warning if not all rows are shown
        sql = "SELECT * FROM tmp_activity_target_parent limit 10000000"
        df = sqlio.read_sql_query(sql, conn)
        print(df.head())
        conn.close()  # Close database connection

        if max_heavy_atoms is not None:
            print(f'{len(df)} cases selected before filtering for NumHeavyAtoms <= {max_heavy_atoms}')
            # Filter out molecules larger than max_heavy_atoms
            molecules = df.canonical_smiles.apply(Chem.MolFromSmiles)
            heavy_atoms = molecules.apply(lambda x: x.GetNumHeavyAtoms())
            df = df[heavy_atoms <= max_heavy_atoms]

        print(f'{len(df)} cases in resulting dataset.', )
        print('Selected organisms:')
        print(df.target_organism.value_counts())
        print('Calculating InChIKeys...')
        inchikeys = []
        parent_inchikeys = []
        for k,v in df.iterrows():
            mol = Chem.MolFromSmiles(v.canonical_smiles)
            inchikey = Chem.inchi.MolToInchiKey(mol)
            inchikeys.append(inchikey)
            mol = Chem.MolFromSmiles(v.parent_canonical_smiles)
            inchikey = Chem.inchi.MolToInchiKey(mol)
            parent_inchikeys.append(inchikey)
        df['inchikey'] = inchikeys
        df['parent_inchikey'] = parent_inchikeys

        columns_a = [
            "activity_id_all",
            "assay_id_all",
            "assay_chembl_id",
            "assay_type",
            "standard_relation",
            "standard_value",
            "standard_units",
            "pchembl_value",
        ]

        columns_b = [
            "target_chembl_id",
            "target_type",
            "target_organism",
            "target_pref_name",
            "protein_accession_class",
        ]

        columns_c = [
            "compound_chembl_id",
            "inchikey",
            "canonical_smiles",
            "parent_inchikey",
            "parent_canonical_smiles"
        ]

        columns_d = [
            "year"
        ]

        return df.copy()
    
if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description='Get pathogen data from ChEMBL database.')
    parser.add_argument('-d', '--output_dir', type=str, required=True, help='Output directory')
    parser.add_argument('-p', '--pathogen_code', type=str, required=False, help='Pathogen code to search for')
    parser.add_argument('-s', '--strict_match', action='store_true', help='Use strict match for pathogen search')
    parser.add_argument('--max_heavy_atoms', type=int, default=100, help='Maximum number of heavy atoms for molecules')
    args = parser.parse_args()
    # Obtain desired pathogen data  
    df_pathogens = pd.read_csv(PATHOGENSPATH)
    list_pathogen_codes = list(df_pathogens.pathogen_code)
    list_pathogen_search_text = list(df_pathogens.search_text)
    print(list_pathogen_codes)
    if args.pathogen_code is not None:
        if args.pathogen_code not in list_pathogen_codes:
            raise ValueError(f'Pathogen code {args.pathogen_code} not found in pathogens.csv')
    for i, patho_code in enumerate(list_pathogen_codes):
        if args.pathogen_code is not None:
            if patho_code != args.pathogen_code:
                continue
        print('------------------------------------------------------------')
        print(f'Creating data for pathogen {patho_code} ({list_pathogen_search_text[i]})')
        print('------------------------------------------------------------')
        patho_getter = PathogenGetter(list_pathogen_search_text[i], strict_match=args.strict_match)
        df = patho_getter.chembl_activity_target(max_heavy_atoms=args.max_heavy_atoms)

        if not os.path.exists(os.path.join(args.output_dir, patho_code)):
            os.makedirs(os.path.join(args.output_dir, patho_code))
        df.to_csv(os.path.join(args.output_dir, patho_code, "00_{}_original.csv".format(patho_code)), index=False)