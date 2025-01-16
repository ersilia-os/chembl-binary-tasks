import psycopg2
import pandas as pd
import pandas.io.sql as sqlio
import os.path
from rdkit import Chem


class ChemblMoleculeSampler:
    """Get samples of molecules from Chembl database"""
    
    # Name of the file that contains the full list of molecules downloaded from Chembl
    FILENAME_ALL_MOLECULES = 'chembl_all_molecules.csv'

    def __init__(self, data_path, db_user, db_password,
                 db_name='chembl_35', db_host='localhost', db_port=5432,
                 max_heavy_atoms=None):
        """
        Load the list of all molecules into dataframe df_all_molecules
        
        Parameters:
            data_path (str): A valid path where we store the csv file with the full list of molecules.
            db_user, db_password, db_name, db_host, db_port: Parameters to connect to the Chembl database
            max_heavy_atoms (int, optional): Select only molecules with at most this
                number of non-hydrogen atoms.
        """
        
        self.data_path = data_path
        
        # Check that data_path exists
        if (not os.path.isdir(data_path)):
            raise ValueError('The parameter data_path must contain a directory that exists. '
                            f'Directory {data_path} does not exist.')
            
        # Check if the list of all molecules is already downloaded
        file_with_path = os.path.join(self.data_path, self.FILENAME_ALL_MOLECULES)
        if not os.path.exists(file_with_path):
            df = self.download_all_molecules(db_user, db_password, db_name, db_host, db_port)
            df.to_csv(file_with_path, index=False)
            print(f'List of molecules downloaded into file {self.FILENAME_ALL_MOLECULES}. Table shape (rows,columns): {df.shape}')

        # Load the list of all molecules from csv into a dataframe
        df_all_molecules1 = pd.read_csv(file_with_path, index_col=0)

        if max_heavy_atoms is not None:
            # For performance reasons, we limit the number of molecules
            MAX_MOLECULES_TO_COUNT_ATOMS = 100000
            if len(df_all_molecules1) > MAX_MOLECULES_TO_COUNT_ATOMS:
                df_all_molecules1 = df_all_molecules1.sample(n=MAX_MOLECULES_TO_COUNT_ATOMS)

            # Filter out molecules larger than max_heavy_atoms
            molecules = df_all_molecules1.canonical_smiles.apply(Chem.MolFromSmiles)
            heavy_atoms = molecules.apply(lambda x: self.GetNumHeavyAtoms_prevent_exception(x))
            self.df_all_molecules = df_all_molecules1[heavy_atoms <= max_heavy_atoms]
            print(f'{len(self.df_all_molecules)} candidate molecules available after'
                  f' filtering for NumHeavyAtoms <= {max_heavy_atoms}')
        else:
            self.df_all_molecules = df_all_molecules1
            print(f'{len(self.df_all_molecules)} candidate molecules available')

    def GetNumHeavyAtoms_prevent_exception(self, x):
        """Encapsulate method GetNumHeavyAtoms to return 1e10 in case of exception

        Background: Occasionally the method GetNumHeavyAtoms fails for specific
        molecules with the error 'Can't kekulize mol'.
        In this case we do not want the program to fail, just report a very
        large value for the molecule so it will not pass the max_heavy_atoms filter.
        """
        try:
            return x.GetNumHeavyAtoms()
        except:
            return 1e10


    def download_all_molecules(self, db_user, db_password, db_name, db_host, db_port):
        """
        Download all molecules from the Chembl database and store them in a csv file
        
        Parameters:
            db_user, db_password, db_name, db_host, db_port: Parameters to connect to the Chembl database
            
        Returns:
            df (DataFrame): Dataframe with molecule data
        """
        
        # Connect to chembl database
        conn = psycopg2.connect(database=db_name, user=db_user, password=db_password, host=db_host, port=db_port)
        cursor = conn.cursor()
        
        print('Downloading all molecules from the Chembl database')
        # Run query on Chembl database: Get all compounds that have at least one activity.
        cursor.execute(
"""
CREATE TEMPORARY TABLE tmp_all_compounds AS
SELECT DISTINCT
  md.chembl_id,
  cnd_s.canonical_smiles
FROM activities act
  JOIN molecule_dictionary md ON act.molregno = md.molregno
  JOIN compound_structures cnd_s ON md.molregno = cnd_s.molregno
ORDER BY chembl_id
"""
        )
        
        # Count rows
        cursor.execute("SELECT count(*) FROM tmp_all_compounds")
        print('Count of rows extracted:', cursor.fetchone()[0])

        # Download into dataframe
        sql = "SELECT * FROM tmp_all_compounds"
        df = sqlio.read_sql_query(sql, conn)
        print('The list of molecules has been obtained from the Chembl database.')

        # Close database connection
        conn.close()
        
        return df

    def negative_sample(self, num_molecules, list_positive_molecules=None):
        """
        Get a random sample of molecules to be used as negative cases for a model.
        If a data frame with known positive cases is provided, these will not be included in the sample.
        
        Parameters:
            num_molecules (int): Number of molecules required in the sample
            list_positive_molecules (list, optional): Positive cases, not to be included in the sample.
                The list must contain chembl_id identifiers (example: 'CHEMBL1430368')
            
        Returns:
            df_negatives_sample (DataFrame) Dataframe with sample of negative cases
        """
        
        print(f'{len(self.df_all_molecules)} total molecules available')

        # If dataframe with positive cases is provided, exclude those cases previous to the sampling
        if list_positive_molecules:
            df_negatives = self.df_all_molecules[~self.df_all_molecules.index.isin(list_positive_molecules)]
            print(f'{len(df_negatives)} negative molecules available'
                f' (after substracting the provided positive cases)')
        else:
            df_negatives = self.df_all_molecules

        num_negatives = len(df_negatives)
        if num_molecules > num_negatives:  # Fewer available than requested: provide all available
            print(f'Warning: The number of requested molecules ({num_molecules}) is higher than the number '
                  f'of available negative molecules ({num_negatives}). A dataframe of size {num_negatives} '
                  f'will be returned.'
                  )
            return df_negatives.reset_index()  # Convert index into column
        else:  # More available than requested: Random sample
            df_negatives_sample = df_negatives.sample(n=num_molecules)
            print(f'{num_molecules} molecules provided')
            return df_negatives_sample.reset_index()  # Convert index into column
        
