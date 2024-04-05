from setuptools import setup, find_packages

setup(name='chemblbinarytasks',
      version='0.0.1',
      description='Tools for obtaining Chembl data for machine learning models',
      long_description=open('README.md').read().strip(),
      author="Ersilia Open Source Initiative",
      author_email="hello@ersilia.io",
      url='https://github.com/ersilia-os/chembl_ml_tools',
      license='GPLv3',
      python_requires='>=3.10',
      install_requires=[
        'pandas==2.2.1',
        'rdkit==2023.9.5',
        'psycopg2-binary==2.9.9'
      ],
      packages=find_packages(exclude=("utilities")),
      keywords='drug-discovery machine-learning ersilia chembl',
      classifiers=[
          "Programming Language :: Python :: 3.10",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
          "Topic :: Scientific/Engineering :: Artificial Intelligence",
      ],
      include_package_data=True
      )