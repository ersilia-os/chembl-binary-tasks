# CONFIG Files

These files contains the descriptions and cutoffs used for the different assays

## Table activity_stds_lookup.csv

This table contains the descriptions of the "standard_type" variable, with their corresponding units of measurement.

Obtained from Chembl 31 by running this command in the postgresql command line:
```
\copy (SELECT * FROM activity_stds_lookup) to 'data/activity_stds_lookup.csv' with csv header;
```
