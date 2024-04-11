# CONFIG Files

These files contains the descriptions and cutoffs used for the different assays

## activity_stds_lookup.csv

This table contains the descriptions of the "standard_type" variable, with their corresponding units of measurement.

Obtained from Chembl 33 by running this command in the postgresql command line:
```
\copy (SELECT * FROM activity_stds_lookup) to 'data/activity_stds_lookup.csv' with csv header;

```

## pathogens.csv
List of pathogens used to create the standard files. If other pathogens are used, consider re-doing the pipeline and adding the necessary fields on the `_manual.csv` files. See notebooks/units.ipynb for more information

## units.csv and ucum.csv
The units file contains all the units that appear in the datasets obtained from ChEMBL_33 using the pathogens specified in pathogens.csv. These units have been processed in the [UCUM server](https://ucum.nlm.nih.gov/ucum-lhc/demo.html) to obtain the Unified Code for Units of Measure.
The ucum.csv file, correspondingly, contains the units and its standardisation to ucum units, as well as the necessary conversions to final units, like umol/L. If a unit does not appear in the file, it will not be processed. The units are strings (case sensitive)

## st_type_summary.csv and st_type_summary_manual.csv
A count of how many standard_type - final_unit combinations exist in the data collected from ChEMBL33 (using the pathogens specified in pathogens.csv). Only combinations with more than 250 instances are considered, and a manual curation has been performed to decide which assays will be used in this pipeline. Simply change the `use` column from 0 to 1 to include a standard_type - final_unit combination, and define its cutoff in the following files. we also include information about the direction of activity (-1 for assays where smaller values are preferred and 1 for assays where larger values are preferred)

## percentiles and cutoffs_manual.csv
The percentiles file contain the percentile distribution of the standard_type - final_unit combinations identified in `st_type_summary_manual.csv`, and whether the cut-off has been defined by expertise or by the percentiles. Final cut-offs can be found on the `cutoffs_manual.csv` file under Low_cut (for a non-stringent activity cut-off) and High_cut (for a more restrictive activity cut-off)