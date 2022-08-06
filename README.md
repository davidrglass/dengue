# dengue
Code to reproduce figures from preprint:

Batch-corrected fcs files are found here:

Instructions for each script written at the top of the file.

Scripts need to be run in the following order:
1. fcs_to_csv.R
  Requires: finalized_fcs.csv
2. cluster.R
3. integrate_and_subsample.R
  Requires: metadata.csv
4. figures.R
  Requires: severeonly.csv; exposure.csv

Note that batch_normalize.R does not need to be run, as fcs files are already batch corrected.
Included for transparency.
