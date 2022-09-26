# dengue
Code to reproduce figures from "Magnitude and kinetics of the human immune cell response associated with severe dengue progression by single-cell proteomics" https://www.biorxiv.org/content/10.1101/2022.09.21.508901v1.full

Batch-corrected fcs files are found here: https://flowrepository.org/id/FR-FCM-Z5MQ

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
