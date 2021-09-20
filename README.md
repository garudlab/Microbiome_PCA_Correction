# Microbiome_PCA_Correction


### Step 1: CLR Transformation

Correction > transformations.R

```
Rscript transformations.R ~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 rds
```
Example output:

```
[1] "Starting CLR transformation of data"
[1] "Adding pseudocount to zeroes"
[1] 6.748107e-06
[1] "Number of zero values in data"
[1] 0
[1] "Proportion of values in data that are zero"
[1] 0
[1] "Transformation done and exported to /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5"
```

### Step 2: Calculate PCs (if doing PCA correction)

Correction > calc_pcs.R

### Step 3: Correction

Correction > correction.R

### Analysis Opt 1: Correlation analyses

Correlation > pc_correlations_plot.R

### Analysis Opt 2: Titration analyses

Titration > titratoins.ipynb

### Analysis Opt 3: Variance partitioning

VariancePartitioning > variance_partitioning.R

VariancePartitioning > variance_partitioning_plotting.R

### Analysis Opt 4: Prediction

For binary phenotype:

RF_Classifier > classifier.py

For continuous phenotype:

LinRegression > regression.py

