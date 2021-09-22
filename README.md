# Code to implement analyses of PCA correction and other corrections in Microbiome data 

Preprint of the paper can be found on [biorxiv](https://www.biorxiv.org/content/10.1101/2021.03.19.436199v1)


## Table of contents
1. [Quick summary of steps](#quick)
2. [Breakdown of steps and expected output](#breakdown)



##<a name =quick> Quick summary of steps </a>

### Needed input data:

1. **metadata.txt** that contains the information regarding phenotype and other experimental variables (if performing correlation analyses)

2. **feature\_table\_rel.rds** or **feature_table\_rel.txt**

These files should be stored in a directory of choice. The path for this directory will be hereafter referred to `<your_directory>`


### Step 1: CLR Transformation

**Input**: feature table |
**Output**: CLR-transformed data in both rds and txt formats 

```
Rscript Correction/transformations.R <your_directory> <input format of feature table: txt or rds>

Rscript Correction/transformations.R ~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 rds
```

### Step 2: Calculate PCs (if doing PCA correction)

**Input**: feature table |
**Output**: eigenvalues and pca scores

```
Rscript Correction/calc_pcs.R <your_directory> <input format of feature table: txt or rds>

Rscript Correction/calc_pcs.R ~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 rds
```

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


##<a name=breakdown> Breakdown of steps and expected output</a>




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
