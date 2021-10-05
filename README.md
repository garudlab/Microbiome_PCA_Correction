# Code to implement analyses of PCA correction and other corrections in Microbiome data 

Preprint of the paper can be found on [biorxiv](https://www.biorxiv.org/content/10.1101/2021.03.19.436199v1)


## Table of contents
1. [Quick summary of steps](#quick)
2. [Breakdown of steps and expected output](#breakdown)



## <a name =quick> Quick summary of steps </a>

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
Rscript Correction/calc_pcs.R <your_directory> <input format of feature table: txt or rds> <data transformation (relative abundance or CLR): rel or rel_clr>

Rscript Correction/calc_pcs.R ~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 rds rel_clr
```

### Step 3: Correction

**Input**: feature table & metadata (for limma, percentilenorm, combat, etc.), & PCA scores (for PCA correction) |
**Output**: eigenvalues and pca scores

*make sure metadata categories are already cleaned up to distinct groups if you're doing to use percentilenorm


```
Rscript Correction/correction.R ~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 rds rel limma study bin_crc_normal
```

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

### Step 2: Calculate PCs (if doing PCA correction)

```
Loading required package: bigstatsr
Time difference of 1.820936 secs
[1] "eigen values"
 [1] 1.58260062 0.21910492 0.12427020 0.09583746 0.08502080 0.07795310
 [7] 0.06726517 0.06550892 0.06077923 0.05249138 0.05210248 0.04988210
[13] 0.04238802 0.03892835 0.03716319
[1] "dim pca scores"
[1] 723  15
```