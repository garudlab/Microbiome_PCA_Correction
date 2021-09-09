## New Revision sequence

Thomasr_complete_otu:
0.87029 % zeroes. Pseudocount = 6.632665e-08
bmc combat percentilenorm limma DCC


### All thomas jobs

Rscript pc_correlations.R Thomasr_max_k6

```
for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 16 -t 24 -hp -v 3.6.0 -arg Thomasr_max_k"$k"; done
 
for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 12 -t 24 -hp -v 3.6.0 -arg Thomasr_max_k"$k" -arg rel; done
 
 
/u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 12 -t 24 -hp -v 3.6.0 -arg Thomasr_complete_otu -arg rel_clr
 
 Rscript pc_correlations.R Thomasr_complete_otu logCPM

 /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations_plot.R -m 4 -t 24 -hp -v 3.6.0 -arg Thomasr_complete_otu -arg rel_clr_scale -arg dataset_name
 
 Rscript pc_correlations_plot.R Thomasr_complete_otu rel_clr_scale dataset_name

for k in 6; do for p in 1 2 3 4; do for r in rel_clr; do for method in pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -hp -t 24 -v 3.6.0 -arg Thomasr_max_k"$k" -arg $r -arg $p -arg $method; done; done; done; done

for p in 2; do for r in vsd; do for method in pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 18 -hp -t 24 -v 3.6.0 -arg Thomasr_complete_otu -arg $r -arg $p -arg $method; done; done; done;

Rscript correction.R Thomasr_complete_otu vsd 1 pca


for k in 7; do for p in 1; do for method in bmc combat percentilenorm limma DCC pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 18 -hp -t 24 -v 3.6.0 -arg Thomasr_max_k"$k" -arg rel -arg $p -arg $method; done; done; done


 

for k in 1 2 3 4 5; do for method in pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 12 -hp -t 24 -v 3.6.0 -arg Thomasr_max_k7 -arg rel_clr -arg $k -arg $method; done; done

# 180 210 240 270 300 330 360)
indexs=(30 60 90 120 150)  
methods=(bmc combat percentilenorm limma DCC)

indexs=(30 60 90 120 150)
methods=(pca1counts pca2counts pca3counts pca4counts pca5counts)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 500));
	method="${methods[$i]}"
	./qsub_classifier.sh $index Thomasr_complete_otu vsd $method 1 bin_crc_normal
done;


indexs=(30 60 90 120 150 180)
methods=(clr_pca1roundcounts clr_pca1 clr_pca2roundcounts clr_pca2 clr_pca3roundcounts clr_pca3)

 
indexs=(30 60 90 120 150 180 210 240 270 300 330)
methods=(bmc combat percentilenorm limma DCC nocorrection clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts)

indexs=(30 60 90 120 150 180)
methods=(clr clr_pca1 clr_pca2 clr_pca3 clr_pca4 clr_pca5)

indexs=(30 60)
methods=(clr clr_pca1)
rel=
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 1200));
	method="${methods[$i]}"
	./qsub_classifier.sh $index Thomasr_complete_otu $rel $method 1 bin_crc_normal
done;


indexs=(30 60)
methods=(logCPM vsd)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 1200));
	method="${methods[$i]}"
	./qsub_classifier.sh $index Thomasr_complete_otu $method nocorrection 1 bin_crc_normal
done;




```
   8963886 0.50500 RF         briscoel     qw    08/02/2021 08:59:26                                                                   1 1274-1286:1
   8963890 0.50500 RF         briscoel     qw    08/02/2021 08:59:27                                                                   1 1290-1316:1
   8963892 0.50500 RF         briscoel     qw    08/02/2021 08:59:28                                                                   1 1320-1346:1
   8963895 0.50500 RF         briscoel     qw    08/02/2021 08:59:29                                                                   1 1350-1376:1
   8963897 0.00000 RF         briscoel     qw    08/02/2021 08:59:30                                                                   1 1380-1406:1
   
```


metric=val_auc
for method in clr clr_pca1 clr_pca2 clr_pca3 clr_pca4 clr_pca5; do python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction $method --lodo 0 --phenotype bin_crc_normal --metric $metric; done


metric=val_auc
for method in pca1counts pca2counts pca3counts pca4counts pca5counts; do python process_rf_result.py --folder Thomasr_complete_otu --trans vsd --correction $method --lodo 1 --phenotype bin_crc_normal --metric $metric; done



```

## Process

```

for p in 1 2 3; do python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction clr_pca"$p" --lodo 1 --phenotype bin_crc_normal; done


for p in 1 2 3; do python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction clr_pca"$p"roundcounts --lodo 1 --phenotype bin_crc_normal; done

for method in nocorrection; do python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction $method --lodo 0 --phenotype bin_crc_normal; done

for method in nocorrection bmc combat limma DCC clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts; do python process_rf_result.py --folder Thomasr_max_k7 --trans rel --correction $method --lodo 0 --phenotype bin_crc_normal --metric val_auc; done


```


### AGP jobs

```
 /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k8
 
/u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k8
 
 /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu -arg rel -arg 2 -arg pca
 

for method in DCC; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu -arg rel -arg 2 -arg $method; done

for k in 1 3 4 5; do for trans in rel_clr rel_clr_scale; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 24 -t 24 -v 3.6.0 -arg AGPr_complete_otu -arg $trans -arg $k -arg pca; done; done

```


## All AGP kmer jobs
```
for data in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data; done

for data in 6; do for trans in rel rel_clr; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 12 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu -arg $trans; done; done

for data in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data; done


20 for k7, 15 for k6, 13 for k5


=====
for data in 8; do for method in combat limma bmc DCC; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data -arg rel -arg 2 -arg $method; done; done

for data in 8; do for p in 1 2 3 4 5; do for trans in rel_clr; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -hp -t 24 -v 3.6.0 -arg AGPr_max_k"$data" -arg $trans -arg $p -arg pca; done; done; done


for data in 5; do for p in 1; do for trans in rel_clr; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 12 -hp -t 24 -v 3.6.0 -arg AGPr_max_k"$data" -arg $trans -arg $p -arg pca; done; done; done







=====

for data in 6; do for method in combat; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 18 -hp -t 24 -v 3.6.0 -arg AGPr_max_k$data -arg rel -arg 2 -arg $method; done; done

for k in 5; do for p in 5; do for method in pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 24 -t 24 -v 3.6.0 -arg AGPr_max_k"$k" -arg rel_clr -arg $p -arg $method; done; done; done

```

=====

```
python regression.py --folder AGPr_complete_otu --trans rel --correction nocorrection --lodo 0 --phenotype bmi_corrected 



indexs=(30 60 90 120 150 180 210 240 270 300 330 360 390 420 450)
methods=(nocorrection bmc limma DCC clr clr_pca1 clr_pca2 clr_pca3 clr_pca4 clr_pca5 clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 1));
	method="${methods[$i]}";
	./qsub_classifier.sh $index AGPr_max_k8 rel $method 0 bin_antibiotic_last_year;
done;

indexs=(90)
methods=(clr)

indexs=(120 150 180 210 240 270 300)
methods=(DCC clr clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 1));
	method="${methods[$i]}";
	./qsub_classifier.sh $index AGPr_max_k8 rel $method 0 bin_antibiotic_last_year;
done;


indexs=(90)
methods=(clr)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 559));
	method="${methods[$i]}";
	./qsub_classifier.sh $index AGPr_max_k7 rel $method 0 bin_antibiotic_last_year;
done;

indexs=(60 90 120 150 180)
methods=(clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts)

indexs=(30 60 90 120 150 180 210 240 270)
methods=(nocorrection bmc limma DCC clr clr_pca1 clr_pca2 clr_pca3 clr_pca4 clr_pca5)


metric=val_auc
for method in nocorrection bmc combat limma DCC clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts; do python process_rf_result.py --folder AGPr_max_k6 --trans rel --correction $method --lodo 0 --phenotype bin_antibiotic_last_year --metric $metric; done


metric=val_auc
for method in method in nocorrection bmc combat limma DCC clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts; do python process_rf_result.py --folder AGPr_max_k8 --trans rel --correction $method --lodo 0 --phenotype bin_antibiotic_last_year --metric $metric; done




```

_italics_

**bold**

[link](www.google.com)

__double__

<a> Title </a>




# Gibbons jobs
bmc limma combat DCC
```
for k in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 16 -t 24 -hp -v 3.6.0 -arg Gibbonsr_max_k"$k"; done

for k in 6; do for trans in rel rel_clr; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 12 -t 24 -hp -v 3.6.0 -arg Gibbonsr_complete_otu -arg $trans;  done; done

for k in 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 16 -t 24 -hp -v 3.6.0 -arg Gibbonsr_max_k"$k";  done

 ## all methods
  


Rscript correction.R Gibbonsr_max_k5 rel 2 combat

for k in 6; do for method in limma; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 12 -t 24 -v 3.6.0 -arg Gibbonsr_max_k"$k" -arg rel -arg $k -arg $method; done; done

for k in 8; do for method in combat; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -t 24 -v 3.6.0 -arg Gibbonsr_max_k"$k" -arg rel -arg 2 -arg $method; done; done



# pc methods

for k in 7; do for p in 4 5; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -t 24 -v 3.6.0 -arg Gibbonsr_max_k"$k" -arg rel_clr -arg $p -arg pca; done; done


for k in 8; for method in bmc; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -t 24 -v 3.6.0 -arg Gibbonsr_max_k"$k" -arg rel -arg 2 -arg pca; done; done

 
 Rscript correction.R Gibbonsr_complete_otu rel 2 DCC
 
 Rscript pc_correlations.R Gibbonsr_complete_otu 
 
 
indexs=(30 60 90 120 150 180 210 240 270 300 330)
methods=(bmc combat percentilenorm limma DCC nocorrection clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 1));
	method="${methods[$i]}"
	./qsub_classifier.sh $index Gibbonsr_max_k8 $method 0 bin_crc_normal
done;



indexs=(130 160 190 220 250 280 )
methods=(clr clr_pca1 clr_pca2 clr_pca3 clr_pca4 clr_pca5)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 400));
	method="${methods[$i]}"
	./qsub_classifier.sh $index Gibbonsr_max_k8 rel $method 1 bin_crc_normal
done;


metric=val_auc
for method in nocorrection bmc combat limma DCC clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts; do python process_rf_result.py --folder Gibbonsr_max_k8 --trans rel --correction $method --lodo 0 --phenotype bin_crc_normal --metric $metric; done


metric=val_auc
for method in clr clr_pca1 clr_pca2 clr_pca3 clr_pca4 clr_pca5; do python process_rf_result.py --folder Gibbonsr_max_k8 --trans rel --correction $method --lodo 0 --phenotype bin_crc_normal --metric $metric; done

metric=val_auc
for method in clr; do python process_rf_result.py --folder Gibbonsr_complete_otu --trans rel --correction $method --lodo 1 --phenotype bin_crc_normal --metric $metric; done




 ```
 
 
 needed:
 all Pca methods for AGP_otu
 
 all mehods for AGP_max_k5
 all metiods for AGP+max_k6
  all metiods for AGP+max_k7
  
  remmebr lodo
  
  
indexs=(30 60 90 120 150 180 210 240 270 300 330)
methods=(nocorrection bmc combat percentilenorm limma DCC clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts)


indexs=(30 60 90 120 150 180 210 240 270 300 330)
methods=(nocorrection bmc combat percentilenorm limma DCC clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts clr clr_pca1 clr_pca2 clr_pca3 clr_pca4 clr_pca5)

indexs=(30 60 90 120 150 180 210 240 270 300 330)
methods=(nocorrection bmc combat percentilenorm limma DCC clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 1));
	method="${methods[$i]}"
	./qsub_classifier.sh $index AGPr_complete_otu $method 0 bin_antibiotic_last_year
done;


## KAPLAN

```
for k in 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 20 -t 24 -hp -v 3.6.0 -arg Kaplanr_max_k"$k"; done

/u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 20 -t 24 -hp -v 3.6.0 -arg Kaplanr_complete_otu

for k in 6; do for trans in rel rel_clr; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 12 -t 24 -hp -v 3.6.0 -arg Kaplanr_complete_otu -arg $trans; done; done

/u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 16 -t 24 -hp -v 3.6.0 -arg Kaplanr_complete_otu

Rscript pc_correlations.R Kaplanr_complete_otu

Rscript correction.R Kaplanr_max_k5 rel 2 combat
NOTE:8698987

for method in DCC limma combat bmc; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -hp -m 8 -t 24 -v 3.6.0 -arg Kaplanr_max_k5 -arg rel -arg 2 -arg $method; done
NOTE: 8699247 - 50




for k in 6; do for method in combat; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -hp -t 24 -v 3.6.0 -arg Kaplanr_max_k"$k" -arg rel -arg 2 -arg $method; done; done

for k in 8; do for method in DCC limma combat bmc; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg Kaplanr_max_k"$k" -arg rel -arg 2 -arg $method; done; done


------
for p in 1 2 3 4 5; do for trans in rel_clr_scale; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -hp -t 24 -v 3.6.0 -arg Kaplanr_complete_otu -arg $trans -arg $p -arg pca; done; done;
NOTE:8699259 to 63


for k in 5; do for p in 1 2 3 4 5; do for trans in rel_clr_scale; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -hp -t 24 -v 3.6.0 -arg Kaplanr_max_k"$k" -arg $trans -arg $p -arg pca; done; done; done

for k in 8; do for p in 1 2 3; do for trans in rel_clr; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -hp -t 24 -v 3.6.0 -arg Kaplanr_max_k"$k" -arg $trans -arg $p -arg pca; done; done; done



NOTE: two anove: 8698969 to 8698985

```


### prediction
```
diabetes_self_v2

indexs=(30)
methods=(clr_pca3counts)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 2000));
	method="${methods[$i]}"
	./qsub_classifier.sh $index Kaplanr_max_k8 $method 0 diabetes_self_v2
done;

indexs=(30 60 90 120 150 180 210 240 270 300 330)
methods=(bmc combat percentilenorm limma DCC nocorrection clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts)



indexs=(30 60 90 120 150 180 210 240 270 300 330 360 390 420 450 480 510 540)
methods=(bmc combat percentilenorm limma DCC nocorrection clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts clr_scale clr clr_pca1 clr_pca2 clr_pca3 clr_pca4 clr_pca5)


indexs=(400 430 460 490)
methods=(bmc combat limma DCC)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 2500));
	method="${methods[$i]}"
	./qsub_regression.sh $index Kaplanr_max_k8 $method 0 bmi_v2
done;
```

### process prediction

metric=val_auc
for method in nocorrection bmc combat limma DCC clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts; do python process_rf_result.py --folder Thomasr_max_k6 --trans rel --correction $method --lodo 0 --phenotype bin_crc_normal --metric $metric; done

metric=val_auc
for method in nocorrection bmc combat limma DCC clr_pca1counts clr_pca2counts clr_pca3counts clr_pca4counts clr_pca5counts; do python process_rf_result.py --folder Kaplanr_max_k6 --trans rel --correction $method --lodo 1 --phenotype diabetes_self_v2 --metric $metric; done


 


python regression.py --folder AGPr_max_k5 --trans rel --correction nocorrection --lodo 0 --phenotype bmi_corrected 




# R WORK

for k in 8; do scp -r hoffy:/u/home/b/briscoel/project-halperin/MicroBatch/data/Kaplanr_max_k"$k"/*PRED* Kaplanr_max_k"$k"/; cp Kaplanr_max_k"$k"/PRED_OUTPUT_rel_clr_pca3counts_lodo_False.csv Kaplanr_max_k"$k"/PRED_OUTPUT_rel_clr_pca33counts_lodo_False.csv; done


for k in 5 6 7; do scp -r hoffy:/u/home/b/briscoel/project-halperin/MicroBatch/data/AGPr_max_k"$k"/*PRED* AGPr_max_k"$k"/; cp AGPr_max_k"$k"/PRED_OUTPUT_rel_clr_pca3counts_lodo_False.csv AGPr_max_k"$k"/PRED_OUTPUT_rel_clr_pca33counts_lodo_False.csv; done

for k in 5 6 7; do scp -r hoffy:/u/home/b/briscoel/project-halperin/MicroBatch/data/AGPr_max_k"$k"/*GRID* AGPr_max_k"$k"/; cp AGPr_max_k"$k"/GRID_AUC_OUTPUT_rel_clr_pca3counts_lodo_False.csv AGPr_max_k"$k"/GRID_AUC_OUTPUT_rel_clr_pca33counts_lodo_False.csv; done


for k in 5; do scp -r hoffy:/u/home/b/briscoel/project-halperin/MicroBatch/data/AGPr_complete_otu/*GRID* AGPr_complete_otu/; cp AGPr_complete_otu/GRID_AUC_OUTPUT_rel_clr_pca3counts_lodo_False.csv AGPr_complete_otu/GRID_AUC_OUTPUT_rel_clr_pca33counts_lodo_False.csv; done



for k in 6 7; do scp -r hoffy:/u/home/b/briscoel/project-halperin/MicroBatch/data/Thomasr_max_k"$k"/*GRID* Thomasr_max_k"$k"/; cp Thomasr_max_k"$k"/GRID_AUC_OUTPUT_rel_clr_pca3counts_lodo_False.csv Thomasr_max_k"$k"/GRID_AUC_OUTPUT_rel_clr_pca33counts_lodo_False.csv; done

for k in 6; do scp -r hoffy:/u/home/b/briscoel/project-halperin/MicroBatch/data/Thomasr_complete_otu/*GRID* Thomasr_complete_otu/; cp Thomasr_complete_otu/GRID_AUC_OUTPUT_rel_clr_pca3counts_lodo_False.csv Thomasr_complete_otu/GRID_AUC_OUTPUT_rel_clr_pca33counts_lodo_False.csv; done




# ESGE R

for k in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations_edgeR_deSeq2.R -m 16 -t 24 -hp -v 3.6.0 -arg Gibbonsr_max_k"$k"; done

for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations_edgeR_deSeq2.R -m 16 -t 24 -hp -v 3.6.0 -arg Kaplanr_max_k"$k"; done

for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations_edgeR_deSeq2.R -m 16 -t 24 -hp -v 3.6.0 -arg Thomasr_max_k"$k"; done

/u/local/apps/submit_scripts/R_job_submitter.sh -n transformations_edgeR_deSeq2.R -m 16 -t 16 -hp -v 3.6.0 -arg Thomasr_complete_otu

k=5
Rscript transformations_edgeR_deSeq2.R Gibbonsr_max_k"$k"

 ./qsub_regression.sh $index Kaplanr_complete_otu $method 0 bmi_v2
 
python regression.py --folder Kaplanr_max_k5 --trans vsd --correction nocorrection --lodo 1 --phenotype bmi_v2 




