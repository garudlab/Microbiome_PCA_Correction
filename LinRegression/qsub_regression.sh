#!/bin/bash


# ./qsub_classifier.sh 1 Thomasr_complete_otu nocorrection 0 bin_crc_normal 
#  ./qsub_classifier.sh 500 Thomasr_complete_otu nocorrection 1 bin_crc_normal
# ./qsub_classifier.sh 1000 Thomasr_complete_otu clr_pca 0 bin_crc_normal
# ./qsub_classifier.sh 1500 Thomasr_complete_otu clr_pcacounts 0 bin_crc_normal

# ./qsub_classifier.sh 1 AGPr_complete_otu nocorrection 0 in_antibiotic_last_year






first_count_input=$1
dataset_input=$2
rel=$3
corr_input=$4
lodo_input=$5
phen_input=$6




COUNTER=$first_count_input-1
echo $COUNTER
COUNTER=$((COUNTER + 1)); 
echo $COUNTER; 

# if [[ "$corr_input" == *"vst"* || "$corr_input" == *"cpm"* ]] ; then
# 	echo "$dataset_input $corr_input nocorrection $lodo_input $phen_input" > inputs/data_$COUNTER.in; 
# else
# 	echo "$dataset_input rel $corr_input $lodo_input $phen_input" > inputs/data_$COUNTER.in; 	
# fi
echo "$dataset_input $rel $corr_input $lodo_input $phen_input" > inputs/data_$COUNTER.in; 	




if [[ "$dataset_input" == *"AGPr_max_k5"* ]] ; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N Lin -l h_data=4G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_regression.sh"

elif [[ "$dataset_input" == *"AGPr_max_k6"* ]] ; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N Lin -l h_data=12G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./run_regression.sh"

elif [[ "$dataset_input" == *"AGPr_"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N Lin -l h_data=18G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_regression.sh"

elif [[ "$dataset_input" == *"k5"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N Lin -l h_data=6G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./run_regression.sh"


elif [[ "$dataset_input" == *"k7"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N Lin -l h_data=14G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./run_regression.sh"


elif [[ "$dataset_input" == *"k8"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N Lin -l h_data=20G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./run_regression.sh"

elif [[ "$dataset_input" == *"k6"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N Lin -l h_data=12G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./run_regression.sh"

else
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N Lin -l h_data=5G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./run_regression.sh"
fi




