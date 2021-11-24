#!/bin/bash
. /u/local/Modules/default/init/modules.sh
module load python

while read -r arg_1 arg_2 arg_3 arg_4 arg_5 arg_6 arg_7 arg_8 arg_9 arg_10 arg_11; do
	python classifier.py --folder $arg_1 --trans $arg_2 --correction $arg_3 --lodo $arg_4 --phenotype $arg_5 --n_estimators $arg_6 --criterion $arg_7 --max_depth $arg_8 --min_samples_split $arg_9 --min_samples_leaf $arg_10 --max_features $arg_11;
done < /u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/RF_Classifier/inputs/data_$SGE_TASK_ID.in

