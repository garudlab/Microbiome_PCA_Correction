#!/bin/bash
. /u/local/Modules/default/init/modules.sh
module load python/anaconda3

while read -r arg_1 arg_2 arg_3 arg_4 arg_5; do
	python regression.py --folder $arg_1 --trans $arg_2 --correction $arg_3 --lodo $arg_4 --phenotype $arg_5;
done < inputs/data_$SGE_TASK_ID.in

