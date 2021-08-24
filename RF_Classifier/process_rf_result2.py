
##python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction nocorrection --lodo 0 --phenotype bin_crc_normal
##python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction clr_pcacounts --lodo 0 --phenotype bin_crc_normal
##python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction clr_pca --lodo 0 --phenotype bin_crc_normal
##python process_rf_result2.py --folder AGPr_max_k7 --trans rel --correction nocorrection --lodo 0 --phenotype bin_antibiotic_last_year

import argparse,sys
import classifier_utils
import data_utils
import os
import pandas as pd
import numpy as np
import random
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn import model_selection 
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score
from collections import Counter
import os.path
from os import path
parser=argparse.ArgumentParser()

parser.add_argument('--folder', help='Name of dataset folder',type=str)
parser.add_argument('--trans', help='Which transformation of data',type=str)
parser.add_argument('--correction', help='Which correction method',type=str)
parser.add_argument('--lodo', help='Whether to use lodo (1 or 0)',type=int)
parser.add_argument('--phenotype', help='What phenotype are you predicting',type =str)
args=parser.parse_args()
folder = args.folder # "AGPr_max_k5" #"AGPr_complete_otu" 
trans = args.trans #"rel"
correction = args.correction
lodo = bool(args.lodo)
phenotype = args.phenotype


main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"
local = False

if local:
	main_dir = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
	script_folder = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"



data_dir = main_dir + folder + "/"
out_dir = data_dir + "RF/"

metrics = pd.DataFrame(columns=["nest", "crit","maxd","miss","misl","maf","val_auc","val_f1","test_auc","test_f1"]) 

for nest in [100,1000, 1500]: 
	for crit in ["entropy", "gini"] : 
		for maxd in [1,2, 3]: 
			for miss in [2 , 5,  10]:
				for misl in [1, 5, 10]:
					for maf in [0.1, 0.3, 0.5]: 
						output_string  = "_lodo_" + str(lodo) + "_nest_" + str(nest) + "_crit_" + str(crit) + "_maxd_" + str(maxd) + "_miss_" + str(miss) + \
						"_misl_" + str(misl) + "_maf_" + str(maf) 
						file_path = out_dir  + "GRID_" + trans + "_" + correction  + output_string  + ".pkl"
						if not path.isfile(file_path):
							print("Path missing: " + str(file_path))
						else:
							
							model_result = pickle.load( open( file_path ,"rb"))

							if lodo:
								lodo_dict = {"nest":nest,"crit":crit,"maxd":maxd,"miss":miss,"misl":misl, \
									"maf":maf,"val_auc":np.mean(model_result['val_auc']),"val_f1":np.mean(model_result['val_f1']), \
									"test_auc":np.mean(model_result['test_auc']),"test_f1":np.mean(model_result['test_f1'])}
								for test in range(len(model_result["test_dataset"])):
									lodo_dict["AUC" + model_result["test_dataset"][test] ] = model_result['test_auc'][test]
									lodo_dict["F1" + model_result["test_dataset"][test] ] = model_result['test_f1'][test]

								metrics = metrics.append(lodo_dict,ignore_index=True)


							else:
								reg_dict = {"nest":nest,"crit":crit,"maxd":maxd,"miss":miss,"misl":misl, \
									"maf":maf,"val_auc":np.mean(model_result['val_auc']),"val_f1":np.mean(model_result['val_f1']), "test_auc":np.mean(model_result['test_auc']),"test_f1":np.mean(model_result['test_f1']) }

								for test in range(len(model_result["test_auc"])):
									reg_dict["AUC_" + str(test) ] = model_result['test_auc'][test]
									reg_dict["F1"+ str(test) ] = model_result['test_f1'][test]

								metrics = metrics.append(reg_dict,ignore_index=True)
print(metrics)
best_model_index = list(metrics["val_auc"]).index(max(list(metrics["val_auc"])))
print("best_model_index" + str(best_model_index))
print(metrics.iloc[best_model_index])

metrics.to_csv(data_dir + "GRID_OUTPUT_" + trans + "_" + correction  + "_lodo_" + str(lodo)  + ".csv",index=False)
