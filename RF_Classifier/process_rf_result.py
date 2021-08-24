
##python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction nocorrection --lodo 1 --phenotype bin_crc_normal
##python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction clr_pca4counts --lodo 1 --phenotype bin_crc_normal
##python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction clr_scale_pca --lodo 1--phenotype bin_crc_normal
##python process_rf_result.py --folder AGPr_max_k6 --trans rel --correction nocorrection --lodo 0 --phenotype bin_antibiotic_last_year --metric val_auc

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
parser.add_argument('--metric', help='What metric do you want to use',type =str)
args=parser.parse_args()
folder = args.folder # "AGPr_max_k5" #"AGPr_complete_otu" 
trans = args.trans #"rel"
correction = args.correction
lodo = bool(args.lodo)
phenotype = args.phenotype
meas = args.metric


main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"
local = False

if local:
	main_dir = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
	script_folder = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"



studies_dict = {"Thomasr_complete_otu": ['FengQ_2015','HanniganGD_2017', 'ThomasAM_2018a', 'ThomasAM_2018b', 'VogtmannE_2016', 'YuJ_2015', 'ZellerG_2014'],
"Thomasr_max_k6": ['FengQ_2015','HanniganGD_2017', 'ThomasAM_2018a', 'ThomasAM_2018b', 'VogtmannE_2016', 'YuJ_2015', 'ZellerG_2014'],
"Gibbonsr_complete_otu": ["crc_zeller","crc_baxter","crc_zackular"],
"Kaplanr_complete_otu": ["HOWE_KF1", "HOWE_KF2", "HOWE_KF3", "HOWE_KF4"],
"Kaplanr_max_k6": ["HOWE_KF1", "HOWE_KF2", "HOWE_KF3", "HOWE_KF4"]}
data_dir = main_dir + folder + "/"
out_dir = data_dir + "RF/"

metrics_dict = dict()
if lodo:
	for s in studies_dict[folder]:
		metrics_dict[s] = pd.DataFrame(columns=["nest", "crit","maxd","miss","misl","maf","val_auc","val_f1","test_auc","test_f1"]) 

else:
	for s in range(50):
		metrics_dict[s] = pd.DataFrame(columns=["nest", "crit","maxd","miss","misl","maf","val_auc","val_f1","test_auc","test_f1"]) 


#metrics = pd.DataFrame(columns=["nest", "crit","maxd","miss","misl","maf","val_auc","val_f1","test_auc","test_f1"]) 






for nest in [100,1000, 1500]: 
	for crit in ["entropy"] : 
		for maxd in ["None"]: 
			for miss in [2 , 5,  10]:
				for misl in [1, 5, 10]:
					for maf in ["auto"]: 
						output_string  = "_lodo_" + str(lodo) + "_nest_" + str(nest) + "_crit_" + str(crit) + "_maxd_" + str(maxd) + "_miss_" + str(miss) + \
						"_misl_" + str(misl) + "_maf_" + str(maf) 
						file_path = out_dir  + "GRID_" + trans + "_" + correction  + output_string  + ".pkl"
						if not path.isfile(file_path):
							print("Path missing: " + str(file_path))
						else:
							
							model_result = pickle.load( open( file_path ,"rb"))
							#print(model_result)
							#sys.exit()

							if lodo:
								
								for test in range(len(model_result["test_dataset"])):
									metrics_dict[model_result["test_dataset"][test]] = metrics_dict[model_result["test_dataset"][test]].append({"nest":nest,"crit":crit,"maxd":maxd,"miss":miss,"misl":misl, \
									"maf":maf,"val_auc":model_result['val_auc'][test],"val_f1":model_result['val_f1'][test],
									"test_auc": model_result['test_auc'][test],"test_f1":model_result['test_f1'][test]},ignore_index=True)
							else:
								for test in range(len(model_result["test_auc"])):
									metrics_dict[test] = metrics_dict[test].append({"nest":nest,"crit":crit,"maxd":maxd,"miss":miss,"misl":misl, \
									"maf":maf,"val_auc":model_result['val_auc'][test],"val_f1":model_result['val_f1'][test],
									"test_auc": model_result['test_auc'][test],"test_f1":model_result['test_f1'][test]},ignore_index=True)

									

final_metrics = pd.DataFrame(columns=["fold", "nest", "crit","maxd","miss","misl","maf","val_auc","val_f1","test_auc","test_f1"]) 

print(metrics_dict)

for k in metrics_dict.keys():
	print("max")
	print(max(list(metrics_dict[k][meas])))
	best_model_index = list(metrics_dict[k][meas]).index(max(list(metrics_dict[k][meas])))
	temp = metrics_dict[k].iloc[best_model_index]
	#temp.insert(0,"fold",k, True)
	temp.set_value("fold",k)
	final_metrics = final_metrics.append(temp)


print(final_metrics)

if meas == "val_auc":
	final_metrics.to_csv(data_dir + "GRID_AUC_OUTPUT_" + trans + "_" + correction  + "_lodo_" + str(lodo)  + ".csv",index=False)
else:
	final_metrics.to_csv(data_dir + "GRID_F1_OUTPUT_" + trans + "_" + correction  + "_lodo_" + str(lodo)  + ".csv",index=False)





