# python classifier.py --folder Gibbonsr_max_k5 --trans rel --correction nocorrection --lodo 1 --phenotype bin_crc_normal --n_estimators 100 --criterion entropy --max_depth None --min_samples_split 2 --min_samples_leaf 5 --max_features auto


# python classifier.py --folder Thomasr_max_k7 --trans rel --correction nocorrection --lodo 0 --phenotype bin_crc_normal --n_estimators 100 --criterion gini --max_depth None --min_samples_split 5 --min_samples_leaf 1 --max_features auto

# python classifier.py --folder Kaplanr_complete_otu --trans rel --correction clr --lodo 1 --phenotype diabetes_self_v2 --n_estimators 100 --criterion entropy --max_depth None --min_samples_split 5 --min_samples_leaf 1 --max_features auto
# python classifier.py --folder Thomasr_complete_otu --trans rel --correction nocorrection --lodo 1 --phenotype bin_crc_normal --n_estimators 100 --criterion entropy --max_depth 1 --min_samples_split 5 --min_samples_leaf 1 --max_features 0.1
# python classifier.py --folder AGPr_max_k5 --trans rel --correction clr_pca2counts --lodo 0 --phenotype bin_antibiotic_last_year --n_estimators 100 --criterion entropy --max_depth 1 --min_samples_split 5 --min_samples_leaf 1 --max_features 0.1

# python regression.py --folder AGPr_max_k5 --trans rel --correction nocorrection --lodo 0 --phenotype bmi_corrected 
# python regression.py --folder Kaplanr_max_k5 --trans rel --correction nocorrection --lodo 1 --phenotype bmi_v2 

import argparse,sys

import os
import pandas as pd
import numpy as np
import random
import pickle
from sklearn.linear_model import LinearRegression
from sklearn import model_selection 
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, f1_score
from collections import Counter
parser=argparse.ArgumentParser()

main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"
local = False

if local:
	main_dir = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
	script_folder = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"

groups = {"AGPr_complete_otu":"Instrument","Thomasr_complete_otu":"dataset_name","AGPr_max_k7":"Instrument", "AGPr_max_k5":"Instrument", \
"AGPr_max_k6":"Instrument", "Gibbonsr_complete_otu":"study", "Thomasr_max_k6":"dataset_name","Thomasr_max_k7":"dataset_name",\
 "Gibbonsr_max_k5":"study", "Gibbonsr_max_k6":"study", "Gibbonsr_max_k7":"study","Gibbonsr_max_k8":"study", "Kaplanr_complete_otu": "extraction_robot..exp.",
 "Kaplanr_max_k6": "extraction_robot..exp.",'Kaplanr_max_k5':"extraction_robot..exp."}

parser.add_argument('--folder', help='Name of dataset folder',type=str)
parser.add_argument('--trans', help='Which transformation of data',type=str)
parser.add_argument('--correction', help='Which correction method',type=str)
parser.add_argument('--lodo', help='Whether to use lodo (1 or 0)',type=int)
parser.add_argument('--phenotype', help='What phenotype are you predicting',type =str)


args=parser.parse_args()
print(args)



folder = args.folder # "AGPr_max_k5" #"AGPr_complete_otu" 
trans = args.trans #"rel"
correction = args.correction
lodo = bool(args.lodo)
phenotype = args.phenotype


data_dir = main_dir + folder + "/"
out_dir = data_dir + "PRED/"
if not os.path.isdir(out_dir):
	os.makedirs(out_dir)
	print("created folder : ", out_dir)

else:
	print(out_dir, "folder already exists.")

### IMPORT DATA
metadata_table = pd.read_csv(data_dir + "metadata.txt", delimiter = "\t",header=0)
if correction == "nocorrection":
	feature_table = pd.read_csv(data_dir  + "feature_table_" + trans + ".txt" ,delimiter = "\t",header=0)
else:
	feature_table = pd.read_csv(data_dir  +"feature_table_" + trans + "_" + correction + ".txt" ,delimiter = "\t",header=0)

### FIXED PARAM

n_splits = 5
n_repeats = 10
random.seed(567)
rskf = model_selection.RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=123)
logo = LeaveOneGroupOut()

# ## CHECK
print("before filter")
print(feature_table.shape)
print(len(metadata_table[phenotype]))

### POLISH OUTCOME VARIABLE
#non_nan_samples = metadata_table.index[np.invert(np.isnan(metadata_table[phenotype]))]
#non_nan_samples = intersection(feature_table.columns,non_nan_samples)
#feature_table = feature_table[non_nan_samples]
#metadata_labels = (metadata_table.loc[non_nan_samples])[phenotype]


### PREPARE FOR PREDICION
feature_table = np.array(feature_table).transpose()
labels = np.array(metadata_table[phenotype])
if "bmi" in phenotype:
	labels = np.array([np.nan if lab in ["not provided","not applicable","Not provided","Unspecified"] else float(lab) for lab in labels])

na_mask = pd.isna(labels)
feature_table = feature_table[~na_mask,:] # get rid of samples with na labels
labels = labels[~na_mask]
metadata_table = metadata_table.loc[~na_mask,:]  

print("after filter")
print(feature_table.shape)
print(len(labels))
print("metadata shape")
print(metadata_table.shape)

#### CROSS VALIDTION


if lodo:
	logo = LeaveOneGroupOut()

	print("group cheeck")
	print(Counter(metadata_table[groups[folder]]))

	splitter = logo.split(feature_table, labels, np.array(metadata_table[groups[folder]]))   
else:
	splitter = rskf.split(feature_table, labels)

results_dict = dict()
  

train_iter = 0
results_dict= dict()

if lodo:
	results_dict['val_score'] = []
	results_dict['test_score'] = []
	results_dict['val_pearson'] = []
	results_dict['test_pearson'] = []
	results_dict["test_dataset"] = []  

else:
	results_dict['val_score'] = []
	results_dict['test_score'] = []
	results_dict['val_pearson'] = []
	results_dict['test_pearson'] = []
	results_dict["train_shape"] = []  
	results_dict["val_shape"] = []  
	results_dict["test_shape"] = []  
for train_index, test_index in splitter: 
	# results_dict[train_iter] = dict()
	# results_dict[train_iter]['train_best_params'] = dict()
	# results_dict[train_iter]['train_auc_trained'] = []
	# results_dict[train_iter]['mean_train_cv_auc'] = []
	# results_dict[train_iter]['mean_test_cv_auc'] = []
	# results_dict[train_iter]['test_auc_trained'] = []
	# results_dict[train_iter]['val_score_trained']= [] 
	# results_dict[train_iter]["number samples"] = []  

	X_train, X_test = feature_table[train_index,], feature_table[test_index,]

	
	y_train, y_test = labels[train_index], labels[test_index]

	## SUBSAMNPLE
	

	#sys.exit()
	X_train, X_val, y_train, y_val = train_test_split(X_train, y_train,test_size=0.30, random_state=1) 


	if lodo:
		print("get test dataset name if lodo")
		print(metadata_table.iloc[test_index][groups[folder]][0])
		results_dict["test_dataset"].append(metadata_table.iloc[test_index][groups[folder]][0])
		
	else:
		results_dict["train_shape"].append(X_train.shape[0])
		results_dict["val_shape"].append(X_val.shape[0])
		results_dict["test_shape"].append(X_test.shape[0])

	# perform grid search on train
	print("DIM train")
	print(X_train.shape)

	clf =  LinearRegression()
	print("collections")
	#print(Counter(list(y_train)))
	print(y_train)
	clf.fit(X_train, y_train)


	
	## TESTONLY 
	results_dict['val_score'].append(clf.score(X_val,y_val))
	results_dict['test_score'].append(clf.score(X_test,y_test))
	results_dict['val_pearson'].append(np.corrcoef(y_val, clf.predict(X_val))[0,1])
	results_dict['test_pearson'].append(np.corrcoef(y_test, clf.predict(X_test))[0,1])


print(results_dict)	



output_string  = "_linreg_" 

metrics_dict = pd.DataFrame(columns=["test_dataset", "val_score","val_pearson","test_score","test_pearson"]) 



for k in range(len(results_dict['test_pearson'])):
	if lodo:

		metrics_dict= metrics_dict.append({"test_dataset":results_dict["test_dataset"][k], "val_score":results_dict['val_score'][k],
			"test_score":results_dict['test_score'][k],"val_pearson":results_dict['val_pearson'][k],
			"test_pearson":results_dict['test_pearson'][k]},ignore_index=True)
	else:

		metrics_dict= metrics_dict.append({"test_dataset":k, "val_score":results_dict['val_score'][k],
			"test_score":results_dict['test_score'][k],"val_pearson":results_dict['val_pearson'][k],
			"test_pearson":results_dict['test_pearson'][k]},ignore_index=True)
print(metrics_dict)

metrics_dict.to_csv(data_dir + "PRED_OUTPUT_" + trans + "_" + correction  + "_lodo_" + str(lodo)  + ".csv",index=False)





   