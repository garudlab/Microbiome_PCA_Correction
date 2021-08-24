# python classifier.py --folder Gibbonsr_max_k5 --trans rel --correction nocorrection --lodo 1 --phenotype bin_crc_normal --n_estimators 100 --criterion entropy --max_depth None --min_samples_split 2 --min_samples_leaf 5 --max_features auto


# python classifier.py --folder Thomasr_max_k7 --trans rel --correction nocorrection --lodo 0 --phenotype bin_crc_normal --n_estimators 100 --criterion gini --max_depth None --min_samples_split 5 --min_samples_leaf 1 --max_features auto

# python classifier.py --folder Kaplanr_complete_otu --trans rel --correction clr --lodo 1 --phenotype diabetes_self_v2 --n_estimators 100 --criterion entropy --max_depth None --min_samples_split 5 --min_samples_leaf 1 --max_features auto
# python classifier.py --folder Thomasr_complete_otu --trans rel --correction nocorrection --lodo 1 --phenotype bin_crc_normal --n_estimators 100 --criterion entropy --max_depth 1 --min_samples_split 5 --min_samples_leaf 1 --max_features 0.1
# python classifier.py --folder AGPr_max_k5 --trans rel --correction clr_pca2counts --lodo 0 --phenotype bin_antibiotic_last_year --n_estimators 100 --criterion entropy --max_depth 1 --min_samples_split 5 --min_samples_leaf 1 --max_features 0.1

# python classifier.py --folder AGPr_complete_otu --trans rel --correction nocorrection --lodo 0 --phenotype bin_antibiotic_last_year --n_estimators 100 --criterion entropy --max_depth 1 --min_samples_split 5 --min_samples_leaf 1 --max_features 0.1


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
 "Kaplanr_max_k6": "extraction_robot..exp."}

parser.add_argument('--folder', help='Name of dataset folder',type=str)
parser.add_argument('--trans', help='Which transformation of data',type=str)
parser.add_argument('--correction', help='Which correction method',type=str)
parser.add_argument('--lodo', help='Whether to use lodo (1 or 0)',type=int)
parser.add_argument('--phenotype', help='What phenotype are you predicting',type =str)

parser.add_argument('--n_estimators', help='What is the minimum allele frequency for non mothers', type = int)
parser.add_argument('--criterion', help='Strain for this job', type = str)
parser.add_argument('--max_depth', help='Study to study',type=str)
parser.add_argument('--min_samples_split', help='End index', type = int)
parser.add_argument('--min_samples_leaf', help='Minimum read depth for a valid allele frequency',type = int)
parser.add_argument('--max_features', help='Max features per tree', type = str)
args=parser.parse_args()
print(args)



folder = args.folder # "AGPr_max_k5" #"AGPr_complete_otu" 
trans = args.trans #"rel"
correction = args.correction
lodo = bool(args.lodo)
phenotype = args.phenotype
param_n_estimators = args.n_estimators
param_criterion = args.criterion
if args.max_depth == "None":
	param_max_depth = None
else:
	param_max_depth = int(args.max_depth)
param_min_samples_split = args.min_samples_split
param_min_samples_leaf = args.min_samples_leaf
if args.max_features== "auto":
	param_max_features = args.max_features
else:
	param_max_features = float(args.max_features)


data_dir = main_dir + folder + "/"
out_dir = data_dir + "RF/"
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
rskf = model_selection.RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=123)
logo = LeaveOneGroupOut()
parameter_dict = {'n_estimators':[param_n_estimators],'criterion': [param_criterion],\
'min_samples_leaf': [param_min_samples_leaf],'max_features':[param_max_features],\
'min_samples_split': [param_min_samples_split],'max_depth':[param_max_depth]}

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
if "Kaplan" in folder:
	labels = np.array([np.nan if lab=="not provided" else int(lab) for lab in labels])
print(labels)

na_mask = pd.isna(labels)
feature_table = feature_table[~na_mask,:] # get rid of samples with na labels
labels = labels[~na_mask]
metadata_table = metadata_table.loc[~na_mask,:]  
print(Counter(labels))
## BINARIZE
if phenotype == "bin_antibiotic_last_year":
	labels = np.array([1 if lab=="Yes" else 0 for lab in labels])

if (phenotype == "bin_crc_normal") & (("Gibbonsr_complete" in folder) | ("Thomasr_max" in folder) ):
	labels = np.array([1 if lab=="CRC" else 0 for lab in labels])


print(Counter(labels))
#### CHECK
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
	results_dict['val_auc'] = []
	results_dict['test_auc'] = []
	results_dict['val_f1'] = []
	results_dict['test_f1'] = []
	results_dict["test_dataset"] = []  

else:
	results_dict['val_auc'] = []
	results_dict['test_auc'] = []
	results_dict['val_f1'] = []
	results_dict['test_f1'] = []
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
	# results_dict[train_iter]['val_auc_trained']= [] 
	# results_dict[train_iter]["number samples"] = []  

	X_train, X_test = feature_table[train_index,], feature_table[test_index,]

	
	y_train, y_test = labels[train_index], labels[test_index]

	## SUBSAMNPLE
	if phenotype in ["diabetes_self_v2","bin_antibiotic_last_year"]:
		print("impose balance")
		print(Counter(y_train))


		X_train,y_train = classifier_utils.balanced_subsample(X_train,y_train,subsample_size=1.0)
		print("after impose balance")
		print(Counter(y_train))

	#sys.exit()
	X_train, X_val, y_train, y_val = train_test_split(X_train, y_train,test_size=0.30, random_state=1) 


	if lodo:
		print("get test dataset name if lodo")
		results_dict["test_dataset"].append(metadata_table.iloc[test_index][groups[folder]][0])
		
	else:
		results_dict["train_shape"].append(X_train.shape[0])
		results_dict["val_shape"].append(X_val.shape[0])
		results_dict["test_shape"].append(X_test.shape[0])

	# perform grid search on train
	print("DIM train")
	print(X_train.shape)

	clf = RandomForestClassifier(random_state=0,n_estimators = param_n_estimators, criterion = param_criterion, max_depth = param_max_depth, 
		min_samples_split = param_min_samples_split, min_samples_leaf = param_min_samples_leaf, max_features = param_max_features)

	print("collections")
	#print(Counter(list(y_train)))
	print(y_train)
	clf.fit(X_train, y_train)

	test = [(est.get_depth(), est.tree_.max_depth, est.max_depth) for est in clf.estimators_]
	print(test)

	
	## TESTONLY 

	results_dict['val_auc'].append(roc_auc_score(y_val, clf.predict_proba(X_val)[:, 1]))
	results_dict['test_auc'].append(roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1]))
	results_dict['val_f1'].append(f1_score(y_val, clf.predict(X_val)))
	results_dict['test_f1'].append(f1_score(y_test, clf.predict(X_test)))


print(results_dict)	



output_string  = "_lodo_" + str(lodo) + "_nest_" + str(param_n_estimators) + "_crit_" + str(param_criterion) + "_maxd_" + str(param_max_depth) +  "_miss_" + str(param_min_samples_split) + \
"_misl_" + str(param_min_samples_leaf) + "_maf_" + str(param_max_features) 

pickle.dump(results_dict , open( out_dir  + "GRID_" + trans + "_" + correction  + output_string  + ".pkl", "wb" ) )





   