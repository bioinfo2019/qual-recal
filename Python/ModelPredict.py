#!/home/eliot/anaconda2/bin/python

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.externals import joblib
import itertools
import math
from sklearn import metrics

from collections import Counter

from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import ensemble
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn import preprocessing
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.metrics import confusion_matrix
from sklearn.isotonic import IsotonicRegression as IR
from sklearn.metrics import brier_score_loss, f1_score
from sklearn.metrics import roc_curve
from sklearn.metrics import recall_score, precision_score
from sklearn.model_selection import RandomizedSearchCV
from sklearn.utils import check_random_state, safe_indexing
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.metrics import average_precision_score

#from sklearn.externals.six import string_types
from sklearn.model_selection import cross_val_predict
from imblearn.under_sampling import TomekLinks
from imblearn.metrics import geometric_mean_score

# Set to be the same as the folder set in ModelTraining.py
BASE_DIR = '/media/eliot/WD4TB/eliot_files/MyFiles/PhD'

MAPQ = 3
INSERT_SIZE = 4
GAP_EXTENSIONS = 5
EDIT_DISTANCE = 6
MISMATCHES = 7
GAP_OPENS = 8
ALIGNMENT_SCORE = 9
SECONDARY_ALIGNMENT_SCORE = 10
READ_COMPLEXITY_SCORE = 11
READ_GC_CONTENT = 12
PAIR_ORIENTATION = 13
PAIR_ALIGNMENT_TYPE = 14
INTERCEPT = 15
SLOPE = 16
R_VALUE = 17
DEPTH = 18
READ_EIGENVAL_1 = 19
READ_EIGENVAL_2 = 20
READ_TRACE = 21
READ_DET = 22
REF_EIGENVAL_1 = 23
REF_EIGENVAL_2 = 24
REF_TRACE = 25
REF_DET = 26
REF_COMPLEXITY_SCORE = 27
REF_GC_CONTENT = 28
N_LOW_QUAL_BASES = 29
AVG_LOW_QUAL_BASE_SCORE = 30
MATE_ALIGNMENT_SCORE = 31
ALIGNMENT_SCORES_DIFF = 32
CLASS = 33


def parallel_predict(X, est, n_chunks=8):

    # Split data array in 5 chunks

    n_samples = len(X)
    
    slices = [(int(n_samples*i/n_chunks), int(n_samples*(i+1)/n_chunks)) for i in range(n_chunks)]
    data_chunks = [X[i[0]:i[1]] for i in slices]
    
    # Setup 5 parallel jobs 
    jobs = (joblib.delayed(est.predict)(array) for array in data_chunks)
    parallel = joblib.Parallel(n_jobs=n_chunks)
    
    # Run jobs: works
    results = parallel(jobs)
    return np.asarray([y for x in results for y in x])


def parallel_predict_proba(X, est, n_chunks=8):

    # Split data array in 5 chunks

    n_samples = len(X)
    
    slices = [(int(n_samples*i/n_chunks), int(n_samples*(i+1)/n_chunks)) for i in range(n_chunks)]
    data_chunks = [X[i[0]:i[1]] for i in slices]
    
    # Setup 5 parallel jobs 
    jobs = (joblib.delayed(est.predict_proba)(array) for array in data_chunks)
    parallel = joblib.Parallel(n_jobs=n_chunks)
    
    # Run jobs: works
    results = parallel(jobs)
    return np.asarray([y for x in results for y in x])

def parallel_decision_function(X, est, n_chunks=8):

    # Split data array in 5 chunks

    n_samples = len(X)
    
    slices = [(int(n_samples*i/n_chunks), int(n_samples*(i+1)/n_chunks)) for i in range(n_chunks)]
    data_chunks = [X[i[0]:i[1]] for i in slices]
    
    # Setup 5 parallel jobs 
    jobs = (joblib.delayed(est.decision_function)(array) for array in data_chunks)
    parallel = joblib.Parallel(n_jobs=n_chunks)
    
    # Run jobs: works
    results = parallel(jobs)
    return np.asarray([y for x in results for y in x])

def write_scores_file(read_names, out_file_name, prob):

    data = []
   
    for i in range(len(prob)):

        if prob[i] < 0.00000001:
            prob[i] = 0.00001
            
        phred = int(-10*math.log(prob[i], 10))

#        if y_pred[i] == 1 and phred > mapqs[i]:
#            phred = int(mapqs[i])
#        if y_pred[i] == 0 and phred < mapqs[i]:
#            phred = int(mapqs[i])
                
        line = []
        line.append(read_names[i][0])
        line.append(phred)

        
        data.append(line)
        
    #out.write(dataset[i] + '\t' + str(phred) + '\t' + str(prob[i]) + '\t' + str(1.0 - prob[i]) + '\t' + str(y_true[i]) + '\t' + str(y_pred[i]) + '\n')
    
    df = pd.DataFrame(data)
    df.to_csv(out_file_name, sep='\t', header=False)
    
################## Run Parameters ######################################

feats_file = BASE_DIR + '/varsim_run/xaa'

#cols = [MAPQ,INSERT_SIZE,GAP_EXTENSIONS,EDIT_DISTANCE,MISMATCHES,GAP_OPENS,ALIGNMENT_SCORE,SECONDARY_ALIGNMENT_SCORE,PAIR_ORIENTATION,PAIR_ALIGNMENT_TYPE,INTERCEPT,SLOPE,N_LOW_QUAL_BASES,CLASS]
cols = (SECONDARY_ALIGNMENT_SCORE,ALIGNMENT_SCORES_DIFF,ALIGNMENT_SCORE,MAPQ,MATE_ALIGNMENT_SCORE,MISMATCHES,EDIT_DISTANCE,CLASS)

est_type = 'ADB'

########################################################################



print "loading data file...."

dataset = pd.read_csv(feats_file, usecols=cols, delimiter="\t").values
read_names = pd.read_csv(feats_file, usecols=(0,), delimiter="\t").values

print "done\n"        

last_col = len(cols) - 1
X = dataset[:,0:last_col]
y = dataset[:,last_col]
y = y.astype(int)

#mapqs = X[:, 0]
ir = IR( out_of_bounds = 'clip' )


print "scaling input data...."
min_max_scaler = preprocessing.MinMaxScaler()
Xs = min_max_scaler.fit_transform(X)
print "done\n"

X = 0

print "Predicting and calibrating probabilities...\n"

clf = joblib.load(BASE_DIR + '/models/' + est_type + '.pkl')
ir = joblib.load(BASE_DIR + '/models/' + est_type + '_ISO.pkl')
ir = IR( out_of_bounds = 'clip' )
#class_pred = parallel_predict(Xs, clf, 2)
class_pred = clf.predict(Xs)

if est_type != 'RUSVM' and est_type != 'RUSVMS':
    #proba_pred = parallel_predict_proba(Xs, clf, 2)[:,1]
    proba_pred = clf.predict_proba(Xs)[:,1]
    p_cal = ir.transform(proba_pred)
else:
    distance = clf.decision_function(Xs)
    isr_fit = ir.fit(distance, class_pred)
    p_cal = ir.transform(distance)
    proba_pred = p_cal
    
cm_preds_no_cal = confusion_matrix(y, class_pred)



test_preds_cal = p_cal > 0.5
test_preds_cal = test_preds_cal.astype(int)
cm_preds_cal = confusion_matrix(y, test_preds_cal)


print "Test stats (no recal):"
print "Confusion matrix:"
print cm_preds_no_cal
print "precision:"
print precision_score(y, class_pred, average='micro')
print "Average precision:"
if est_type == 'RUSVM' or est_type == 'RUSVMS':
    print average_precision_score(y, distance, average='weighted')
else:
    print average_precision_score(y, proba_pred, average='weighted') 


print "recall:"
print recall_score(y, class_pred, average='micro')
print "gmean:"
print geometric_mean_score(y, class_pred, average='weighted')
print "F1 score:"
print f1_score(y, class_pred)
print "Brier score:"
print brier_score_loss(y, proba_pred)


print ""

print "Test stats (recal):"
print "Confusion matrix:"
print cm_preds_cal
print "precision:"
print precision_score(y, test_preds_cal, average='micro')
print "Average precision:"
print average_precision_score(y, p_cal, average='weighted') 
print "recall:"
print recall_score(y, test_preds_cal, average=None)
print "gmean:"
print geometric_mean_score(y, test_preds_cal, average=None)
print "gmean (weighted):"
print geometric_mean_score(y, test_preds_cal, average='weighted')
print "F1 score:"
print f1_score(y, class_pred)
print "Brier score:"
print brier_score_loss(y, p_cal)

print "\nWriting scores file..."

if est_type == 'RUSVM' or est_type == 'RUSVMS':
    write_scores_file(read_names, BASE_DIR + '/pepper_scores_' + est_type + '_1.tab', proba_pred)
else:
    write_scores_file(read_names, BASE_DIR + '/pepper_scores_' + est_type + '_1.tab', p_cal) 
    
print "Done!"
































