#!/home/eliot/anaconda2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 18:34:04 2018

@author: eliot
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.externals import joblib
import itertools
import math
from sklearn import metrics

from collections import Counter
from sklearn.kernel_approximation import Nystroem
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import ensemble
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC, LinearSVC
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
from imblearn.ensemble import BalancedBaggingClassifier, RUSBoostClassifier

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


def invert_phred(mapq):
    
    p = []

    for i in range(len(mapq)):
        if mapq[i] == 0:
            mapq[i] = 0.000001
            
        p.append(10.0**(-1*float(mapq[i])/10.0))
        
    return np.array(p)
    
    
def plot_calibration_curve(name, name2, name3, y, y2, y3, prob, prob2, prob3, fig_index, invert_phred):

    
    fig = plt.figure(fig_index, figsize=(10, 10))
    ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    ax2 = plt.subplot2grid((3, 1), (2, 0))
    
    ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
    
    if invert_phred == 1:
        for i in range(len(prob)):
            if prob[i] == 0:
                prob[i] = 0.00001
            prob[i] = 10.0**(-1*float(prob[i])/10.0)
    
    fraction_of_positives, mean_predicted_value = calibration_curve(y, prob, n_bins=10)
    
    clf_score = brier_score_loss(y, prob, pos_label=1)
    
    ax1.plot(mean_predicted_value, fraction_of_positives, marker="s", c="r", ls="-",
        label="%s (%1.5f)" % (name, clf_score))

    fraction_of_positives, mean_predicted_value = calibration_curve(y2, prob2, n_bins=10)
    
    clf_score = brier_score_loss(y2, prob2, pos_label=1)
    
    ax1.plot(mean_predicted_value, fraction_of_positives, marker="s", c="b", ls="-",
        label="%s (%1.5f)" % (name2, clf_score))

    fraction_of_positives, mean_predicted_value = calibration_curve(y3, prob3, n_bins=10)
    
    clf_score = brier_score_loss(y3, prob3, pos_label=1)
    
    ax1.plot(mean_predicted_value, fraction_of_positives, marker="s", c="g", ls="-",
        label="%s (%1.5f)" % (name3, clf_score))
        
    ax2.hist(prob, range=(0, 1), bins=10, label="Probability Histogram",
        histtype="step", lw=2)
        
    ax1.set_ylabel("Fraction of positives")
    ax1.set_ylim([-0.05, 1.05])
    ax1.legend(loc="upper left")
    ax1.set_title('Calibration plots  (reliability curve)')
    
    ax2.set_xlabel("Mean predicted value")
    ax2.set_ylabel("Count")
    ax2.legend(loc="upper center", ncol=2)
    
    plt.tight_layout()
    



def predict_tp_fp(X, y, prob):
    
    pred_pos_idxs = prob > 0.5
    
    est = LogisticRegression(C=100, class_weight='balanced')
    tl = TomekLinks(return_indices=True, n_jobs=24)
    X_r, y_r, idx_resampled = tl.fit_sample(X[pred_pos_idxs], y[pred_pos_idxs])
    #p_grid = {'C': [100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']}
    #est = SVC(kernel="rbf", probability=True)
    
    #est2 = SVC(kernel="rbf")
    grid_search = GridSearchCV(est, param_grid=p_grid, cv=5, n_jobs=24, scoring='f1')
    grid_search.fit(X_r, y_r)
    model = grid_search.best_estimator_.fit(X_r, y_r)

    pred_train = model.predict(X[pred_pos_idxs])
    
    print confusion_matrix(y[pred_pos_idxs], pred_train)
    
    return model.predict_proba(X[pred_pos_idxs])[:,1], pred_pos_idxs
    
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

    
def fit(estimator, n_max_subset, X, y, random_state=None, return_indices=True, under_sampling_factor=1, parameter_grid=None, refit=False):
    """Resample the dataset.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Matrix containing the data which have to be sampled.

    y : array-like, shape (n_samples,)
        Corresponding label for each sample in X.

    Returns
    -------
    X_resampled : {ndarray, sparse matrix}, shape \
(n_subset, n_samples_new, n_features)
        The array containing the resampled data.

    y_resampled : ndarray, shape (n_subset, n_samples_new)
        The corresponding label of `X_resampled`

    idx_under : ndarray, shape (n_subset, n_samples, )
        If `return_indices` is `True`, a boolean array will be returned
        containing the which samples have been selected.

    """
    
    models = []
    f1_scores = []
    feature_counts = []
    
    random_state = check_random_state(random_state)

    # array to know which samples are available to be taken
    samples_mask = np.ones(y.shape, dtype=bool)
    
    
    # where the different set will be stored
    idx_under = []

    n_subsets = 0
    b_subset_search = True
    while b_subset_search:
        target_stats = Counter(safe_indexing(
            y, np.flatnonzero(samples_mask)))
        ratio_ = {0:target_stats[1] * under_sampling_factor}
        # store the index of the data to under-sample
        index_under_sample = np.empty((0, ), dtype=y.dtype)
        # value which will be picked at each round
        index_constant = np.empty((0, ), dtype=y.dtype)
        for target_class in target_stats.keys():
            if target_class in ratio_.keys():
                n_samples = ratio_[target_class]
                # extract the data of interest for this round from the
                # current class
                index_class = np.flatnonzero(y == target_class)
                index_class_interest = index_class[samples_mask[
                    y == target_class]]
                y_class = safe_indexing(y, index_class_interest)
                # select randomly the desired features
                index_target_class = random_state.choice(
                    range(y_class.size), size=n_samples, replace=False)
                index_under_sample = np.concatenate(
                    (index_under_sample,
                     index_class_interest[index_target_class]),
                    axis=0)
            else:
                
                index_constant = np.concatenate(
                    (index_constant,
                     np.flatnonzero(y == target_class)),
                    axis=0)

        # store the set created
        n_subsets += 1
        subset_indices = np.concatenate((index_under_sample,
                                         index_constant), axis=0)
        idx_under.append(subset_indices)

        # fit and predict using cross validation
        X_subset = safe_indexing(X, subset_indices)
        
        y_subset = safe_indexing(y, subset_indices)
        
        grid_search = GridSearchCV(estimator, param_grid=parameter_grid, cv=5, scoring='f1_weighted', n_jobs=24)
        grid_search.fit(X_subset, y_subset)
        model = grid_search.best_estimator_.fit(X_subset, y_subset)                
        
        pred = model.predict(X_subset)
        
        if refit == True:
            
            fp_idxs = np.where(np.logical_and(y_subset == 0, pred == 1))[0]
            tp_idxs = np.where(np.logical_and(y_subset == 1, pred == 1))[0]
            tn_idxs = np.where(np.logical_and(y_subset == 0, pred == 0))[0]
            fn_idxs = np.where(np.logical_and(y_subset == 1, pred == 0))[0]
            
            
            n_tn = 1 #len(fn_idxs) - 1
            n_fn = 1 #len(fn_idxs) - 1
            n_fp = len(fp_idxs) - 1
            n_tp = (len(fp_idxs) * under_sampling_factor) - 1
            
            tn_pos_idxs = np.random.randint(0, len(tn_idxs)-1, size=n_tn)
            tp_pos_idxs = np.random.randint(0, len(tp_idxs)-1, size=n_tp)
            fp_pos_idxs = np.random.randint(0, len(fp_idxs)-1, size=n_fp)
            fn_pos_idxs = np.random.randint(0, len(fn_idxs)-1, size=n_fn)
        
        
            X_fp = X_subset[fp_idxs[fp_pos_idxs]]
            X_tp = X_subset[tp_idxs[tp_pos_idxs]]
            X_fn = X_subset[fn_idxs[fn_pos_idxs]]
            X_tn = X_subset[tn_idxs[tn_pos_idxs]]
        
            y_fp = y_subset[fp_idxs[fp_pos_idxs]]
            y_tp = y_subset[tp_idxs[tp_pos_idxs]]
            y_fn = y_subset[fn_idxs[fn_pos_idxs]]
            y_tn = y_subset[tn_idxs[tn_pos_idxs]]
        
            
            X_new = np.concatenate((X_fn, X_tn, X_fp, X_tp))
            y_new = np.concatenate((y_fn, y_tn, y_fp, y_tp))
                       
            LogisticRegression(class_weight='balanced', solver='liblinear')
    
            #est = SVC(kernel='rbf', class_weight='balanced', gamma='scale')
            grid_search = GridSearchCV(est, param_grid=parameter_grid, cv=5, scoring='f1_weighted', n_jobs=24)
            grid_search.fit(X_new, y_new)
            model = grid_search.best_estimator_.fit(X_new, y_new)
            feature_counts.append((n_fp, n_tp))
            
        else:
            
            feature_counts.append((np.count_nonzero(y_subset == 0), np.count_nonzero(y_subset == 1)))
            
            
        pred_train = parallel_predict(X_train, model, 8)
        f1_scores.append(f1_score(y_train, pred_train))

        #coeff_tally += model.coef_
        models.append(model)
        
        
        #print coeff_tally.shape
        #print model.coef_.shape
        
        #pred = cross_val_predict(estimator, X_subset, y_subset)
        # extract the prediction about the targeted classes only
        pred_target = pred[:index_under_sample.size]
        index_classified = index_under_sample[
            pred_target == safe_indexing(y_subset,
                                         range(index_under_sample.size))]
        #print np.unique(y_subset)
        #print pred_target.size
        
        samples_mask[index_classified] = False

        # check the stopping criterion
        if n_max_subset is not None:
            if n_subsets == n_max_subset:
                b_subset_search = False
        # check that there is enough samples for another round
        target_stats = Counter(safe_indexing(
            y, np.flatnonzero(samples_mask)))
        for target_class in ratio_.keys():
            if target_stats[target_class] < ratio_[target_class]:
                b_subset_search = False

                
    X_resampled, y_resampled = [], []
    for indices in idx_under:
        X_resampled.append(safe_indexing(X, indices))
        y_resampled.append(safe_indexing(y, indices))


    max_f1 = max(f1_scores)
    max_f1_idx = f1_scores.index(max_f1)
    model = models[max_f1_idx]
    
    if return_indices:
        return (np.array(X_resampled), np.array(y_resampled),
                np.array(idx_under), model, feature_counts[max_f1_idx])
    else:
        return np.array(X_resampled), np.array(y_resampled), model, feature_counts[max_f1_idx]



def adaboost_clf(X, y):

    bdt = AdaBoostClassifier(DecisionTreeClassifier(max_depth=1),
                         algorithm="SAMME",
                         n_estimators=200)

    #bdt.fit(X, y)
    
    return bdt


def gradientboost_clf(X, y):


    gbes = ensemble.GradientBoostingClassifier(n_estimators=1000,
                                               validation_fraction=0.2,
                                               n_iter_no_change=5, tol=0.01,
                                               random_state=0)

    gbes.fit(X, y)
    
    return gbes

        
def write_scores_file(dataset, out_file_name, prob, y_true, y_pred, mapqs):

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
        line.append(dataset[i][0])
        line.append(phred)
        line.append(prob[i])
        line.append(1.0 - prob[i])
        line.append(y_true[i])
        line.append(y_pred[i])
        
        data.append(line)
        
    #out.write(dataset[i] + '\t' + str(phred) + '\t' + str(prob[i]) + '\t' + str(1.0 - prob[i]) + '\t' + str(y_true[i]) + '\t' + str(y_pred[i]) + '\n')
    
    df = pd.DataFrame(data)
    df.to_csv(out_file_name, sep='\t', header=False)


def run_ru_boost(X, y):
    
    dt_stump = DecisionTreeClassifier(max_depth=1, min_samples_leaf=1)

    ada_real = RUSBoostClassifier(base_estimator=dt_stump, algorithm="SAMME.R")
    
    p_grid = {"n_estimators":[50, 100, 200, 500], "learning_rate":[1, 0.1, 0.01, 0.001]}
    p_grid = {"n_estimators":[750, 1000, 2000], "learning_rate":[1, 0.1, 0.01, 0.001]}
    
    grid_search = GridSearchCV(ada_real, param_grid=p_grid, cv=5, n_jobs=8, scoring='f1_weighted')
    grid_search.fit(X, y)
    model = grid_search.best_estimator_.fit(X, y)

    n_pos = np.count_nonzero(ada_real.classes_)
    n_neg = len(ada_real.classes_) - n_pos
    
    return model, (n_neg, n_pos) 


def run_balanced_bagging(X, y):
    
    
    params = {"n_estimators": [30, 40, 50],
                  "base_estimator__min_samples_split": [2, 10, 20],
                  "base_estimator__max_depth": [None, 5, 10],
                  "base_estimator__min_samples_leaf": [1, 5],
                  "base_estimator__max_leaf_nodes": [None, 5, 10],
                  "base_estimator__max_features":[2, 3]
                  }


    dt = DecisionTreeClassifier()
    bc = BalancedBaggingClassifier(base_estimator=dt, oob_score=False, random_state=None) #n_estimators=70, random_state=1)

    # Grid Search to determine best parameters
    grid_search = GridSearchCV(estimator=bc, param_grid=params, scoring='f1_weighted', cv=5, n_jobs=8)

    grid_search.fit(X, y)
    model = grid_search.best_estimator_.fit(X, y)
    
    n_pos = np.count_nonzero(bc.classes_)
    n_neg = len(bc.classes_) - n_pos
    
    return model, (n_neg, n_pos)
    
    
def run_adaboost(X, y):
    
    dt_stump = DecisionTreeClassifier(max_depth=1, min_samples_leaf=1)

    ada_real = AdaBoostClassifier(base_estimator=dt_stump, algorithm="SAMME.R")
    
    p_grid = {"n_estimators":[50, 100, 200, 500], "learning_rate":[1, 0.1, 0.01, 0.001]}
    p_grid = {"n_estimators":[750, 1000, 2000], "learning_rate":[1, 0.1, 0.01, 0.001]}
    
    grid_search = GridSearchCV(ada_real, param_grid=p_grid, cv=5, n_jobs=24, scoring='f1_weighted')
    grid_search.fit(X, y)
    model = grid_search.best_estimator_.fit(X, y)
    
    return model
    
    
################## Run Parameters ######################################

n_cores = 8

feats_file = BASE_DIR + '/models/tom_all_feats.tab'

#cols = [MAPQ,INSERT_SIZE,GAP_EXTENSIONS,EDIT_DISTANCE,MISMATCHES,GAP_OPENS,ALIGNMENT_SCORE,SECONDARY_ALIGNMENT_SCORE,PAIR_ORIENTATION,PAIR_ALIGNMENT_TYPE,INTERCEPT,SLOPE,N_LOW_QUAL_BASES,CLASS]
cols = (SECONDARY_ALIGNMENT_SCORE,ALIGNMENT_SCORES_DIFF,ALIGNMENT_SCORE,MAPQ,MATE_ALIGNMENT_SCORE,MISMATCHES,EDIT_DISTANCE,CLASS)

est_type = 'lr'
n_subsets = 300

########################################################################


# Set up possible values of parameters to optimize over
if est_type == "lr":
    p_grid = {"C": [1, 10, 100, 1000]} #{"C": [1, 10, 100, 1000]}
    est = LogisticRegression(class_weight='balanced', solver='liblinear')
else:
    p_grid = {'C': [0.001, 0.1, 1.0, 10]}
    est = SVC(kernel="rbf", class_weight='balanced', gamma='scale')


print "loading data file...."

dataset = pd.read_csv(feats_file, usecols=cols, delimiter="\t").values
read_names = pd.read_csv(feats_file, usecols=(0,), delimiter="\t").values

print "done\n"        

last_col = len(cols) - 1
X = dataset[:,0:last_col]
y = dataset[:,last_col]
y = y.astype(int)

mapqs = X[:, 0]
ir = IR( out_of_bounds = 'clip' )


print "scaling input data...."
min_max_scaler = preprocessing.MinMaxScaler()
Xs = min_max_scaler.fit_transform(X)
print "done\n"


print "splitting input into training and test sets...."
X_train, X_test, y_train, y_test = train_test_split(Xs, y, random_state=0, stratify=y, test_size=0.8, train_size=0.2)
print "done\n"

if est_type == 'RULRS' or est_type == 'RUSVMS':
    X_res, y_res, idxs, model, sample_counts = fit(est, n_subsets, X_train, y_train, random_state=None, under_sampling_factor=1, parameter_grid=p_grid, refit=True)
elif est_type == 'RULR' or est_type == 'RUSVM':
    X_res, y_res, idxs, model, sample_counts = fit(est, n_subsets, X_train, y_train, random_state=None, under_sampling_factor=1, parameter_grid=p_grid, refit=False)
elif est_type == 'DTBG':
    model, sample_counts = run_balanced_bagging(X_train, y_train)
elif est_type == 'RUBST':
    model, sample_counts = run_ru_boost(X_train, y_train)
    
#################################################################


joblib.dump(model, BASE_DIR + '/models/' + est_type + '.pkl')


print "calibrating probabilities...."

if est_type != 'svm':
    train_proba = model.predict_proba(X_train)[:,1]
    test_proba = model.predict_proba(X_test)
    isr_fit = ir.fit(train_proba, y_train)
    test_proba = test_proba[:,1]
    p_cal = ir.transform(test_proba)
else:    
    train_distance = parallel_decision_function(X_train, model, n_cores)
    isr_fit = ir.fit(train_distance, y_train)
    test_distance = parallel_decision_function(X_test, model, n_cores)
    p_cal = ir.transform(test_distance)
    train_proba = ir.transform(train_distance)
    test_proba = p_cal

print "done"


joblib.dump(isr_fit, BASE_DIR + '/models/' + est_type + '_ISO.pkl')

train_preds = parallel_predict(X_train, model, n_cores)

cm_train_preds = confusion_matrix(y_train, train_preds)

TN = cm_train_preds[0,0]
FN = cm_train_preds[1,0]
TP = cm_train_preds[1,1]
FP = cm_train_preds[0,1]

print "CM training set preds:"
print cm_train_preds

print "done\n"

output_file_header = '\tTP\tFP\tTN\tFN\tPrecision\tAvg Precision\tRecall\tF1 Score\tBrier Score\tNo Pos Class\tNo Neg Class\n'

train_precision = precision_score(y_train, train_preds, average=None)
train_avg_precision = average_precision_score(y_train, train_preds, average='weighted') 
train_recall = recall_score(y_train, train_preds, average=None)
train_f1_score = f1_score(y_train, train_preds)
train_brier_score = brier_score_loss(y_train, train_proba)

train_txt = 'Train' + '\t' + str(TP) + '\t' + str(FP) + '\t' + str(TN) + '\t' + str(FN) + str(train_precision)
train_txt = train_txt + '\t' + str(train_avg_precision) + '\t' + str(train_recall) + '\t' + str(train_f1_score)
train_txt = train_txt + '\t' + str(train_brier_score)

n_pos_class = sample_counts[0]
n_neg_class = sample_counts[1]

train_txt = train_txt + '\t' + str(n_pos_class) + '\t' + str(n_neg_class) + '\n'


test_preds_no_cal = parallel_predict(X_test, model, n_cores)

cm_preds_no_cal = confusion_matrix(y_test, test_preds_no_cal)

test_preds_cal = p_cal > 0.5
test_preds_cal = test_preds_cal.astype(int)
cm_preds_cal = confusion_matrix(y_test, test_preds_cal)


print "Training stats:"
print "Confusion matrix:"
print cm_train_preds
print "precision:"
print train_precision
print "Average precision:"
print train_avg_precision
print "recall:"
print train_recall
print "gmean:"
print geometric_mean_score(y_train, train_preds, average=None)
print "gmean (weighted):"
print geometric_mean_score(y_train, train_preds, average='weighted')
print "F1 score:"
print train_f1_score
print "Brier score:"
print brier_score_loss(y_train, train_proba)
print "No pos class:"
print sample_counts[0]
print "No neg class:"
print sample_counts[1]

print ""

TN = cm_preds_no_cal[0,0]
FN = cm_preds_no_cal[1,0]
TP = cm_preds_no_cal[1,1]
FP = cm_preds_no_cal[0,1]

test_precision = precision_score(y_test, test_preds_no_cal, average=None)
test_avg_precision = average_precision_score(y_test, test_preds_no_cal, average='weighted') 
test_recall = recall_score(y_test, test_preds_no_cal, average=None)
test_f1_score = f1_score(y_test, test_preds_no_cal)
test_brier_score = brier_score_loss(y_test, test_proba)

test_txt = 'Test' + '\t' + str(TP) + '\t' + str(FP) + '\t' + str(TN) + '\t' + str(FN) + str(train_precision)
test_txt = test_txt + '\t' + str(train_avg_precision) + '\t' + str(train_recall) + '\t' + str(train_f1_score)
test_txt = test_txt + '\t' + str(train_brier_score)

test_txt = test_txt + '\t\t\n'

print "Test stats (no recal):"
print "Confusion matrix:"
print cm_preds_no_cal
print "precision:"
print test_precision
print "Average precision:"
print test_avg_precision
print "recall:"
print test_recall
print "gmean:"
print geometric_mean_score(y_test, test_preds_no_cal, average=None)
print "gmean (weighted):"
print geometric_mean_score(y_test, test_preds_no_cal, average='weighted')
print "F1 score:"
print test_f1_score
print "Brier score:"
print test_brier_score

print ""

TN = cm_preds_cal[0,0]
FN = cm_preds_cal[1,0]
TP = cm_preds_cal[1,1]
FP = cm_preds_cal[0,1]

test_cal_precision = precision_score(y_test, test_preds_cal, average=None)
test_cal_avg_precision = average_precision_score(y_test, test_preds_cal, average='weighted')
test_cal_recall = recall_score(y_test, test_preds_cal, average=None)
test_cal_f1_score = f1_score(y_test, test_preds_no_cal)
test_cal_brier_score = brier_score_loss(y_test, p_cal)

test_cal_txt = 'Test' + '\t' + str(TP) + '\t' + str(FP) + '\t' + str(TN) + '\t' + str(FN) + str(train_precision)
test_cal_txt = test_cal_txt + '\t' + str(train_avg_precision) + '\t' + str(train_recall) + '\t' + str(train_f1_score)
test_cal_txt = test_cal_txt + '\t' + str(train_brier_score)

test_cal_txt = test_cal_txt + '\t\t\n'


print "Test stats (recal):"
print "Confusion matrix:"
print cm_preds_cal
print "precision:"
print test_cal_precision
print "Average precision:"
print test_cal_avg_precision
print "recall:"
print test_cal_recall
print "gmean:"
print geometric_mean_score(y_test, test_preds_cal, average=None)
print "gmean (weighted):"
print geometric_mean_score(y_test, test_preds_cal, average='weighted')
print "F1 score:"
print test_cal_f1_score
print "Brier score:"
print test_cal_brier_score


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        