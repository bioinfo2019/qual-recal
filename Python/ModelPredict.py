#!/path/to/anaconda3/bin/python

import argparse
import numpy as np
import pandas as pd
import joblib
import math

from sklearn import preprocessing
from sklearn.metrics import confusion_matrix
from sklearn.isotonic import IsotonicRegression as IR
from sklearn.metrics import brier_score_loss, f1_score
from sklearn.metrics import recall_score, precision_score
from sklearn.metrics import average_precision_score

feats = {'MAPQ':3,
        'INSERT_SIZE':4,
        'GAP_EXTENSIONS':5,
        'EDIT_DISTANCE':6,
        'MISMATCHES':7,
        'GAP_OPENS':8,
        'ALIGNMENT_SCORE':9,
        'SECONDARY_ALIGNMENT_SCORE':10,
        'READ_COMPLEXITY_SCORE':11,
        'READ_GC_CONTENT':12,
        'PAIR_ORIENTATION':13,
        'PAIR_ALIGNMENT_TYPE':14,
        'INTERCEPT':15,
        'SLOPE':16,
        'R_VALUE':17,
        'DEPTH':18,
        'READ_EIGENVAL_1':19,
        'READ_EIGENVAL_2':20,
        'READ_TRACE':21,
        'READ_DET':22,
        'REF_EIGENVAL_1':23,
        'REF_EIGENVAL_2':24,
        'REF_TRACE':25,
        'REF_DET':26,
        'REF_COMPLEXITY_SCORE':27,
        'REF_GC_CONTENT':28,
        'N_LOW_QUAL_BASES':29,
        'AVG_LOW_QUAL_BASE_SCORE':30,
        'MATE_ALIGNMENT_SCORE':31,
        'ALIGNMENT_SCORES_DIFF':32,
        'CLASS':33}


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
    
    # Run jobs: works150bp
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
        line.append(read_names[i])
        line.append(phred)

        
        data.append(line)

    #out.write(dataset[i] + '\t' + str(phred) + '\t' + str(prob[i]) + '\t' + str(1.0 - prob[i]) + '\t' + str(y_true[i]) + '\t' + str(y_pred[i]) + '\n')
    
    df = pd.DataFrame(data)
    df.to_csv(out_file_name, sep='\t', header=False, index=False)

################## Parse command line args #############################
    
parser = argparse.ArgumentParser(description='Train model and write it to disk')

parser.add_argument('-b','--base-dir', help='Top level folder containing scripts and data', required=True)
parser.add_argument('-m','--ml-model', help='Type of model to train', required=True)
parser.add_argument('-f','--features', help='Comma separated list of features to be used for prediction', required=True)
parser.add_argument('-v','--feats-file', help='TSV file containing feature vectors for prediction', required=True)

args = parser.parse_args()

BASE_DIR = args.base_dir
est_type = args.ml_model
feats_file = BASE_DIR + '/models/' + args.feats_file
feature_list = args.features

cols = tuple([feats[key] for key in feature_list.split(',')])

########################################################################



print("loading data file....")

d = pd.read_csv(feats_file, usecols=cols, delimiter="\t")
#wrong_idxs = pd.read_csv('/hdd1/eliot_files/MyFiles/PhD/models/wrong_idxs.txt').iloc[:,0].to_list()
#fps = d.iloc[wrong_idxs, :]
#fps.to_csv(BASE_DIR + '/models/M82_1200K.wrong_feats.tab', header=False, index=False)
#sys.exit(0)
#d = d[d.MAPQ == 60]
#d.drop(['MAPQ'], axis=1, inplace=True)

read_names = pd.read_csv(feats_file, usecols=(0,), delimiter="\t").values

dataset = d.values

print("done\n")    

last_col = len(d.columns) - 1
X = dataset[:,0:last_col]
y = dataset[:,last_col]
y = y.astype(int)

#mapqs = X[:, 0]
ir = IR()


print("scaling input data....")
min_max_scaler = preprocessing.MinMaxScaler()
Xs = min_max_scaler.fit_transform(X)
print("done\n")

X = 0

print("Predicting and calibrating probabilities...\n")

clf = joblib.load(BASE_DIR + '/models/' + est_type + '.TOM046A.FULL.pkl')
ir = joblib.load(BASE_DIR + '/models/' + est_type + '_ISO.TOM046A.FULL.pkl')

class_pred = parallel_predict(Xs, clf, 8)
#class_pred = clf.predict(Xs)

if est_type != 'RUSVM' and est_type != 'RUSVMS':
    proba_pred = parallel_predict_proba(Xs, clf, 8)[:,1]
    #proba_pred = clf.predict_proba(Xs)[:,1]
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

#wrong = (test_preds_cal != y).nonzero()
#np.savetxt("/hdd1/eliot_files/MyFiles/PhD/models/wrong_idxs.txt", wrong, fmt='%d', delimiter='\n')

print("Test stats (no recal):")
print("Confusion matrix:")
print(cm_preds_no_cal)
print("precision:")
print(precision_score(y, class_pred, average='micro'))
print("Average precision:")
if est_type == 'RUSVM' or est_type == 'RUSVMS':
    print(average_precision_score(y, distance, average='weighted'))
else:
    print(average_precision_score(y, proba_pred, average='weighted'))


print("recall:")
print(recall_score(y, class_pred, average='micro'))
print("gmean:")
#print(geometric_mean_score(y, class_pred, average='weighted'))
print("F1 score:")
print(f1_score(y, class_pred))
print("Brier score:")
print(brier_score_loss(y, proba_pred))


print("")

print("Test stats (recal):")
print("Confusion matrix:")
print(cm_preds_cal)
print("precision:")
print(precision_score(y, test_preds_cal, average='micro'))
print("Average precision:")
print(average_precision_score(y, p_cal, average='weighted'))
print("recall:")
print(recall_score(y, test_preds_cal, average=None))
print("gmean:")
#print(geometric_mean_score(y, test_preds_cal, average=None))
print("gmean (weighted):")
#print(geometric_mean_score(y, test_preds_cal, average='weighted'))
print("F1 score:")
print(f1_score(y, class_pred))
print("Brier score:")
print(brier_score_loss(y, p_cal))

print("\nWriting scores file...")


if est_type == 'RUSVM' or est_type == 'RUSVMS':
    write_scores_file(read_names, BASE_DIR + '/models/scores.tab', proba_pred)
else:
    write_scores_file(read_names, BASE_DIR + '/models/scores.tab', p_cal) 

print("Done!")
