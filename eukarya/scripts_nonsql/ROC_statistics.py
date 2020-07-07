#!/hosts/linuxhome/scarab/eva2/Programs/miniconda3/bin/python
#python3
import sys
import random
import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score, auc

#https://stackoverflow.com/questions/52373318/how-to-compare-roc-auc-scores-of-different-binary-classifiers-and-assess-statist
#https://stats.stackexchange.com/questions/165033/how-to-interpret-95-confidence-interval-for-area-under-curve-of-roc
#https://datascience.stackexchange.com/questions/14369/how-to-bootstrap-the-auc-on-a-data-set-with-50-000-entries/14376

#KEEP IN MIND: we are working with distances. This means 0 is 'good' and 1 is 'bad'
#This also means our values in the roc calculations are reversed (normally 1 is 'good' and 0 'bad')
def roc_auc_dists(y_true, y_score, pos_label = "BioGrid"):
    #reverse the scores to get roc for distances
    max_score = max(y_score)
    #fpr, tpr, _ = roc_curve(y_true, [max_score - x for x in y_score], pos_label = pos_label)
    #roc_auc = auc(fpr, tpr) #trapezoidal rule, again, because of distances
    roc_auc = roc_auc_score(y_true, y_score)
    return roc_auc

#Test statistics for ROC curves, testing different aspects between and within ROC
# 1. Bootstrap for 95% confidence interval for one ROC curve
##Get a confidence interval, calculate different AUC values by random sampling the data
def bootstrap_auc(y_true, y_score, n_boots=1000): #not really a bootstrap but permutation
    #random sample and calculate again auc n_boots times
    auc_values = []
    for i in range(n_boots):
        idx = np.random.randint(len(y_true),size = len(y_true)) #get random index to sample the same size of the true data (with replacement)
        y_true = y_true.ravel()
        y_score = y_score.ravel()
        auc_values += [roc_auc_dists(y_true[idx], y_score[idx])] # [roc_auc]
    return np.percentile(auc_values, (2.5, 97.5))

# 2. Permutation test to test against chance performance
##This function takes classifier and returns the AUC value on
##the unshuffled data, mean auc of the shuffled data and a p-value
def permutation_auc(y_true, y_score, n_perms=1000):
    #shuffle the labels and recalculate the auc values
    #aka the shuffled vs unshuffeled data
    y_true = y_true.ravel()
    y_score = y_score.ravel()
    idx = np.arange(len(y_score))
    auc_values = np.empty(n_perms)
    for i in range(n_perms):
        np.random.shuffle(idx)
        auc_values[i] = roc_auc_dists(y_true[idx], y_score) #auc values for all the randomized labeled
    #compare now to the original
    roc_auco = roc_auc_dists(y_true, y_score) #original
    #the p-value (i.e. probability to observe an AUC value larger
    #than or equal to what you have in the unshuffled data).
    #Not exact due to sample size (1 vs. n_perms) but indication
    #p-value = fraction of AUC_r greater than or equal to your observed AUC.
    return roc_auco, np.mean(auc_values), np.mean(auc_values >= roc_auco)

# 3. Permutation test for difference between different classifyers
def permutation_different_auc(y_true1, y_score1, y_true2, y_score2, n_perms=1000):
    #shuffle between the two classifyers and recalculate the auc values and
    #get the differences before and after
    y_true1 = y_true1.ravel()
    y_score1 = y_score1.ravel()
    y_true2 = y_true2.ravel()
    y_score2 = y_score2.ravel()
    #get differences between the auc of the original data
    auc_diff = []
    auc1 = roc_auc_dists(y_true1, y_score1)
    auc2 = roc_auc_dists(y_true2, y_score2)
    observed_diff = abs(auc1 - auc2)
    #shortest length of the data, because we need same length data
    shortest_len = len(min([y_score1, y_score2], key=len))
    for i in range(n_perms):
        #shuffle the scores, because of different data lengths we want to have a
        #random pick of distances, and then get the same length of the data
        pi1 = np.random.randint(shortest_len, size = shortest_len)
        pi2 = np.random.randint(shortest_len, size = shortest_len)

        #shuffle scores randomly from one or the other
        mask = np.random.randint(2, size=shortest_len) #randomly get out of score 1 or score 2
        p1 = np.where(mask, y_score1[pi1][0:shortest_len], y_score2[pi2][0:shortest_len])
        p2 = np.where(mask, y_score2[pi2][0:shortest_len], y_score1[pi1][0:shortest_len])
        auc1 = roc_auc_dists(y_true1[pi1], p1)
        auc2 = roc_auc_dists(y_true2[pi2], p2)
        auc_diff += [abs(auc1 - auc2)]
        #the p-value (i.e. probability to observe an AUC value larger
        #than or equal to what you have in the unshuffled data).
        #Not exact due to sample size (1 vs. n_perms) but indication
        #p-value = fraction of AUC_r greater than or equal to your observed AUC.
    return observed_diff,  np.mean(auc_diff), np.mean(auc_diff >= observed_diff)
