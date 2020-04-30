from __future__ import print_function
import sys, os
import math
import numpy as np


def get_scores_distribution(scores, nbins, scorecolumn, hist_ofile):
    
    H, xedge = np.histogram(scores, bins=nbins)
    with open(hist_ofile, 'w+') as f1:
        for i in range(nbins):
            print(xedge[i], H[i], file=f1)
    return

def get_top_scorings_statistics(scores, scorecolumn, systemname):

    list_threshold = []
    [list_threshold.append( int((factor / 10.0) *len(scores))) for factor in range(1, 11)]
    with open("%s.Top_Score_Conv.txt" % systemname, 'w+') as f1:
        print("Getting top scoring models at subsets of size:",list_threshold)

        for t in list_threshold:
            samples = np.array([np.random.choice(
                        scores, t, replace=False).min() for i in range(150)])
            print(t, samples.mean(), samples.std(), file=f1)

def get_scores_distributions_KS_Stats(score_A, score_B, nbins, systemname):
    import scipy as sp
    from scipy.stats import mannwhitneyu, ks_2samp
    d_stat, p_value = ks_2samp(score_A, score_B)

    get_scores_distribution(score_A, nbins, 0, "%s.Score_Hist_A.txt" % systemname)
    get_scores_distribution(score_B, nbins, 0, "%s.Score_Hist_B.txt" % systemname) 
    
    with open("%s.KS_Test.txt" % systemname, 'w+') as f1:
        print(d_stat, p_value, file=f1)
    return d_stat, p_value
