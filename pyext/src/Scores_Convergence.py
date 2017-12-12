import sys, os
import math
import numpy as np


def get_scores_distribution(scores, nbins, scorecolumn, hist_ofile):
    
    H, xedge = np.histogram(scores, bins=nbins)
    f1=open(hist_ofile, 'w+')
    
    for i in range(nbins):
        print >>f1, xedge[i], H[i]
    return

def get_top_scorings_statistics(scores, scorecolumn, systemname):

    list_threshold = []
    [list_threshold.append( int((factor / 10.0) *len(scores))) for factor in range(1, 11)]
    f1=open("%s.TS.txt" % systemname, 'w+')
    print list_threshold

    for t in list_threshold:
        samples = np.array([np.random.choice(scores, t, replace=False).min() for i in range(150)])
        print >>f1, t, samples.mean(), samples.std()

    return

def get_scores_distributions_KS_Stats(score_A, score_B, nbins, systemname):
    import scipy as sp
    from scipy.stats import mannwhitneyu, ks_2samp
    d_stat, p_value = ks_2samp(score_A, score_B)

    get_scores_distribution(score_A, nbins, 0, "%s.HA.txt" % systemname)
    get_scores_distribution(score_B, nbins, 0, "%s.HB.txt" % systemname) 
    
    f1=open("%s.KS.txt" % systemname, 'w+')
    print >>f1, d_stat, p_value
    return d_stat, p_value

