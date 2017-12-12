import os,sys, shutil
import random
import numpy
import math
import pickle

import scipy as sp
from scipy import spatial
import scipy.stats

from Scores_Convergence import *
from Clustering_RMSD import *
from RMSD_Calculation import *
from Precision_RMSD import *

import IMP
import IMP.rmf
import RMF

import argparse


###########################################################
#Scripts written by Shruthi Viswanath and Ilan E. Chemmama#
#             in Andrej Sali Lab at UCSF.                 #
#  Based on Viswanath, Chemmama et al. Biophys. J. (2017) #
#                                                         #
###########################################################


parser = argparse.ArgumentParser(description="First stages of analysis for assessing sampling convergence")
parser.add_argument('--sysname', '-n', dest="sysname", help='name of the system', default="")
parser.add_argument('--path', '-p', dest="path", help='path to the good-scoring models', default="./")
parser.add_argument('--extension', '-e', dest="extension", help='extension of the file; rmf or pdb', default="rmf")
parser.add_argument('--mode', '-m', dest="mode", help='cuda, cpu_omp, or cpu_serial', default="cuda")
parser.add_argument('--cores', '-c', dest="cores", help='number of cores for RMSD matrix calculations; only for  cpu_omp', default="1")
parser.add_argument('--align', '-a', dest="align", help='boolean flag to allow superposition of models', default=False, action='store_true')
parser.add_argument('--scoreA', '-sa', dest="scoreA", help='name of the file having the good-scoring scores for sample A', default="Scores_A.txt")
parser.add_argument('--scoreB', '-sb', dest="scoreB",help='name of the file having the good-scoring scores for sample B', default="Scores_B.txt")
parser.add_argument('--voxel', '-v', dest="voxel", help='voxel size for the localization densities', default=5.0)
parser.add_argument('--threshold', '-t', dest="threshold", help='threshold for localization densities', default=20.0)
parser.add_argument('--density', '-d', dest="density", help='dictionary of density custom ranges', default=None)
parser.add_argument('--gridsize', '-g', dest="gridsize", help='clustering grid size', default=10.0)
parser.add_argument('--gnuplot', '-gp', dest="gnuplot", help="plotting automatically with gnuplot", default=False, action='store_true')
args = parser.parse_args()

idfile_A = "Identities_A.txt"
idfile_B = "Identities_B.txt"
'''
#Step 0: Compute Score convergence
score_A = []
score_B = []

with open(args.path + args.scoreA, 'r') as f:
    for line in f:
        score_A.append(float(line.strip("\n")))

with open(args.path + args.scoreB, 'r') as f:
    for line in f:
        score_B.append(float(line.strip("\n")))

scores = score_A + score_B
get_top_scorings_statistics(scores, 0, args.sysname)
get_scores_distributions_KS_Stats(score_A, score_B, 100, args.sysname)
'''

#Step 1: Compute RMSD matrix
if args.extension == "pdb":
    conforms, masses, models_name = get_pdbs_coordinates(args.path, idfile_A, idfile_B)
else:
    ps_names, masses, radii, conforms, models_name = get_rmfs_coordinates(args.path, idfile_A, idfile_B)
print conforms.shape

inner_data = get_rmsds_matrix(conforms, args.mode, args.align, args.cores)
print inner_data.shape

import pyRMSD.RMSDCalculator
from pyRMSD.matrixHandler import MatrixHandler
mHandler = MatrixHandler()
mHandler.loadMatrix("Distances_Matrix.data")

rmsd_matrix = mHandler.getMatrix()
distmat = rmsd_matrix.get_data()

distmat_full = sp.spatial.distance.squareform(distmat)
print distmat_full.shape

# Step 2. Clustering starts here:
gridSize=args.gridsize

#GetModel Lists
run1_all_models,run2_all_models=get_run_identity(idfile_A, idfile_B)
total_num_models=len(run1_all_models)+len(run2_all_models)
all_models=run1_all_models+run2_all_models

print len(run1_all_models), len(run2_all_models), total_num_models

#GetCutOffs
cutoffs_list=get_cutoffs_list(distmat, gridSize)
print cutoffs_list

#Do Clustering on a Grid
pvals, cvs, percents = get_clusters(cutoffs_list, distmat_full, all_models, total_num_models, run1_all_models, run2_all_models, args.sysname)

# Now apply the rule for selecting the right precision and the pvalue/cramersv
sampling_precision = get_sampling_precision(cutoffs_list, pvals, cvs, percents)

# Redo the clustering at the required precision of sampling convergence
#TODO could have saved time by storing this beforehand
cluster_centers,cluster_members=precision_cluster(distmat_full, total_num_models, sampling_precision)

# We need to know which clusters to ignore and which to focus on!
ctable,retained_clusters=get_contingency_table(len(cluster_centers),cluster_members,all_models,run1_all_models,run2_all_models)
print ctable
(pval,cramersv)=test_sampling_convergence(ctable,total_num_models)
percent_explained= percent_ensemble_explained(ctable,total_num_models)


fpv=open("%s.PV.txt" % args.sysname, 'w+')
print >>fpv, sampling_precision, pval, cramersv, percent_explained

fcp=open("%s.CP.txt" % args.sysname, 'w+')
for rows in range(len(ctable)):
    print >>fcp, rows, ctable[rows][0], ctable[rows][1]

# Output models to files
fl = open(args.path + args.density, 'r')
density_custom_ranges= fl.readlines()[0].strip()
exec(density_custom_ranges)
fl.close()

fpc=open("%s.PC.txt" % args.sysname, 'w+')

for i in range(len(retained_clusters)):
    clus=retained_clusters[i]

    # create a directory for the cluster 
    if not os.path.exists("./cluster.%s" %i):
        os.mkdir("./cluster.%s" %i)
        os.mkdir("./cluster.%s/Sample_1/" % i)
        os.mkdir("./cluster.%s/Sample_2/" % i)
    else:
        shutil.rmtree("./cluster.%s" %i)
        os.mkdir("./cluster.%s" %i)
        os.mkdir("./cluster.%s/Sample_1/" % i)
        os.mkdir("./cluster.%s/Sample_2/" % i)       
    
    # Now that we have cluster precision we can create densities with the same resolution as cluster precision
    gmd1 = GetModelDensity(custom_ranges=density_custom_ranges,resolution=args.threshold, voxel=args.voxel, molnames=ps_names)
    gmd2 = GetModelDensity(custom_ranges=density_custom_ranges,resolution=args.threshold, voxel=args.voxel, molnames=ps_names)
    gmdt = GetModelDensity(custom_ranges=density_custom_ranges,resolution=args.threshold, voxel=args.voxel, molnames=ps_names)
    
    # Add the superposed particles to the respective density maps
    # Also output the identities of cluster members
    both_file=open('cluster.'+str(i)+'.all.txt','w')
    run1_file=open('cluster.'+str(i)+'.run1.txt','w')
    run2_file=open('cluster.'+str(i)+'.run2.txt','w')
    

    cluster_precision = 0.0
    conform_0 = conforms[all_models[cluster_members[clus][0]]]

    for mem in cluster_members[clus]:
            
        model_index=all_models[mem]
        rmsd, model, superposed_ps = get_particles_from_superposed(conforms[model_index], conform_0, masses, radii, args.align)        
        cluster_precision+=rmsd

        gmdt.add_subunits_density(superposed_ps)
        print >>both_file,model_index

        if model_index in run1_all_models:
            gmd1.add_subunits_density(superposed_ps)
            print >>run1_file, model_index
        else:
            gmd2.add_subunits_density(superposed_ps)
            print >>run2_file, model_index

    cluster_precision /= float(len(cluster_members[clus]) - 1.0)

    print ""
    print >> fpc, "Cluster precision of cluster ", str(i), " is ", cluster_precision, "A"
    print ""
            
    both_file.close()
    run1_file.close()
    run2_file.close()

    # Finally, output densities for the cluster
    gmdt.write_mrc(path="./cluster.%s" %i, file_prefix = "LPD")
    gmd1.write_mrc(path="./cluster.%s/Sample_1/" % i, file_prefix = "LPD")
    gmd2.write_mrc(path="./cluster.%s/Sample_2/" % i, file_prefix = "LPD")
 
if args.gnuplot:
    from os import system
    import glob
    
    for filename in sorted(glob.glob("./gnuplot_scripts/*.plt")):
        system('gnuplot -c %s %s' % (filename, args.sysname))
