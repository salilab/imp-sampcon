from __future__ import print_function
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

import IMP.pmi.analysis

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
parser.add_argument('--extension', '-e', dest="extension",
                    help='extension of the file', choices=['rmf', 'pdb'],
                    default="rmf")
parser.add_argument('--mode', '-m', dest="mode", help='pyRMSD calculator',
                    choices=['cuda', 'cpu_omp', 'cpu_serial'], default="cuda")
parser.add_argument('--cores', '-c', dest="cores", type=int,
                    help='number of cores for RMSD matrix calculations; '
                         'only for  cpu_omp', default=1)
parser.add_argument('--subunit','-su',dest="subunit",help='calculate RMSD/sampling and cluster precision/densities etc over this subunit only',default=None)
parser.add_argument('--align', '-a', dest="align", help='boolean flag to allow superposition of models', default=False, action='store_true')
parser.add_argument('--ambiguity', '-amb', dest="ambiguity", help='boolean flag for multiple protein copies', default=False, action='store_true')
parser.add_argument('--scoreA', '-sa', dest="scoreA", help='name of the file having the good-scoring scores for sample A', default="scoresA.txt")
parser.add_argument('--scoreB', '-sb', dest="scoreB",help='name of the file having the good-scoring scores for sample B', default="scoresB.txt")
parser.add_argument('--rmfA', '-ra', dest="rmf_A", help='RMF file with conformations from Sample A', default=None)
parser.add_argument('--rmfB', '-rb', dest="rmf_B",help='RMF file with conformations from Sample B', default=None)
parser.add_argument('--gridsize', '-g', dest="gridsize", type=float,help='grid size for calculating sampling precision', default=10.0)
parser.add_argument('--skip','-s',dest="skip_sampling_precision",help="This option will bypass the calculation of sampling precision. This option needs to be used with the clustering threhsold option.Otherwise by default, sampling precision is calculated and the clustering threshold is the calculated sampling precision.",default=False,action='store_true')
parser.add_argument('--cluster_threshold','-ct',dest="cluster_threshold",type=float,help='final clustering threshold to visualize clusters. Assumes that the user has previously calculated sampling precision and wants clusters defined at a threshold higher than the sampling precision for ease of analysis (lesser number of clusters).',default=30.0)
parser.add_argument('--voxel', '-v', dest="voxel", type=float,help='voxel size for the localization densities', default=5.0)
parser.add_argument('--density_threshold', '-dt', type=float,dest="density_threshold", help='threshold for localization densities', default=20.0)
parser.add_argument('--density', '-d', dest="density", help='file containing dictionary of density custom ranges', default=None)
parser.add_argument('--gnuplot', '-gp', dest="gnuplot", help="plotting automatically with gnuplot", default=False, action='store_true')
args = parser.parse_args()

idfile_A = "Identities_A.txt"
idfile_B = "Identities_B.txt"

#Step 0: Compute Score convergence
score_A = []
score_B = []

with open(os.path.join(args.path, args.scoreA), 'r') as f:
    for line in f:
        score_A.append(float(line.strip("\n")))

with open(os.path.join(args.path, args.scoreB), 'r') as f:
    for line in f:
        score_B.append(float(line.strip("\n")))

scores = score_A + score_B

# Get the convergence of the best score
get_top_scorings_statistics(scores, 0, args.sysname)

# Check if the two score distributions are similar
get_scores_distributions_KS_Stats(score_A, score_B, 100, args.sysname)

#Step 1: Compute RMSD matrix
if args.extension == "pdb":
    ps_names = [] # bead names are not stored in PDB files
    conforms, masses, radii, models_name = get_pdbs_coordinates(args.path, idfile_A, idfile_B)
else:
    args.extension = "rmf3"
    # If we have a single RMF file, read conformations from that
    if args.rmf_A is not None:
        if args.ambiguity is True:
            ps_names, masses, radii, conforms, symm_groups, models_name, n_models = get_rmfs_coordinates_one_rmf_amb(args.path, args.rmf_A, args.rmf_B, args.subunit)
        else:
        	print("here")
        	ps_names, masses, radii, conforms, models_name, n_models = get_rmfs_coordinates_one_rmf(args.path, args.rmf_A, args.rmf_B, args.subunit)

    # If not, default to the Identities.txt file
    else:
        ps_names, masses, radii, conforms, models_name = get_rmfs_coordinates(args.path, idfile_A, idfile_B, args.subunit)

print("Size of conformation matrix",conforms.shape)

if not args.skip_sampling_precision:
    # get_rmsds_matrix modifies conforms, so save it to a file and restore
    # afterwards (so that we retain the original IMP orientation)
    numpy.save("conforms", conforms)
    inner_data = get_rmsds_matrix(conforms, args.mode, args.align, args.cores)
    print("Size of RMSD matrix (flattened):",inner_data.shape)
    del conforms
    conforms = numpy.load("conforms.npy")
    os.unlink('conforms.npy')
import pyRMSD.RMSDCalculator
from pyRMSD.matrixHandler import MatrixHandler
mHandler = MatrixHandler()
mHandler.loadMatrix("Distances_Matrix.data")

rmsd_matrix = mHandler.getMatrix()
distmat = rmsd_matrix.get_data()

distmat_full = sp.spatial.distance.squareform(distmat)
print("Size of RMSD matrix (unpacked, N x N):",distmat_full.shape)

# Get model lists
if args.rmf_A is not None:
    sampleA_all_models=range(n_models[0])
    sampleB_all_models=range(n_models[0],n_models[1]+n_models[0])
    print(sampleA_all_models,sampleB_all_models)
    total_num_models = n_models[1]+n_models[0]
else:
    sampleA_all_models,sampleB_all_models=get_sample_identity(idfile_A, idfile_B)
    total_num_models=len(sampleA_all_models)+len(sampleB_all_models)
all_models=list(sampleA_all_models)+list(sampleB_all_models)
print("Size of Sample A:",len(sampleA_all_models)," ; Size of Sample B: ",len(sampleB_all_models),"; Total", total_num_models)
    
if not args.skip_sampling_precision:
    
    print("Calculating sampling precision")
    
    # Step 2: Cluster at intervals of grid size to get the sampling precision
    gridSize=args.gridsize

    # Get cutoffs for clustering
    cutoffs_list=get_cutoffs_list(distmat, gridSize)
    print("Clustering at thresholds:",cutoffs_list)

    # Do clustering at each cutoff
    pvals, cvs, percents = get_clusters(cutoffs_list, distmat_full, all_models, total_num_models, sampleA_all_models, sampleB_all_models, args.sysname)

    # Now apply the rule for selecting the right precision based on population of contingency table, pvalue and cramersv
    sampling_precision,pval_converged,cramersv_converged,percent_converged = get_sampling_precision(cutoffs_list, pvals, cvs, percents)
        
    # Output test statistics 
    with open("%s.Sampling_Precision_Stats.txt" % args.sysname, 'w+') as fpv:
        print("The sampling precision is defined as the largest allowed RMSD between the cluster centroid and a ",args.sysname,"model within any cluster in the finest clustering for which each sample contributes models proportionally to its size (considering both significance and magnitude of the difference) and for which a sufficient proportion of all models occur in sufficiently large clusters. The sampling precision for our ",args.sysname," modeling is %.3f" %(sampling_precision)," A.", file=fpv)

        print("Sampling precision, P-value, Cramer's V and percentage of clustered models below:", file=fpv)
        print("%.3f\t%.3f\t%.3f\t%.3f" %(sampling_precision, pval_converged, cramersv_converged, percent_converged), file=fpv)
        print("", file=fpv)
    
    final_clustering_threshold = sampling_precision
    
else:
    final_clustering_threshold = args.cluster_threshold
    
# Perform final clustering at the required precision 
print("Clustering at threshold %.3f" %(final_clustering_threshold))
cluster_centers,cluster_members=precision_cluster(distmat_full, total_num_models, final_clustering_threshold)

ctable,retained_clusters=get_contingency_table(len(cluster_centers),cluster_members,all_models,sampleA_all_models,sampleB_all_models)
print("Contingency table:",ctable)
# Output the number of models in each cluster and each sample 
with open("%s.Cluster_Population.txt" % args.sysname, 'w+') as fcp:
    for rows in range(len(ctable)):
        print(rows, ctable[rows][0], ctable[rows][1], file=fcp)

# Obtain the subunits for which we need to calculate densities
density_custom_ranges = parse_custom_ranges(args.density)

# Output cluster precisions
fpc=open("%s.Cluster_Precision.txt" % args.sysname, 'w+')

# For each cluster, output the models in the cluster
# Also output the densities for the cluster models
for i in range(len(retained_clusters)):
    clus=retained_clusters[i]

    # The cluster centroid is the first conformation.
    # We use this as to align and compute RMSD/precision
    conform_0 = conforms[all_models[cluster_members[clus][0]]]

    # create a directory for the cluster 
    if not os.path.exists("./cluster.%s" %i):
        os.mkdir("./cluster.%s" %i)
        os.mkdir("./cluster.%s/Sample_A/" % i)
        os.mkdir("./cluster.%s/Sample_B/" % i)
    else:
        shutil.rmtree("./cluster.%s" %i)
        os.mkdir("./cluster.%s" %i)
        os.mkdir("./cluster.%s/Sample_A/" % i)
        os.mkdir("./cluster.%s/Sample_B/" % i)       

    # Create densities for all subunits for both sample A and sample B as well as separately.
    gmd1 = GetModelDensity(custom_ranges=density_custom_ranges,resolution=args.density_threshold, voxel=args.voxel, bead_names=ps_names)
    gmd2 = GetModelDensity(custom_ranges=density_custom_ranges,resolution=args.density_threshold, voxel=args.voxel, bead_names=ps_names)
    gmdt = GetModelDensity(custom_ranges=density_custom_ranges,resolution=args.density_threshold, voxel=args.voxel, bead_names=ps_names)

    # Also output the identities of cluster members
    both_file=open('cluster.'+str(i)+'.all.txt','w')
    sampleA_file=open('cluster.'+str(i)+'.sample_A.txt','w')
    sampleB_file=open('cluster.'+str(i)+'.sample_B.txt','w')

    # Add the cluster center model RMF to the cluster directory
    cluster_center_index = cluster_members[clus][0]
    if args.rmf_A is not None:
        cluster_center_model_id = cluster_center_index
        if cluster_center_index < n_models[0]:
            os.system('rmf_slice -q '+args.path+args.rmf_A+ " ./cluster."+str(i)+"/cluster_center_model.rmf --frame "+str(cluster_center_index) )
        else:
            os.system('rmf_slice -q '+args.path+args.rmf_B+ " ./cluster."+str(i)+"/cluster_center_model.rmf --frame "+str(cluster_center_index-n_models[0]) )
    else:
        cluster_center_model_id = all_models[cluster_center_index] # index to Identities file.
        shutil.copy(models_name[cluster_center_model_id],os.path.join("./cluster."+str(i),"cluster_center_model."+args.extension))
  
    # Create a model with just the cluster_member particles
    model = IMP.Model()
    ps = [] # particle list to be updated by each RMF frame
    for pi in range(len(conform_0)):
        p = IMP.Particle(model, "%s" % str(pi))
        IMP.core.XYZ.setup_particle(p, (0,0,0))
        IMP.core.XYZR.setup_particle(p, float(radii[pi]))
        IMP.atom.Mass.setup_particle(p, float(masses[pi]))
        ps.append(p)

    # Obtain cluster precision by obtaining average RMSD of each model to the cluster center
    cluster_precision = 0.0

    # transformation from internal pyRMSD orientation
    trans = None
    # for each model in the cluster
    for mem in cluster_members[clus]:
            
        model_index=all_models[mem]
        
        # get superposition of each model to cluster center and the RMSD between the two
        if args.ambiguity is True:
        	rmsd, superposed_ps, trans = get_particles_from_superposed_amb(conforms[model_index], conform_0, args.align, ps, trans, symm_groups)
        else:
        	rmsd, superposed_ps, trans = get_particles_from_superposed(conforms[model_index], conform_0, args.align, ps, trans)

        model.update() # why not?

        cluster_precision+=rmsd

        # Add the superposed particles to the respective density maps
        gmdt.add_subunits_density(superposed_ps) # total density map
        print(model_index, file=both_file)

        if model_index in sampleA_all_models:
            gmd1.add_subunits_density(superposed_ps) # density map for sample A
            print(model_index, file=sampleA_file)
        else:
            gmd2.add_subunits_density(superposed_ps) # density map for sample B
            print(model_index, file=sampleB_file)
         
    cluster_precision /= float(len(cluster_members[clus]) - 1.0)

    print("Cluster precision (average distance to cluster centroid) of cluster ", str(i), " is %.3f" %(cluster_precision), "A", file=fpc)

    both_file.close()
    sampleA_file.close()
    sampleB_file.close()

    # Finally, output density files for the cluster
    gmdt.write_mrc(path="./cluster.%s" %i,file_prefix = "LPD")
    gmd1.write_mrc(path="./cluster.%s/Sample_A/" % i,file_prefix = "LPD")
    gmd2.write_mrc(path="./cluster.%s/Sample_B/" % i,file_prefix = "LPD")
fpc.close()

# generate plots for the score and structure tests
if args.gnuplot:
    import subprocess
    import glob
    
    thisdir = os.path.dirname(__file__)
    gnuplotdir = os.path.join(thisdir, "..", "..", "utilities",
                              "gnuplot_scripts")
    for filename in sorted(glob.glob(os.path.join(gnuplotdir, "*.plt"))):
        cmd = ['gnuplot', '-e', 'sysname="%s"' % args.sysname, filename]
        print(" ".join(cmd))
        subprocess.check_call(cmd)