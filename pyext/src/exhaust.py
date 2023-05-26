from __future__ import print_function
from IMP import ArgumentParser
import os
import json

__doc__ = "Perform analysis to determine sampling convergence."

############################################################
# Scripts written by Shruthi Viswanath and Ilan E. Chemmama#
#              in Andrej Sali Lab at UCSF.                 #
#   Based on Viswanath, Chemmama et al. Biophys. J. (2017) #
#                                                          #
############################################################


def parse_args():
    parser = ArgumentParser(
        description="First stages of analysis for assessing sampling "
                    "convergence")
    parser.add_argument(
        '--sysname', '-n', dest="sysname",
        help='name of the system', default="")
    parser.add_argument(
        '--path', '-p', dest="path",
        help='path to the good-scoring models', default="./")
    parser.add_argument(
        '--extension', '-e', dest="extension",
        help='extension of the file', choices=['rmf', 'pdb'], default="rmf")
    parser.add_argument(
        '--mode', '-m', dest="mode", help='pyRMSD calculator',
        choices=['cuda', 'cpu_omp', 'cpu_serial'], default="cuda")
    parser.add_argument(
        '--matrix-cores', '-c', dest="cores", type=int,
        help='number of cores for parallel RMSD matrix calculations; '
             'only for cpu_omp', default=1)
    parser.add_argument(
        '--cluster-cores', '-cc', dest="cores2", type=int,
        help='number of cores for clustering at different thresholds'
             ' and parallel IO; only for cpu_omp', default=1)
    parser.add_argument(
        '--resolution', '-r', dest="resolution", type=int,
        help='resolution at which to select proteins in a multiscale system',
        default=1)
    parser.add_argument(
        '--subunit', '-su', dest="subunit",
        help='calculate RMSD/sampling and cluster precision/densities '
             'etc over this subunit only', default=None)
    parser.add_argument(
        '--align', '-a', dest="align",
        help='boolean flag to allow superposition of models',
        default=False, action='store_true')
    parser.add_argument(
        '--ambiguity', '-amb', dest="symmetry_groups",
        help='file containing symmetry groups', default=None)
    parser.add_argument(
        '--scoreA', '-sa', dest="scoreA",
        help='name of the file having the good-scoring scores for sample A',
        default="scoresA.txt")
    parser.add_argument(
        '--scoreB', '-sb', dest="scoreB",
        help='name of the file having the good-scoring scores for sample B',
        default="scoresB.txt")
    parser.add_argument(
        '--rmfA', '-ra', dest="rmf_A",
        help='RMF file with conformations from Sample A', default=None)
    parser.add_argument(
        '--rmfB', '-rb', dest="rmf_B",
        help='RMF file with conformations from Sample B', default=None)
    parser.add_argument(
        '--gridsize', '-g', dest="gridsize", type=float,
        help='grid size for calculating sampling precision', default=10.0)
    parser.add_argument(
        '--skip', '-s', dest="skip_sampling_precision",
        help="This option will bypass the calculation of sampling "
             "precision. This option needs to be used with the clustering "
             "threshold option. Otherwise by default, sampling precision "
             "is calculated and the clustering threshold is the "
             "calculated sampling precision.", default=False,
        action='store_true')
    parser.add_argument(
        '--cluster_threshold', '-ct', dest="cluster_threshold", type=float,
        help='final clustering threshold to visualize clusters. Assumes '
             'that the user has previously calculated sampling precision '
             'and wants clusters defined at a threshold higher than the '
             'sampling precision for ease of analysis (lesser number of '
             'clusters).', default=30.0)
    parser.add_argument(
        '--voxel', '-v', dest="voxel", type=float,
        help='voxel size for the localization densities', default=5.0)
    parser.add_argument(
        '--density_threshold', '-dt', type=float,
        dest="density_threshold",
        help='threshold for localization densities', default=20.0)
    parser.add_argument(
        '--density', '-d', dest="density",
        help='file containing dictionary of density custom ranges',
        default=None)
    parser.add_argument(
        '--gnuplot', '-gp', dest="gnuplot",
        help="plotting automatically with gnuplot", default=False,
        action='store_true')
    parser.add_argument(
        '--selection', '-sn', dest="selection",
        help='file containing dictionary'
        'of selected subunits and residues'
        'for RMSD and clustering calculation'
        "each entry in the dictionary takes the form"
        "'selection name': [(residue_start, residue_end, protein name)",
        default=None)
    parser.add_argument(
        '--prism', '-pr', dest="prism",
        help="Save input files for PrISM", default=False,
        action='store_true')
    return parser.parse_args()


def make_cluster_centroid(infname, frame, outfname, cluster_index,
                          cluster_size, precision, metadata_fname, path):

    import RMF
    # If we have new enough IMP/RMF, do our own RMF slicing with provenance
    if hasattr(RMF.NodeHandle, 'replace_child'):
        print(infname, outfname)
        inr = RMF.open_rmf_file_read_only(infname)
        outr = RMF.create_rmf_file(outfname)
        cpf = RMF.ClusterProvenanceFactory(outr)
        RMF.clone_file_info(inr, outr)
        RMF.clone_hierarchy(inr, outr)
        RMF.clone_static_frame(inr, outr)
        inr.set_current_frame(RMF.FrameID(frame))
        outr.add_frame("f0")
        RMF.clone_loaded_frame(inr, outr)
        rn = outr.get_root_node()
        children = rn.get_children()
        if len(children) == 0:
            return
        rn = children[0]  # Should be the top-level IMP node
        prov = [c for c in rn.get_children() if c.get_type() == RMF.PROVENANCE]
        if not prov:
            return
        prov = prov[0]
        # Add cluster-provenance info
        newp = rn.replace_child(
            prov, "cluster.%d" % cluster_index, RMF.PROVENANCE)
        cp = cpf.get(newp)
        cp.set_members(cluster_size)
        cp.set_precision(precision)
        cp.set_density(os.path.abspath(metadata_fname))
    else:
        # Otherwise, fall back to RMF's command line tool
        import subprocess
        print(infname, frame, outfname)
        subprocess.call(['rmf_slice', path + infname, '-f', str(frame),
                         outfname])


def main():
    args = parse_args()

    import os
    import shutil
    import numpy

    import scipy as sp

    import IMP.sampcon
    from IMP.sampcon import scores_convergence, clustering_rmsd
    from IMP.sampcon import rmsd_calculation, precision_rmsd

    import IMP

    # Write output metadata in JSON format
    metadata = {}
    metadata_fname = "%s.output.json" % args.sysname
    metadata['producer'] = {'name': 'IMP.sampcon',
                            'version': IMP.sampcon.__version__}

    idfile_A = "Identities_A.txt"
    idfile_B = "Identities_B.txt"

    # Step 0: Compute Score convergence
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
    scores_convergence.get_top_scorings_statistics(scores, 0, args.sysname)

    # Check if the two score distributions are similar
    scores_convergence.get_scores_distributions_KS_Stats(
            score_A, score_B, 100, args.sysname)

    # Step 1: Compute RMSD matrix
    if args.extension == "pdb":
        ps_names = []  # bead names are not stored in PDB files
        symm_groups = None
        conforms, masses, radii, models_name = \
            rmsd_calculation.get_pdbs_coordinates(
                args.path, idfile_A, idfile_B)
        metadata['input_frames'] = models_name
    else:
        args.extension = "rmf3"
        if args.selection is not None:
            rmsd_custom_ranges = \
                precision_rmsd.parse_custom_ranges(args.selection)
        else:
            rmsd_custom_ranges = None
        # If we have a single RMF file, read conformations from that
        if args.rmf_A is not None:
            metadata['input_files'] = {'A': args.rmf_A, 'B': args.rmf_B}
            (ps_names, masses, radii, conforms, symm_groups, models_name,
                n_models) = rmsd_calculation.get_rmfs_coordinates_one_rmf(
                     args.path, args.rmf_A, args.rmf_B,
                     args.subunit,
                     args.symmetry_groups,
                     rmsd_custom_ranges,
                     args.resolution,
                     args.cores2)

        # If not, default to the Identities.txt file
        else:
            symm_groups = None
            (ps_names, masses, radii, conforms,
             models_name) = rmsd_calculation.get_rmfs_coordinates(
                     args.path, idfile_A, idfile_B, args.subunit,
                     selection=rmsd_custom_ranges,
                     resolution=args.resolution)
            metadata['input_frames'] = models_name

    print("Size of conformation matrix", conforms.shape)

    if not args.skip_sampling_precision:
        # get_rmsds_matrix modifies conforms, so save it to a file and restore
        # afterwards (so that we retain the original IMP orientation)
        numpy.save("conforms", conforms)
        inner_data = rmsd_calculation.get_rmsds_matrix(
                conforms, args.mode, args.align, args.cores, symm_groups)
        print("Size of RMSD matrix (flattened):", inner_data.shape)
        del conforms
        conforms = numpy.load("conforms.npy")
        os.unlink('conforms.npy')

    from pyRMSD.matrixHandler import MatrixHandler
    mHandler = MatrixHandler()
    mHandler.loadMatrix("Distances_Matrix.data")

    rmsd_matrix = mHandler.getMatrix()
    distmat = rmsd_matrix.get_data()

    distmat_full = sp.spatial.distance.squareform(distmat)
    print("Size of RMSD matrix (unpacked, N x N):", distmat_full.shape)

    # Get model lists
    if args.rmf_A is not None:
        sampleA_all_models = list(range(n_models[0]))
        sampleB_all_models = list(range(n_models[0],
                                        n_models[1] + n_models[0]))
        total_num_models = n_models[1] + n_models[0]
    else:
        (sampleA_all_models,
         sampleB_all_models) = clustering_rmsd.get_sample_identity(
                idfile_A, idfile_B)
        total_num_models = len(sampleA_all_models) + len(sampleB_all_models)
    all_models = list(sampleA_all_models) + list(sampleB_all_models)
    print("Size of Sample A:", len(sampleA_all_models),
          " ; Size of Sample B: ", len(sampleB_all_models),
          "; Total", total_num_models)

    if not args.skip_sampling_precision:

        print("Calculating sampling precision")

        # Step 2: Cluster at intervals of grid size to get the
        # sampling precision
        gridSize = args.gridsize

        # Get cutoffs for clustering
        cutoffs_list = clustering_rmsd.get_cutoffs_list(distmat, gridSize)
        print("Clustering at thresholds:", cutoffs_list)

        # Do clustering at each cutoff
        pvals, cvs, percents = clustering_rmsd.get_clusters(
            cutoffs_list, distmat_full, all_models, total_num_models,
            sampleA_all_models, sampleB_all_models, args.sysname,
            args.cores2)
        metadata['chi_square_grid_stats'] = {
            'cutoffs': list(cutoffs_list),
            'p_value': pvals,
            'cramers_v': cvs,
            'percent_clustered': percents}

        # Now apply the rule for selecting the right precision based
        # on population of contingency table, pvalue and cramersv
        (sampling_precision, pval_converged, cramersv_converged,
         percent_converged) = clustering_rmsd.get_sampling_precision(
                 cutoffs_list, pvals, cvs, percents)

        # Output test statistics
        with open("%s.Sampling_Precision_Stats.txt"
                  % args.sysname, 'w+') as fpv:
            print("The sampling precision is defined as the largest allowed "
                  "RMSD between the cluster centroid and a ", args.sysname,
                  "model within any cluster in the finest clustering for "
                  "which each sample contributes models proportionally to "
                  "its size (considering both significance and magnitude of "
                  "the difference) and for which a sufficient proportion of "
                  "all models occur in sufficiently large clusters. The "
                  "sampling precision for our ", args.sysname,
                  " modeling is %.3f" % (sampling_precision), " A.", file=fpv)

            print("Sampling precision, P-value, Cramer's V and percentage "
                  "of clustered models below:", file=fpv)
            print("%.3f\t%.3f\t%.3f\t%.3f"
                  % (sampling_precision, pval_converged, cramersv_converged,
                     percent_converged), file=fpv)
            print("", file=fpv)

        final_clustering_threshold = sampling_precision
        metadata['precision'] = {
            'sampling_precision': sampling_precision,
            'p_value': pval_converged,
            'cramers_v': cramersv_converged,
            'percent_clustered': percent_converged}

    else:
        final_clustering_threshold = args.cluster_threshold

    metadata['clustering_threshold'] = final_clustering_threshold

    # Perform final clustering at the required precision
    print("Clustering at threshold %.3f" % final_clustering_threshold)
    (cluster_centers, cluster_members) = clustering_rmsd.precision_cluster(
            distmat_full, total_num_models, final_clustering_threshold)

    (ctable, retained_clusters) = clustering_rmsd.get_contingency_table(
            len(cluster_centers), cluster_members, all_models,
            sampleA_all_models, sampleB_all_models)
    print("Contingency table:", ctable)
    # metadata['contingency_table'] = ctable
    # Output the number of models in each cluster and each sample
    with open("%s.Cluster_Population.txt" % args.sysname, 'w+') as fcp:
        for rows in range(len(ctable)):
            print(rows, ctable[rows][0], ctable[rows][1], file=fcp)

    # Obtain the subunits for which we need to calculate densities
    density_custom_ranges = precision_rmsd.parse_custom_ranges(args.density)
    metadata['density_custom_ranges'] = density_custom_ranges

    # Output cluster precisions
    fpc = open("%s.Cluster_Precision.txt" % args.sysname, 'w+')

    metadata['clusters'] = []
    # For each cluster, output the models in the cluster
    # Also output the densities for the cluster models
    for i in range(len(retained_clusters)):
        cmeta = {'name': 'cluster.%d' % i}
        clus = retained_clusters[i]

        # The cluster centroid is the first conformation.
        # We use this as to align and compute RMSD/precision
        conform_0 = conforms[all_models[cluster_members[clus][0]]]

        # create a directory for the cluster
        if not os.path.exists("./cluster.%s" % i):
            os.mkdir("./cluster.%s" % i)
            os.mkdir("./cluster.%s/Sample_A/" % i)
            os.mkdir("./cluster.%s/Sample_B/" % i)
        else:
            shutil.rmtree("./cluster.%s" % i)
            os.mkdir("./cluster.%s" % i)
            os.mkdir("./cluster.%s/Sample_A/" % i)
            os.mkdir("./cluster.%s/Sample_B/" % i)
        # File for saving input to PrISM
        if args.prism:
            prism_file = 'cluster.'+str(i)+'.prism.npz'
            superposed_coords_cluster = []
        # Create densities for all subunits for both sample A and sample B
        # as well as separately.
        gmd1 = precision_rmsd.GetModelDensity(
                custom_ranges=density_custom_ranges,
                resolution=args.density_threshold, voxel=args.voxel,
                bead_names=ps_names)
        gmd2 = precision_rmsd.GetModelDensity(
                custom_ranges=density_custom_ranges,
                resolution=args.density_threshold, voxel=args.voxel,
                bead_names=ps_names)
        gmdt = precision_rmsd.GetModelDensity(
                custom_ranges=density_custom_ranges,
                resolution=args.density_threshold, voxel=args.voxel,
                bead_names=ps_names)

        # Also output the identities of cluster members
        both_file = open('cluster.'+str(i)+'.all.txt', 'w')
        sampleA_file = open('cluster.'+str(i)+'.sample_A.txt', 'w')
        sampleB_file = open('cluster.'+str(i)+'.sample_B.txt', 'w')

        # Create a model with just the cluster_member particles
        model = IMP.Model()
        ps = []  # particle list to be updated by each RMF frame
        for pi in range(len(conform_0)):
            p = IMP.Particle(model, "%s" % str(pi))
            IMP.core.XYZ.setup_particle(p, (0, 0, 0))
            IMP.core.XYZR.setup_particle(p, float(radii[pi]))
            IMP.atom.Mass.setup_particle(p, float(masses[pi]))
            ps.append(p)

        # Obtain cluster precision by obtaining average RMSD of each model
        # to the cluster center
        cluster_precision = 0.0

        cmeta['members'] = {
            'A': [int(x) for x in cluster_members[clus]
                  if x in sampleA_all_models],
            'B': [int(x) for x in cluster_members[clus]
                  if x in sampleB_all_models]}

        # transformation from internal pyRMSD orientation
        trans = None
        # for each model in the cluster
        for mem in cluster_members[clus]:

            model_index = all_models[mem]

            # get superposition of each model to cluster center and the
            # RMSD between the two
            rmsd, superposed_ps, trans = \
                precision_rmsd.get_particles_from_superposed(
                        conforms[model_index], conform_0, args.align,
                        ps, trans, symm_groups)

            model.update()  # why not?

            cluster_precision += rmsd

            # Add the superposed particles to the respective density maps
            gmdt.add_subunits_density(superposed_ps)  # total density map
            print(model_index, file=both_file)

            if model_index in sampleA_all_models:
                # density map for sample A
                gmd1.add_subunits_density(superposed_ps)
                print(model_index, file=sampleA_file)
            else:
                # density map for sample B
                gmd2.add_subunits_density(superposed_ps)
                print(model_index, file=sampleB_file)
            if args.prism:
                superposed_coords = \
                    [IMP.core.XYZ(s_ps).get_coordinates()
                        for s_ps in superposed_ps]
                superposed_coords_cluster.append(
                    numpy.array(superposed_coords))
        if args.prism:
            mass = \
                [IMP.atom.Mass(m_p).get_mass() for m_p in superposed_ps]
            radii = \
                [IMP.core.XYZR(r_p).get_radius() for r_p in superposed_ps]
            numpy.savez(
                prism_file,
                numpy.array(superposed_coords_cluster),
                numpy.array(mass),
                numpy.array(radii),
                numpy.array(ps_names))
        cluster_precision /= float(len(cluster_members[clus]) - 1.0)
        cmeta['precision'] = cluster_precision
        print("Cluster precision (average distance to cluster centroid) "
              "of cluster ", str(i), " is %.3f" % cluster_precision, "A",
              file=fpc)

        both_file.close()
        sampleA_file.close()
        sampleB_file.close()

        # Output density files for the cluster
        cmeta['density'] = gmdt.write_mrc(path="./cluster.%s" % i,
                                          file_prefix="LPD")
        cmeta['densityA'] = gmd1.write_mrc(path="./cluster.%s/Sample_A/" % i,
                                           file_prefix="LPD")
        cmeta['densityB'] = gmd2.write_mrc(path="./cluster.%s/Sample_B/" % i,
                                           file_prefix="LPD")

        # Add the cluster center model RMF to the cluster directory
        cluster_center_index = cluster_members[clus][0]
        if args.rmf_A is not None:
            outfname = os.path.join("cluster.%d" % i,
                                    "cluster_center_model.rmf3")
            cluster_center_model_id = cluster_center_index
            if cluster_center_index < n_models[0]:
                make_cluster_centroid(
                    os.path.join(args.path, args.rmf_A),
                    cluster_center_index,
                    outfname, i, len(cluster_members[clus]),
                    cluster_precision, metadata_fname, args.path)
            else:
                make_cluster_centroid(
                    os.path.join(args.path, args.rmf_B),
                    cluster_center_index - n_models[0],
                    outfname, i, len(cluster_members[clus]),
                    cluster_precision, metadata_fname, args.path)
        else:
            # index to Identities file.
            cluster_center_model_id = all_models[cluster_center_index]
            outfname = os.path.join("cluster.%d" % i,
                                    "cluster_center_model." + args.extension)
            if 'rmf' in args.extension:
                make_cluster_centroid(
                        models_name[cluster_center_model_id], 0, outfname,
                        i, len(cluster_members[clus]),
                        cluster_precision, metadata_fname, args.path)
            else:
                shutil.copy(models_name[cluster_center_model_id], outfname)
        cmeta['centroid'] = {'index': cluster_center_index,
                             'file': outfname}
        metadata['clusters'].append(cmeta)
    with open(metadata_fname, 'w') as jfh:
        json.dump(metadata, jfh)
    # generate plots for the score and structure tests
    if args.gnuplot:
        import subprocess
        import glob

        gnuplotdir = IMP.sampcon.get_data_path("gnuplot_scripts")
        for filename in sorted(glob.glob(os.path.join(gnuplotdir, "*.plt"))):
            cmd = ['gnuplot', '-e', 'sysname="%s"' % args.sysname, filename]
            print(" ".join(cmd))
            subprocess.check_call(cmd)


if __name__ == '__main__':
    main()
