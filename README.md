\brief Sampling exhaustiveness protocol

[![Build Status](https://github.com/salilab/imp-sampcon/workflows/build/badge.svg?branch=master)](https://github.com/salilab/imp-sampcon/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/salilab/imp-sampcon/branch/master/graph/badge.svg)](https://codecov.io/gh/salilab/imp-sampcon)

This module implements the sampling exhaustiveness test described in
[Viswanath et al, 2017](https://www.ncbi.nlm.nih.gov/pubmed/29211988).
The protocol is primarily designed to work with models generated by
the [Integrative Modeling Platform (IMP)](https://integrativemodeling.org)
(and more specifically the IMP::pmi module), but could probably be adapted
for other systems.

# Dependencies:

[pyRMSD](https://github.com/salilab/pyRMSD) is needed. (This is a fork of the
original pyRMSD - which is no longer maintained - to fix bugs and add
Python 3 support.)

In the Sali lab, pyRMSD is already built, so can be used with
`module load python2/pyrmsd` or `module load python3/pyrmsd`.

# imp_sampcon: sampling exhaustiveness test {#imp_sampcon}

The protocol is typically run using the `imp_sampcon` command line tool:

 - `imp_sampcon show_stat` to show the available fields (e.g. scoring function
   terms, restraint satisfaction) in an IMP::pmi stat file.
 - `imp_sampcon select_good` to select a subset of good-scoring models from
   a set of IMP::pmi trajectories.
 - `imp_sampcon plot_score` to plot the score distributions of the selected
   models.
 - `imp_sampcon exhaust` to analyze the selected models and determine whether
   sampling was exhaustive.

For a full demonstration of the protocol, see its usage in
IMP's [actin modeling tutorial](https://integrativemodeling.org/tutorials/actin/analysis.html).

# Selecting models for exhaustiveness test

[PMI analysis](https://github.com/salilab/PMI_analysis/) is the current way to select models for the exhaustiveness test. It selects models based on MCMC equilibration and HDBSCAN density-based clustering of model scores.

In this module, we also have scripts for an earlier method to select models based on thresholds on score terms.

`select_good`  selects a set of good-scoring models based on user-specified score thresholds. The options are described in detail in the [actin modeling tutorial](https://salilab.org/pdf/Saltzberg_MethodsMolBiol_2019.pdf).

A few things to note:

1. How does one set score thresholds? `select_good` has two modes to help with this.
First, one can just guess the values of lower and upper thresholds of score terms based on individual stat files from sampling and use `select_good` in FILTER mode (i.e. without -e option) to see how many models one can get with the current thresholds. One can also plot these scores using `plot_score` to aid in refining the thresholds.
Once the final thresholds are obtained, one can use `select_good` in EXTRACT mode (i.e. with -e option) to extract the corresponding RMFs.

2. `-alt` and `-aut` (aggregate lower and upper thresholds) are the lower threshold and upper thresholds to be set for most score terms like `Total_Score`. `-mlt` and `-mut` (member lower and upper thresholds) only need to be set for crosslinks when specifying distance thresholds for each type of crosslink.

3. `-sl` (selection list) specifies which scores are selected based on specified thresholds and `-pl`(print list) is the extra set of scores that will be printed for each selected model. If a score term is already in the selection list it need not be specified in the print list. Currently, empty print lists are not allowed.

4. The script assumes the output stat files from sampling are in `$run_directory/$run_prefix*/output` where $run_directory and $run_prefix are specified by `-rd` and `-rp`.

# Running the exhaustiveness test

See the usage of `exhaust` in IMP's [actin modeling tutorial](https://integrativemodeling.org/tutorials/actin/analysis.html).

## Outputs

The output includes

* *Convergence of top score file* : `*.Top_Score_Conv.txt` contains the average and standard deviation of the top score (second and third columns) obtained by random subsets of size 10%, 20% …. 100% of the total number of models (first column).
* *Difference in score distributions* : `*.Score_Hist_*.txt` stores the histogram of score distributions and `*.KS_test.txt` contains the effect size D and p-value from the KS test.
* *Chi-square test* : `rnapol.ChiSquare_Grid_Stats.txt` contains for each clustering threshold (first column), the p-value, Cramer’s V and population of models in the contingency table (second, third and last columns). `rnapol.Sampling_Precision_Stats.txt` contains the value of sampling precision and the p-value, Cramer’s V and population of models at the sampling precision.

PDF files are generated for results of the above tests if the gnuplot option is specified (`-gp`).

* *Clusters* :
A directory is created for each cluster (`cluster.0`, `cluster.1` and so on). The indices of models belonging to each cluster x are listed in `cluster.x.all.txt`, and listed by sample  in `cluster.x.sample_A.txt` and `cluster.x.sample_B.txt`.
Cluster populations are in `rnapol.Cluster_Population.txt`, showing the number of models in samples A (second column) and B (third column) in each cluster (first column).
Cluster precisions are in `rnapol.Cluster_Precision.txt`, with the precision defined by the average RMSD to the cluster centroid.
The individual cluster directories contain representative bead models and localization densities.

## Special input options
Here a few special cases of `exhaust` are mentioned.

### Clustering without sampling precision calculation
Sometimes one would like to simply obtain clusters given a threshold, without going through the sampling precision calculation. For example, this could be because the user has already run `exhaust` once on the same input and found too many clusters at the sampling precision and would like to visualize a smaller number of clusters by clustering at a threshold worse than the sampling precision. In this case, one can use the `-s -ct <threshold>` options to skip the expensive steps of distance matrix generation and sampling precision calculation and directly cluster models at the given threshold.

### Alignment
This is usually not necessary if the sampling uses an EM map to place the complex. One must not align when specifying the subunit or a selected set of particles (`-su` or `-sn` option), since the reference frame will be specified by the fixed particles.

### Doing the clustering, RMSD, and precision calculation on selected subunits
By default, `exhaust` considers all subunits for RMSD and clustering. There are a couple of options to specify one or more subunits in particular.
To select a single subunit `-su` option can be used.
To select multiple subunits or domains of subunits, `-sn` option can be used. Protein domains are specified with start and end residue numbers. Each selection is listed as an entry in a dictionary called `density_custom_ranges` in a text file, like in the following example.
```
density_custom_ranges = {"Rpb4":[(1,200,"Rpb4")],"Rpb7":[("Rpb7")],"Rpb11-Rpb14-subcomplex":[("Rpb11"),("Rpb14")]}
```
If there are multiple protein copies, to specify copy number in the protein name the format `prot.copy_number` must be used. For example:
```
density_custom_ranges = {"Rpb4":[(1,200,"Rpb4.0")],"Rpb7.1":[(1,710,"Rpb7.1")]}
```

Note that one usually does not align (`-a`) when specifying select subunits, since the frame of reference is determined by the fixed subunits.

### Getting localization densities
A density file must be specified. The syntax is identical to the selection text file above (`-sn`). The output will contain MRC files, one per element of the dictionary in the density file.

Voxel (`-v`) and density threshold (`-dt`) input options are related to density generation.

### Symmetry and ambiguity

The protocol can also handle systems with ambiguity and symmetry (e.g. multiple protein copies whose positions can be interchanged), where this information needs to be considered while calculating the RMSD between models. The RMSD between two protein models is the minimum RMSD over permutations of equivalent proteins.

**Note**: You only need to use this option if your protein copies are symmetric, i.e. they are interchangeable in the model. If there are multiple copies but they occupy distinct positions in the model, this option will not be useful and will give the same result as regular `exhaust` without ambiguity option.

#### Example
If a system has 2 copies of protein A and 1 copy of protein B, i.e. the proteins are A.0, A.1,B.0. The RMSD between any pair of models m0 and m1, is the minimum RMSD between `RMSD[m0(A.0,A.1,B.0) , m1(A.0,A.1,B.1)]` and `RMSD[m0(A.0,A.1,B.1), m1(A.1,A.0,B.1]`. Note that the copies of A in m1 were interchanged while calculating the second RMSD.

To implement this, pyRMSD takes an additional argument `symm_groups` which is a list of particle indices of equivalent particles. For the above case for instance, `symm_groups` has one symmetric group with the particle indices of A.0 and A.1. `symm_groups=[[[A.0.b0,A.1.b0],[A.0.b1,A.1.b1],[A.0.b2,A.1.b2]..[A.0.bn,A.1.bn]]]`. Here `A.X.bi` is the index of the i'th bead in protein A.X and the ith beads of the two protein copies are considered equivalent particles.

To generate this list of symmetric groups, one needs to pass an additional file with the ambiguity option to the master exhaust script. The file contains one line per symmetric group, and components of symmetric groups are separated by white space. See also the example in `symminput`.
