from __future__ import print_function
import pyRMSD.RMSDCalculator
from pyRMSD.condensedMatrix import CondensedMatrix

import numpy as np
import os
import glob

import IMP
import IMP.atom
import IMP.rmf
import RMF
import multiprocessing as mp
from itertools import combinations


def parse_symm_groups_for_pyrmsd(s):
    output = []
    for grp in s:
        n = len(grp[0])
        if n == 2:
            output.append(grp)
            continue
        for j, k in combinations(np.arange(n), 2):
            sub_grp = np.array(grp)[:, np.array([j, k])]
            output.append(sub_grp.tolist())
    return output


def parse_rmsd_selection(h, selection, resolution=1):
    s0 = None
    for domain_list in selection.values():
        # each element of the dictionary is a list of domains
        # each domain is a tuple like (start,end,protein) or (protein)
        # to add copy number use (start,end,protein.copy_number)
        # or (protein.copy_number)

        for domain in domain_list:

            start_res = int(domain[0])
            end_res = int(domain[1])

            prot_plus_copy = domain[2]

            if "." in prot_plus_copy:
                copy_number = int(prot_plus_copy.split(".")[1])
                prot_name = prot_plus_copy.split(".")[0]

            else:
                copy_number = 0
                prot_name = prot_plus_copy

            s = IMP.atom.Selection(h, resolution=resolution,
                                   molecule=prot_name,
                                   copy_index=copy_number,
                                   residue_indexes=range(start_res, end_res+1))

            if s0:
                s0 |= s
            else:
                s0 = s
    return s0


def get_pdbs_coordinates(path, idfile_A, idfile_B):
    pts = []
    conform = []
    num = 0
    masses = []
    radii = []

    models_name = []

    with open(idfile_A, 'w+') as f1:
        for str_file in sorted(
                glob.glob("%s/sample_A/*.pdb" % path),
                key=lambda x: int(x.split('/')[-1].split('.')[0])):
            print(str_file, num, file=f1)
            models_name.append(str_file)

            m = IMP.Model()
            mh = IMP.atom.read_pdb(str_file, m,
                                   IMP.atom.NonWaterNonHydrogenPDBSelector())
            mps = IMP.core.get_leaves(mh)
            pts = [IMP.core.XYZ(p).get_coordinates() for p in mps]
            if num == 0:
                masses = [IMP.atom.Mass(p).get_mass() for p in mps]
                radii = [IMP.core.XYZR(p).get_radius() for p in mps]
            conform.append(pts)
            pts = []
            num = num + 1

    with open(idfile_B, 'w+') as f2:
        for str_file in sorted(
                glob.glob("%s/sample_B/*.pdb" % path),
                key=lambda x: int(x.split('/')[-1].split('.')[0])):
            print(str_file, num, file=f2)
            models_name.append(str_file)

            m = IMP.Model()
            mh = IMP.atom.read_pdb(str_file, m,
                                   IMP.atom.NonWaterNonHydrogenPDBSelector())
            mps = IMP.core.get_leaves(mh)
            pts = [IMP.core.XYZ(p).get_coordinates() for p in mps]
            conform.append(pts)
            pts = []
            num = num + 1

    return np.array(conform), masses, radii, models_name


def get_rmfs_coordinates(path, idfile_A, idfile_B,
                         subunit_name=None, selection=None, resolution=1):

    conform = []
    num = 0
    masses = []
    radii = []
    ps_names = []

    f1 = open(idfile_A, 'w+')
    f2 = open(idfile_B, 'w+')

    models_name = []

    for sample_name, sample_id_file in zip(['A', 'B'], [f1, f2]):

        for str_file in sorted(
                glob.glob("%s/sample_%s/*.rmf3" % (path, sample_name)),
                key=lambda x: int(x.split('/')[-1].split('.')[0])):
            print(str_file, num, file=sample_id_file)
            models_name.append(str_file)

            m = IMP.Model()
            inf = RMF.open_rmf_file_read_only(str_file)
            h = IMP.rmf.create_hierarchies(inf, m)[0]
            IMP.rmf.load_frame(inf, 0)

            pts = []

            if subunit_name:
                s0 = IMP.atom.Selection(h, resolution=resolution,
                                        molecule=subunit_name)
            elif selection is not None:
                s0 = parse_rmsd_selection(h, selection, resolution)
            else:
                s0 = IMP.atom.Selection(h, resolution=resolution)

            for leaf in s0.get_selected_particles():

                p = IMP.core.XYZR(leaf)
                pts.append([p.get_coordinates()[i] for i in range(3)])

                if num == 0 and sample_name == 'A':
                    masses.append(IMP.atom.Mass(leaf).get_mass())
                    radii.append(p.get_radius())
                    mol_name = \
                        IMP.atom.get_molecule_name(IMP.atom.Hierarchy(leaf))
                    # traverse up the Hierarchy to get copy number of the
                    # molecule that the bead belongs to
                    copy_number = \
                        IMP.atom.get_copy_index(IMP.atom.Hierarchy(leaf))

                    if IMP.atom.Fragment.get_is_setup(leaf):
                        # TODO not tested on non-fragment systems
                        residues_in_bead = \
                            IMP.atom.Fragment(leaf).get_residue_indexes()

                        ps_names.append(
                            mol_name + "_" + str(min(residues_in_bead)) + "_"
                            + str(max(residues_in_bead)) + "_"
                            + str(copy_number))

                    else:
                        residue_in_bead = \
                            str(IMP.atom.Residue(leaf).get_index())

                        ps_names.append(mol_name + "_" + residue_in_bead + "_"
                                        + residue_in_bead + "_"
                                        + str(copy_number))

            conform.append(pts)
            pts = []
            num = num + 1
    f1.close()
    f2.close()

    return ps_names, masses, radii, np.array(conform), models_name


def parse_symmetric_groups_file(symm_groups_file):
    symm_groups = []
    member_to_symm_group = {}
    first_group_member = []
    curr_particle_index_in_group = []

    sgf = open(symm_groups_file, 'r')

    for indx, ln in enumerate(sgf.readlines()):

        symm_groups.append([])  # create new symm group list
        # particle index for new symmetric group
        curr_particle_index_in_group.append(-1)
        fields = ln.strip().split()

        for fld in fields:
            # group that the current protein copy belongs to
            for subfld in fld.split('/'):
                member_to_symm_group[subfld] = indx

        # the first group member is special! We create a symm group list
        # of particles for the first group member
        first_group_member.append(fields[0].split('/'))

    sgf.close()

    return (symm_groups, member_to_symm_group, curr_particle_index_in_group,
            first_group_member)


def get_conforms_per_frame_batch(arg_bundle):
    rmf_file, frames, mod_id_start = arg_bundle[:3]
    resolution, subunit_name, selection, path = arg_bundle[3:]
    from collections import defaultdict
    m = IMP.Model()
    rmf_fh = RMF.open_rmf_file_read_only(os.path.join(path, rmf_file))
    h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
    result = defaultdict(list)
    mod_id = mod_id_start
    for f in frames:
        IMP.rmf.load_frame(rmf_fh, f)
        m.update()
        if subunit_name:
            s0 = IMP.atom.Selection(h, resolution=resolution,
                                    molecule=subunit_name)
        elif selection is not None:
            s0 = parse_rmsd_selection(h, selection, resolution)
        else:
            s0 = IMP.atom.Selection(h, resolution=resolution)

        particles = s0.get_selected_particles()

        # Copy particle coordinates
        for i in range(len(particles)):
            # i is an index over all particles in the system

            leaf = particles[i]
            p = IMP.core.XYZR(leaf)
            pxyz = p.get_coordinates()
            result[mod_id].append(list(pxyz))
        mod_id += 1
    return result


def get_rmfs_coordinates_one_rmf(path, rmf_A, rmf_B,
                                 subunit_name=None,
                                 symm_groups_file=None, selection=None,
                                 resolution=1, n_cores=None):
    '''Modified RMF coordinates function to work with symmetric copies'''
    if n_cores is None:
        n_cores = 1
    # Open RMFs and get total number of models
    rmf_fh = RMF.open_rmf_file_read_only(os.path.join(path, rmf_A))
    n_models = [rmf_fh.get_number_of_frames()]
    rmf_fh = RMF.open_rmf_file_read_only(os.path.join(path, rmf_B))
    n_models.append(rmf_fh.get_number_of_frames())

    masses = []
    radii = []
    ps_names = []

    models_name = []

    # Build hierarchy from the RMF file
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
    IMP.rmf.load_frame(rmf_fh, 0)
    m.update()

    ######
    # Initialize output array
    pts = 0  # number of particles in each model

    # Get selection
    if subunit_name:
        s0 = IMP.atom.Selection(h, resolution=resolution,
                                molecule=subunit_name)
    elif selection is not None:
        s0 = parse_rmsd_selection(h, selection, resolution)
    else:
        s0 = IMP.atom.Selection(h, resolution=resolution)

    # Count particles
    for leaf in s0.get_selected_particles():
        p = IMP.core.XYZR(leaf)
        pts += 1

    # Initialize array
    conform = np.empty([n_models[0]+n_models[1], pts, 3])

    # Initialize the symmetric group particles list, and protein to
    # symmetry group mapping
    if symm_groups_file:
        (symm_groups, group_member_to_symm_group_map,
         curr_particle_index_in_group, first_group_member) = \
                 parse_symmetric_groups_file(symm_groups_file)
    else:
        symm_groups = None

    # add masses, particle names and radii
    rmf_fh = RMF.open_rmf_file_read_only(os.path.join(path, rmf_A))
    h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
    IMP.rmf.load_frame(rmf_fh, 0)
    m.update()
    if subunit_name:
        s0 = IMP.atom.Selection(h, resolution=resolution,
                                molecule=subunit_name)
    elif selection is not None:
        s0 = parse_rmsd_selection(h, selection, resolution)
    else:
        s0 = IMP.atom.Selection(h, resolution=resolution)
    particles = s0.get_selected_particles()

    # Copy particle coordinates
    for i in range(len(particles)):
        # i is an index over all particles in the system

        leaf = particles[i]
        p = IMP.core.XYZR(leaf)
        masses.append(IMP.atom.Mass(leaf).get_mass())
        radii.append(p.get_radius())
        mol_name = IMP.atom.get_molecule_name(
                IMP.atom.Hierarchy(leaf))
        # traverse up the Hierarchy to get copy number of
        # the molecule that the bead belongs to
        copy_number = \
            IMP.atom.get_copy_index(IMP.atom.Hierarchy(leaf))

        # Add to symmetric groups if needed
        if symm_groups_file:

            protein_plus_copy = mol_name+'.'+str(copy_number)
            if protein_plus_copy in group_member_to_symm_group_map:
                # protein copy is in a symmetric group
                group_index = group_member_to_symm_group_map[
                    protein_plus_copy]
                curr_particle_index_in_group[group_index] += 1

                if protein_plus_copy \
                        in first_group_member[group_index]:
                    symm_groups[group_index].append([i])

                else:
                    j = curr_particle_index_in_group[group_index] \
                            % len(symm_groups[group_index])
                    symm_groups[group_index][j].append(i)

        # TODO not tested on non-fragment systems
        if IMP.atom.Fragment.get_is_setup(leaf):
            residues_in_bead = \
                IMP.atom.Fragment(leaf).get_residue_indexes()
            ps_names.append(
                mol_name + "_" + str(min(residues_in_bead)) + "_"
                + str(max(residues_in_bead)) + "_"
                + str(copy_number))
        else:
            residue_in_bead = \
                str(IMP.atom.Residue(leaf).get_index())
            ps_names.append(
                mol_name + "_" + residue_in_bead + "_" +
                residue_in_bead + "_" + str(copy_number))

    mod_id = 0  # index for each model in conform.
    for rmf_file in [rmf_A, rmf_B]:

        models_name.append(rmf_file)
        rmf_fh = RMF.open_rmf_file_read_only(os.path.join(path, rmf_file))
        h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
        n_frames = rmf_fh.get_number_of_frames()
        print("Opening RMF file:", rmf_file, "with",
              n_frames, "frames")

        n_cores = min(n_cores, n_frames)  # to ensure n_per_core > 0
        n_per_core = n_frames // n_cores
        spacing = np.arange(0, n_per_core * n_cores, n_per_core)
        mod_id_starts = spacing + mod_id
        frame_number_starts = spacing
        frame_number_ends = spacing + n_per_core
        frame_number_ends[-1] = n_frames - 1
        frame_lists = []
        for i in range(n_cores):
            a = frame_number_starts[i]
            b = frame_number_ends[i]
            frame_lists.append(np.arange(a, b + 1))
        p = mp.Pool(n_cores)
        args_list = [(rmf_file, frame_lists[i], mod_id_starts[i], resolution,
                     subunit_name, selection, path) for i in range(n_cores)]
        results = p.map(get_conforms_per_frame_batch, args_list)
        for res in results:
            for m_id in res:
                conform[m_id, :, :] = np.array(res[m_id])
        mod_id += n_frames

    if symm_groups_file:
        for grp in symm_groups:
            if len(grp) == 0:
                print("Warning. Symmetry option specified but created "
                      "symmetry group is empty. Cross-check the "
                      "specification of symmetry groups.")

    return ps_names, masses, radii, conform, symm_groups, models_name, n_models


def get_rmsds_matrix(conforms,  mode,  sup,  cores, symm_groups=None):

    if (mode == "cpu_serial" and not sup) or (mode == "cpu_omp" and not sup):
        calculator_name = "NOSUP_OMP_CALCULATOR"

    elif mode == "cpu_omp" and sup:
        calculator_name = "QCP_OMP_CALCULATOR"
        print("we are using QCP_OMP to compute RMSD")
    elif mode == "cuda" and sup:
        calculator_name = "QCP_CUDA_MEM_CALCULATOR"
    else:
        print("Wrong values to pyRMSD ! Please Fix")
        exit()

    if symm_groups:
        print("We have ambiguity.")
        s1 = parse_symm_groups_for_pyrmsd(symm_groups)
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(
            calculator_name,
            fittingCoordsets=conforms,
            calcSymmetryGroups=s1,
            fitSymmetryGroups=s1)

    else:
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(
            calculator_name, conforms)

    # additionally set number of cores for parallel calculator
    if mode == "cpu_omp":
        calculator.setNumberOfOpenMPThreads(int(cores))

    if not symm_groups:
        rmsd = calculator.pairwiseRMSDMatrix()
    else:
        rmsd = []
        for i in range(len(conforms) - 1):
            rmsd += list(calculator.oneVsFollowing(i))
    rmsd_matrix = CondensedMatrix(rmsd)
    inner_data = rmsd_matrix.get_data()
    np.save("Distances_Matrix.data", inner_data)

    return inner_data
