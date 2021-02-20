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


def get_rmfs_coordinates(path, idfile_A, idfile_B, subunit_name):

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
                s0 = IMP.atom.Selection(h, resolution=1, molecule=subunit_name)
            else:
                s0 = IMP.atom.Selection(h, resolution=1)

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
            member_to_symm_group[fld] = indx

        # the first group member is special! We create a symm group list
        # of particles for the first group member
        first_group_member.append(fields[0])

    sgf.close()

    return (symm_groups, member_to_symm_group, curr_particle_index_in_group,
            first_group_member)


def get_rmfs_coordinates_one_rmf(path, rmf_A, rmf_B, subunit_name=None,
                                 symm_groups_file=None):

    '''Modified RMF coordinates function to work with symmetric copies'''

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
        s0 = IMP.atom.Selection(h, resolution=1, molecule=subunit_name)
    else:
        s0 = IMP.atom.Selection(h, resolution=1)

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

    mod_id = 0  # index for each model in conform.

    for rmf_file in [rmf_A, rmf_B]:

        models_name.append(rmf_file)

        rmf_fh = RMF.open_rmf_file_read_only(os.path.join(path, rmf_file))
        h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]

        print("Opening RMF file:", rmf_file, "with",
              rmf_fh.get_number_of_frames(), "frames")

        for f in range(rmf_fh.get_number_of_frames()):

            if f % 100 == 0:
                # pass
                print("  -- Opening frame", f, "of",
                      rmf_fh.get_number_of_frames())

            IMP.rmf.load_frame(rmf_fh, f)

            m.update()

            # Store particle indices and loop over individual protein
            # names for symmetric copies
            if subunit_name:
                s0 = IMP.atom.Selection(h, resolution=1, molecule=subunit_name)
            else:
                s0 = IMP.atom.Selection(h, resolution=1)

            particles = s0.get_selected_particles()

            # Copy particle coordinates
            for i in range(len(particles)):
                # i is an index over all particles in the system

                leaf = particles[i]
                p = IMP.core.XYZR(leaf)
                pxyz = p.get_coordinates()
                conform[mod_id][i][0] = pxyz[0]
                conform[mod_id][i][1] = pxyz[1]
                conform[mod_id][i][2] = pxyz[2]

                # Just for the first model, update the masses and radii
                # and log the particle name in ps_names
                if mod_id == 0 and rmf_file == rmf_A:
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
                                    == first_group_member[group_index]:
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
            mod_id += 1

    return ps_names, masses, radii, conform, symm_groups, models_name, n_models


def get_rmsds_matrix(conforms,  mode,  sup,  cores, symm_groups=None):
    print("Mode:", mode, "Superposition:", sup, "Number of cores:", cores)

    if (mode == "cpu_serial" and not sup) or (mode == "cpu_omp" and not sup):
        calculator_name = "NOSUP_OMP_CALCULATOR"

    elif mode == "cpu_omp" and sup:
        calculator_name = "QCP_OMP_CALCULATOR"

    elif mode == "cuda" and sup:
        calculator_name = "QCP_CUDA_MEM_CALCULATOR"
    else:
        print("Wrong values to pyRMSD ! Please Fix")
        exit()

    if symm_groups:
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(
            calculator_name, fittingCoordsets=conforms,
            calculationCoordsets=conforms, calcSymmetryGroups=symm_groups)
    else:
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(
            calculator_name, conforms)

    # additionally set number of cores for parallel calculator
    if mode == "cpu_omp":
        calculator.setNumberOfOpenMPThreads(int(cores))

    rmsd = calculator.pairwiseRMSDMatrix()
    rmsd_matrix = CondensedMatrix(rmsd)
    inner_data = rmsd_matrix.get_data()
    np.save("Distances_Matrix.data", inner_data)

    return inner_data
