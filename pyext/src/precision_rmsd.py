from __future__ import print_function
import os
import numpy
import IMP
import IMP.em
import pyRMSD.RMSDCalculator
from IMP.sampcon.rmsd_calculation import parse_symm_groups_for_pyrmsd


def parse_custom_ranges(ranges_file):
    if not ranges_file:
        return []
    with open(ranges_file) as fh:
        d = {}
        exec(fh.read(), d)
    return d['density_custom_ranges']


def get_particles_from_superposed(
        cluster_conform_i, cluster_conform_0, align, ps, trans,
        symm_groups=None):
    def _to_vector3ds(numpy_array):
        # No need to fit the whole array - we only need 4 non-coplanar points,
        # so 100 should be plenty
        return [IMP.algebra.Vector3D(c) for c in numpy_array[:100]]

    if align:
        calculator_name = "QCP_SERIAL_CALCULATOR"
    else:
        calculator_name = "NOSUP_SERIAL_CALCULATOR"

    conforms = numpy.array([cluster_conform_0, cluster_conform_i])

    if symm_groups is None:
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(
            calculator_name,
            conforms)
    else:
        if calculator_name == 'NOSUP_SERIAL_CALCULATOR':
            # calc_symm_groups are enough without any fitting
            s1 = parse_symm_groups_for_pyrmsd(symm_groups)
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(
                calculator_name,
                fittingCoordsets=conforms,
                calculationCoordsets=conforms,
                calcSymmetryGroups=s1,
                fitSymmetryGroups=[])
        else:
            calculator = pyRMSD.RMSDCalculator.RMSDCalculator(
                calculator_name,
                fittingCoordsets=conforms,
                calcSymmetryGroups=[],
                fitSymmetryGroups=symm_groups)

    check1 = (calculator_name == 'NOSUP_SERIAL_CALCULATOR')
    check2 = not (symm_groups is None)
    if check1 and check2:
        # superposed calc_coords returned if calc-coords
        # specified while creating the calculator
        rmsd, superposed_fit, _calc_fit = calculator.pairwise(
            0, 1, get_superposed_coordinates=True)
    else:
        rmsd, superposed_fit = calculator.pairwise(
            0, 1, get_superposed_coordinates=True)
    # Get transformation from pyRMSD reference on the first call.
    # This is somewhat inefficient (since we are essentially repeating
    # the pyRMSD calculation) but pyRMSD doesn't appear to make its
    # reference orientation available.
    if trans is None:
        trans = IMP.algebra.get_transformation_aligning_first_to_second(
            _to_vector3ds(superposed_fit[0]), _to_vector3ds(cluster_conform_0))

    for particle_index in range(len(superposed_fit[1])):
        # Transform from pyRMSD back to original reference
        IMP.core.XYZ(ps[particle_index]).set_coordinates(
            trans * IMP.algebra.Vector3D(superposed_fit[1][particle_index]))

    return rmsd, ps, trans


class GetModelDensity(object):
    """Compute mean density maps from structures.
    Keeps a dictionary of density maps,
    keys are in the custom ranges. When you call add_subunits_density, it adds
    particle coordinates to the existing density maps.
    """

    def __init__(self, custom_ranges=None, resolution=20.0, voxel=5.0,
                 bead_names=None):
        """Constructor.
        @param list of particles decorated with mass, radius, and XYZ
        @param resolution The MRC resolution of the output map
               (in Angstrom unit)
        @param voxel The voxel size for the output map (lower is slower)
        """

        self.MRCresolution = resolution
        self.voxel = voxel
        self.count_models = 0.0
        self.densities = {}
        self.bead_names = bead_names
        self.custom_ranges = custom_ranges

        # for each custom range get the particle indices that will be
        # added to the density for that custom range
        self.particle_indices_in_custom_ranges = {}

        for density_name in self.custom_ranges:
            self.particle_indices_in_custom_ranges[density_name] = []

        # go through each bead, put it in the appropriate custom range(s)
        for index, beadname in enumerate(self.bead_names):
            for density_name in self.custom_ranges:
                # each domain in the list custom_ranges[density_name]
                for domain in self.custom_ranges[density_name]:
                    if self._is_contained(beadname, domain):
                        self.particle_indices_in_custom_ranges[
                            density_name].append(index)
                        break  # already added particle to this custom range

    def normalize_density(self):
        pass

    def _create_density_from_particles(self, ps, name,
                                       kernel_type='GAUSSIAN'):
        '''Internal function for adding to densities.
        pass XYZR particles with mass and create a density from them.
        kernel type options are GAUSSIAN, BINARIZED_SPHERE, and SPHERE.'''
        dmap = IMP.em.SampledDensityMap(ps, self.MRCresolution, self.voxel)
        dmap.calcRMS()
        dmap.set_was_used(True)

        if name not in self.densities:
            self.densities[name] = dmap
        else:
            bbox1 = IMP.em.get_bounding_box(self.densities[name])
            bbox2 = IMP.em.get_bounding_box(dmap)
            bbox1 += bbox2
            dmap3 = IMP.em.create_density_map(bbox1, self.voxel)
            dmap3.set_was_used(True)
            dmap3.add(dmap)
            dmap3.add(self.densities[name])
            self.densities[name] = dmap3

    def _is_contained(self, bead_name, domain):
        """ domain can be the name of a single protein or a tuple
            (start_residue,end_residue,protein_name)
        bead is a string of type moleculeName_startResidue_endResidue
        """

        (bead_protein, bead_res_start,
         bead_res_end, bead_copy) = bead_name.split("_")

        # protein name and copy number check
        if isinstance(domain, tuple):
            domain_protein = domain[2]
        else:
            domain_protein = domain
        # A period indicates that we have a copy number
        if "." in domain_protein:
            spl = domain_protein.split(".")
            domain_protein = spl[0]
            domain_copy = int(spl[1])
        else:
            domain_copy = bead_copy = -1

        if bead_protein != domain_protein or int(bead_copy) != domain_copy:
            return False

        # residue range check
        if isinstance(domain, tuple):
            bead_residues = set(range(int(bead_res_start),
                                      int(bead_res_end)+1))
            domain_residues = set(range(int(domain[0]),
                                        int(domain[1])+1))
            return not domain_residues.isdisjoint(bead_residues)
        else:
            return True

    def add_subunits_density(self, ps):
        """Add a frame to the densities.
        @param ps List of particles decorated with XYZR and Mass.
        """
        self.count_models += 1.0
        # initialize custom list of particles
        particles_custom_ranges = {}
        for density_name in self.custom_ranges:
            particles_custom_ranges[density_name] = []

        # add each particle to the relevant custom list
        for density_name in self.custom_ranges:
            for particle_index \
                    in self.particle_indices_in_custom_ranges[density_name]:
                particles_custom_ranges[density_name].append(
                    ps[particle_index])

        # finally, add each custom particle list to the density
        for density_name in self.custom_ranges:
            self._create_density_from_particles(
                particles_custom_ranges[density_name], density_name)

    def get_density_keys(self):
        return list(self.densities.keys())

    def get_density(self, name):
        """Get the current density for some component name"""
        if name not in self.densities:
            return None
        else:
            return self.densities[name]

    def write_mrc(self, path=".", file_prefix=""):
        for density_name in self.densities:
            mrc = os.path.join(path, file_prefix + "_" + density_name + ".mrc")
            self.densities[density_name].multiply(1. / self.count_models)
            IMP.em.write_map(
                self.densities[density_name], mrc,
                IMP.em.MRCReaderWriter())
        if len(self.densities) == 1:
            return mrc
        else:
            return os.path.join(path, file_prefix + "_*.mrc")
