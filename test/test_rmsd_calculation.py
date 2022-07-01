import IMP
import IMP.core
import IMP.atom
import IMP.test
import IMP.rmf
import RMF
import os

try:
    import pyRMSD
    from pyRMSD.matrixHandler import MatrixHandler
    from IMP.sampcon import rmsd_calculation
except ImportError:
    pyRMSD = None

try:
    import scipy as sp
    from scipy import spatial  # noqa
except ImportError:
    sp = None


def get_particles(m, input_rmf):
    # Read input fragment

    inf = RMF.open_rmf_file_read_only(input_rmf)
    mh = IMP.rmf.create_hierarchies(inf, m)[0]
    IMP.rmf.load_frame(inf, 0)

    pts = []
    s0 = IMP.atom.Selection(mh, resolution=1)
    for leaf in s0.get_selected_particles():
        pts.append(IMP.core.XYZR(leaf))

    return pts


def get_particles_swapped(m, input_rmf, subunit1, subunit2):
    # Read input fragment
    inf = RMF.open_rmf_file_read_only(input_rmf)
    mh = IMP.rmf.create_hierarchies(inf, m)[0]
    IMP.rmf.load_frame(inf, 0)

    pts = []
    s1 = IMP.atom.Selection(mh, resolution=1, molecule=subunit1)
    s2 = IMP.atom.Selection(mh, resolution=1, molecule=subunit2)

    for leaf in s1.get_selected_particles():
        pts.append(IMP.core.XYZR(leaf))
    for leaf in s2.get_selected_particles():
        pts.append(IMP.core.XYZR(leaf))

    return pts


class Tests(IMP.test.TestCase):
    """Test ambiguous RMSD calculation"""

    def setUp(self):
        IMP.test.TestCase.setUp(self)
        IMP.set_log_level(IMP.SILENT)
        IMP.set_check_level(IMP.NONE)

        self.m = IMP.Model()
        self.symmetry = self.get_input_file_name(
            "ambiguity.txt")
        self.symmetry_complex1 = self.get_input_file_name(
            "ambiguitycomplex1.txt")
        self.symmetry_complex2 = self.get_input_file_name(
            "ambiguitycomplex2.txt")

        self.pts1 = get_particles(
            self.m,
            self.get_input_file_name("SampledA.rmf3"))

        self.pts2 = get_particles(
            self.m,
            self.get_input_file_name("SampledC.rmf3"))

        self.vec1 = []
        self.vec2 = []

        for xyz in self.pts1:
            self.vec1.append(xyz.get_coordinates())
        for xyz in self.pts2:
            self.vec2.append(xyz.get_coordinates())

        self.pts1sw = get_particles_swapped(
            self.m,
            self.get_input_file_name("SampledA.rmf3"),
            "ProteinA.3", "ProteinA.1")

        self.pts2sw = get_particles_swapped(
            self.m,
            self.get_input_file_name("SampledC.rmf3"),
            "ProteinA.3", "ProteinA.1")

        self.vec1sw = []
        self.vec2sw = []
        for xyz in self.pts1sw:
            self.vec1sw.append(xyz.get_coordinates())
        for xyz in self.pts2sw:
            self.vec2sw.append(xyz.get_coordinates())

    def get_rmsd_matrix(self, align, symmetry, name1='SampledA.rmf3', name2='SampledC.rmf3'):
        (ps, masses, radii,
         conforms, symm_groups, models_name,
         n_models) = rmsd_calculation.get_rmfs_coordinates_one_rmf(
             "./",
             self.get_input_file_name(name1),
             self.get_input_file_name(name2),
             None,
             symmetry,
             None,
             1)

        inner_data = rmsd_calculation.get_rmsds_matrix(  # noqa
            conforms, 'cpu_omp', align, 2, symm_groups)
        del conforms

        mHandler = MatrixHandler()
        mHandler.loadMatrix("Distances_Matrix.data")

        rmsd_matrix = mHandler.getMatrix()
        distmat = rmsd_matrix.get_data()

        distmat_full = sp.spatial.distance.squareform(distmat)
        return distmat_full

    @IMP.test.skipIf(pyRMSD is None, "Requires pyrmsd")
    @IMP.test.skipIf(sp is None, "Requires scipy")
    def test_rmsd_with_neither_alignment_nor_ambiguity(self):
        """Check for rmsd with neither alignment nor ambiguity"""

        rmsd_noali = self.get_rmsd_matrix(False, False)
        rmsd = IMP.atom.get_rmsd(self.pts1, self.pts2)
        os.unlink("./Distances_Matrix.data.npy")
        self.assertAlmostEqual(rmsd, rmsd_noali[0][2], delta=1e-3)

    @IMP.test.skipIf(pyRMSD is None, "Requires pyrmsd")
    def test_rmsd_with_alignment_and_no_ambiguity(self):
        """Check for rmsd with alignment but no ambiguity"""
        rmsd_ali = self.get_rmsd_matrix(True, False)
        align_t = IMP.algebra.get_transformation_aligning_first_to_second(
            self.vec1, self.vec2)

        rmsd = IMP.atom.get_rmsd_transforming_first(
            align_t, self.pts1, self.pts2)
        os.unlink("./Distances_Matrix.data.npy")
        self.assertAlmostEqual(rmsd, rmsd_ali[0][2], delta=1e-3)

    @IMP.test.skipIf(pyRMSD is None, "Requires pyrmsd")
    def test_rmsd_with_alignment_and_with_ambiguity(self):
        """Check for rmsd with alignment and with ambiguity"""
        rmsd_ali_amb = self.get_rmsd_matrix(True, self.symmetry)

        align_t = IMP.algebra.get_transformation_aligning_first_to_second(
            self.vec1, self.vec2)

        rmsd_00 = IMP.atom.get_rmsd_transforming_first(
            align_t, self.pts1, self.pts2)

        align_t = IMP.algebra.get_transformation_aligning_first_to_second(
            self.vec1, self.vec2sw)
        rmsd_01 = IMP.atom.get_rmsd_transforming_first(
            align_t, self.pts1, self.pts2sw)

        align_t = IMP.algebra.get_transformation_aligning_first_to_second(
            self.vec1sw, self.vec2)
        rmsd_10 = IMP.atom.get_rmsd_transforming_first(
            align_t, self.pts1sw, self.pts2)

        align_t = IMP.algebra.get_transformation_aligning_first_to_second(
            self.vec1sw, self.vec2sw)
        rmsd_11 = IMP.atom.get_rmsd_transforming_first(
            align_t, self.pts1sw, self.pts2sw)

        rmsd = min(rmsd_00, rmsd_01, rmsd_10, rmsd_11)
        os.unlink("./Distances_Matrix.data.npy")
        self.assertAlmostEqual(rmsd, rmsd_ali_amb[0][2], delta=1e-3)

    @IMP.test.skipIf(pyRMSD is None, "Requires pyrmsd")
    def test_rmsd_with_symm_group_complexes_no_alignment(self):
        """Check for rmsd with complexes specified versus the original setup"""
        rmsd_complex = self.get_rmsd_matrix(False, self.symmetry_complex1, "SampledAcomplex.rmf3", "SampledCcomplex.rmf3")[0][1]
        rmsd_original = self.get_rmsd_matrix(False, self.symmetry_complex2, "SampledAcomplex.rmf3", "SampledCcomplex.rmf3")[0][1]
        self.assertAlmostEqual(rmsd_complex, rmsd_original, delta=1e-4)
        (ps, masses, radii,
         conforms, symm_groups, models_name,
         n_models) = rmsd_calculation.get_rmfs_coordinates_one_rmf(
             "./",
             self.get_input_file_name("SampledAcomplex.rmf3"),
             self.get_input_file_name("SampledCcomplex.rmf3"),
             None,
             self.symmetry_complex1,
             None,
             1)
        self.assertEqual(len(symm_groups), 1)
        self.assertEqual(len(symm_groups[0]), 443)


if __name__ == '__main__':
    IMP.test.main()
