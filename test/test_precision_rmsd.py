import IMP
import IMP.core
import IMP.atom
import IMP.test
import IMP.algebra
import numpy as np

try:
    import pyRMSD
    import pyRMSD.symmTools
    from IMP.sampcon import precision_rmsd as pr
except ImportError:
    pyRMSD = None


class Tests(IMP.test.TestCase):
    """Test get_particles_from_superposed"""

    def setUp(self):
        IMP.test.TestCase.setUp(self)
        IMP.set_log_level(IMP.SILENT)
        IMP.set_check_level(IMP.NONE)

        self.m = IMP.Model()
        self.ps = []
        for i in range(20):
            self.ps.append(self.setup_p('{' + str(i) + '}'))
        protein_skeleton = np.array([[0, 0, 0], [0, 1, 0],
                                     [0, 0, 1], [1, 0, 0],
                                     [1, 1, 1]])
        thetas = [0, np.pi / 6, np.pi / 2, np.pi]
        proteins = [protein_skeleton.copy()]
        # axis is (1, 0, 0)
        for i in range(4):
            t = thetas[i]
            rot_matrix = np.array([[1, 0, 0],
                                   [0, np.cos(t), -np.sin(t)],
                                   [0, np.sin(t), np.cos(t)]])
            proteins.append(np.matmul(rot_matrix,
                            protein_skeleton.T).T + (i + 1) * 3)
        p0, p1, p2, p3, p4 = proteins
        self.proteins = proteins
        self.protein_skeleton = protein_skeleton
        self.coords_trans_rot = [p0, p2]
        t = np.pi / 3
        rot_matrix = np.array([[1, 0, 0],
                               [0, np.cos(t), -np.sin(t)],
                               [0, np.sin(t), np.cos(t)]])
        self.coords_symm_2_trans_rot = [np.concatenate([p0, p2], axis=0),
                                        np.matmul(rot_matrix,
                                        np.concatenate([p2, p0],
                                                       axis=0).T).T]
        self.coords_symm_n = [np.concatenate([p0, p2, p3, p4], axis=0),
                              np.concatenate([p4, p0, p3, p2], axis=0)]

    def setup_p(self, name):
        p = self.m.add_particle(str(name))
        p = IMP.core.XYZR.setup_particle(self.m, p)
        return p

    def extract_coords(self, n):
        coords = []
        for i in range(n):
            coords.append(np.array(self.ps[i].get_coordinates()))
        return coords

    @IMP.test.skipIf(pyRMSD is None, "Requires pyrmsd")
    def test_get_particles_from_superposed_align(self):
        """Check get_particles_from_superposed with align"""
        c1, c2 = self.coords_trans_rot
        rmsd, ps, trans = pr.get_particles_from_superposed(c2, c1,
                                                           True,
                                                           self.ps[:5],
                                                           None,
                                                           None)
        self.assertAlmostEqual(rmsd, 0, delta=1e-5)
        coords = self.extract_coords(5)
        np.testing.assert_allclose(np.array(coords),
                                   self.protein_skeleton,
                                   rtol=0, atol=1e-5)

    @IMP.test.skipIf(pyRMSD is None, "Requires pyrmsd")
    def test_get_particles_from_superposed(self):
        """Check get_particles_from_superposed without align"""
        c1, c2 = self.coords_trans_rot
        rmsd, ps, trans = pr.get_particles_from_superposed(c2, c1,
                                                           False,
                                                           self.ps[:5],
                                                           None,
                                                           None)
        self.assertGreater(rmsd, 0)

    @IMP.test.skipIf(pyRMSD is None
                     or not hasattr(pyRMSD.symmTools,
                                    'symm_groups_validation_new'),
                     "Requires pyrmsd with support for symmetry groups "
                     "with more than 2 elements")
    def test_get_particles_from_superposed_align_symm(self):
        """Check get_particles_from_superposed with align and symm"""
        c1, c2 = self.coords_symm_2_trans_rot
        rmsd, ps, trans = pr.get_particles_from_superposed(c2, c1,
                                                           True,
                                                           self.ps[:10],
                                                           None,
                                                           [[[0, 5],
                                                             [1, 6],
                                                             [2, 7],
                                                             [3, 8],
                                                             [4, 9]]])
        self.assertAlmostEqual(rmsd, 0, delta=1e-5)
        c1, c2 = self.coords_symm_n
        s = [[[0, 5, 10, 15],
              [1, 6, 11, 16],
              [2, 7, 12, 17],
              [3, 8, 13, 18],
              [4, 9, 14, 19]]]
        rmsd, ps, trans = pr.get_particles_from_superposed(c2, c1,
                                                           True,
                                                           self.ps[:20],
                                                           None,
                                                           s)
        self.assertAlmostEqual(rmsd, 0, delta=1e-5)

    @IMP.test.skipIf(pyRMSD is None, "Requires pyrmsd")
    def test_get_particles_from_superposed_symm(self):
        """Check get_particles_from_superposed without align and symm"""
        c1, c2 = self.coords_symm_2_trans_rot
        rmsd, ps, trans = pr.get_particles_from_superposed(c2, c1,
                                                           False,
                                                           self.ps[:10],
                                                           None,
                                                           [[[0, 5],
                                                             [1, 6],
                                                             [2, 7],
                                                             [3, 8],
                                                             [4, 9]]])
        self.assertGreater(rmsd, 0)
        c1, c2 = self.coords_symm_n
        s = [[[0, 5, 10, 15],
              [1, 6, 11, 16],
              [2, 7, 12, 17],
              [3, 8, 13, 18],
              [4, 9, 14, 19]]]
        rmsd, ps, trans = pr.get_particles_from_superposed(c2, c1,
                                                           False,
                                                           self.ps[:20],
                                                           None,
                                                           s)
        self.assertAlmostEqual(rmsd, 0, delta=1e-5)
        coords = self.extract_coords(20)
        np.testing.assert_allclose(np.array(coords)[5:10, :],
                                   self.protein_skeleton, rtol=0,
                                   atol=1e-5)


if __name__ == '__main__':
    IMP.test.main()
