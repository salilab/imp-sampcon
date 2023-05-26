import subprocess
import os
import shutil
import IMP.atom
import IMP.rmf
import RMF
import IMP.test
import sys
import json
from IMP.sampcon import exhaust, select_good

# Skip gnuplot tests on Python 2.7; on our CI systems we get conflicts between
# system gnuplot and conda's ancient Python 2 packages
TEST_GNUPLOT = sys.version_info[0] >= 3


def make_pdbs_from_rmfs(tmpdir):
    for sample in ('sample_A', 'sample_B'):
        sdir = os.path.join(tmpdir, 'modeling', 'good_scoring_models', sample)
        for rmf in os.listdir(sdir):
            if not rmf.endswith('rmf3'):
                continue
            m = IMP.Model()
            r = RMF.open_rmf_file_read_only(os.path.join(sdir, rmf))
            mhs = IMP.rmf.create_hierarchies(r, m)
            m.update()
            IMP.atom.write_pdb(mhs[0], os.path.join(sdir, rmf[:-4] + 'pdb'))


class Tests(IMP.test.TestCase):
    def test_exhaust_help(self):
        """Test exhaust module help"""
        self.check_runnable_python_module("IMP.sampcon.exhaust")

    def make_models(self, tmpdir, make_rmf=False):
        """Get a set of good-scoring models to use as input"""
        mod_dir = os.path.join(tmpdir, 'modeling')
        shutil.copytree(self.get_input_file_name('modeling'), mod_dir)

        self.run_python_module(
            select_good,
            ['-rd', mod_dir, '-rp', 'run',
             '-sl', 'CrossLinkingMassSpectrometryRestraint_Distance_',
             '-pl', 'ConnectivityRestraint_Rpb1',
             'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen',
             'ExcludedVolumeSphere_None', 'Total_Score',
             '-alt', '1.0', '-aut', '1.0', '-mlt', '0.0', '-mut', '15.0',
             '-e'])
        if make_rmf:
            gsm_dir = os.path.join(mod_dir, 'good_scoring_models')

            def read_sample(subdir):
                return sorted(
                    '%s/%s' % (subdir, rmf)
                    for rmf in os.listdir(os.path.join(gsm_dir, subdir))
                    if rmf.endswith('.rmf3'))
            subprocess.check_call(
                ['rmf_cat'] + read_sample('sample_A') + ['A.rmf3'],
                cwd=gsm_dir)
            subprocess.check_call(
                ['rmf_cat'] + read_sample('sample_B') + ['B.rmf3'],
                cwd=gsm_dir)

    def test_exhaust(self):
        """Test the master sampling exhaustiveness script"""
        try:
            import pyRMSD  # noqa: F401
        except ImportError:
            self.skipTest("this test requires the pyRMSD Python module")
        with IMP.test.temporary_working_directory() as tmpdir:
            self.make_models(tmpdir)
            gsm_dir = os.path.join(tmpdir, 'modeling', 'good_scoring_models')
            gnuplot = ['-gp'] if TEST_GNUPLOT else []
            self.run_python_module(
                exhaust,
                ['-n', 'test', '-p', gsm_dir,
                 '-d', self.get_input_file_name('density_ranges.txt'),
                 '-m', 'cpu_omp', '-c', '8', '-a', '-g', '0.5'] + gnuplot)

            if hasattr(RMF.NodeHandle, 'replace_child'):
                r = RMF.open_rmf_file_read_only(
                        os.path.join(tmpdir, 'cluster.0',
                                     'cluster_center_model.rmf3'))
                clpf = RMF.ClusterProvenanceConstFactory(r)
                cpf = RMF.CombineProvenanceConstFactory(r)
                fpf = RMF.FilterProvenanceConstFactory(r)
                rn = r.get_root_node().get_children()[0]
                # Should be one Provenance node
                prov, = [n for n in rn.get_children()
                         if n.get_type() == RMF.PROVENANCE]
                # Top-level provenance should be ClusterProvenance
                self.assertTrue(clpf.get_is(prov))
                cp = clpf.get(prov)
                self.assertEqual(cp.get_members(), 36)
                self.assertAlmostEqual(cp.get_precision(), 0.42, delta=0.01)
                self.assertEqual(cp.get_density(),
                                 os.path.abspath('test.output.json'))
                # Next provenance should be filter, combine
                prov, = prov.get_children()
                self.assertTrue(fpf.get_is(prov))
                prov, = prov.get_children()
                self.assertTrue(cpf.get_is(prov))
                # Make sure that output JSON is valid
                with open('test.output.json') as fh:
                    _ = json.load(fh)

            # Check for expected files
            expected = [
                'Distances_Matrix.data.npy', 'Identities_A.txt',
                'Identities_B.txt', 'cluster.0.all.txt',
                'cluster.0.sample_A.txt', 'cluster.0.sample_B.txt',
                'test.ChiSquare_Grid_Stats.txt',
                'test.Cluster_Population.txt', 'test.Cluster_Precision.txt',
                'test.KS_Test.txt', 'test.Sampling_Precision_Stats.txt',
                'test.Score_Hist_A.txt', 'test.Score_Hist_B.txt',
                'test.Top_Score_Conv.txt',
                'cluster.0/cluster_center_model.rmf3',
                'cluster.0/LPD_TestAll.mrc',
                'cluster.0/Sample_A/LPD_TestAll.mrc',
                'cluster.0/Sample_B/LPD_TestAll.mrc',
                'test.output.json']
            if TEST_GNUPLOT:
                expected.extend([
                    'test.ChiSquare.pdf', 'test.Cluster_Population.pdf',
                    'test.Score_Dist.pdf', 'test.Top_Score_Conv.pdf'])

            for e in expected:
                os.unlink(os.path.join(tmpdir, e))
            os.rmdir(os.path.join(tmpdir, 'cluster.0', 'Sample_A'))
            os.rmdir(os.path.join(tmpdir, 'cluster.0', 'Sample_B'))
            os.rmdir(os.path.join(tmpdir, 'cluster.0'))

    def test_exhaust_rmf_a_b(self):
        "Test the master sampling exhaustiveness script with rmf_A,B options"
        try:
            import pyRMSD  # noqa: F401
        except ImportError:
            self.skipTest("this test requires the pyRMSD Python module")
        with IMP.test.temporary_working_directory() as tmpdir:
            self.make_models(tmpdir, make_rmf=True)
            gsm_dir = os.path.join(tmpdir, 'modeling', 'good_scoring_models')
            gnuplot = ['-gp'] if TEST_GNUPLOT else []
            self.run_python_module(
                exhaust,
                ['-n', 'test', '-p', gsm_dir,
                 '-ra', 'A.rmf3', '-rb', 'B.rmf3',
                 '-d', self.get_input_file_name('density_ranges.txt'),
                 '-m', 'cpu_omp', '-c', '8', '-a', '-g', '0.5'] + gnuplot)

            # Check for expected files
            expected = [
                'Distances_Matrix.data.npy', 'cluster.0.all.txt',
                'cluster.0.sample_A.txt', 'cluster.0.sample_B.txt',
                'test.ChiSquare_Grid_Stats.txt',
                'test.Cluster_Population.txt', 'test.Cluster_Precision.txt',
                'test.KS_Test.txt', 'test.Sampling_Precision_Stats.txt',
                'test.Score_Hist_A.txt', 'test.Score_Hist_B.txt',
                'test.Top_Score_Conv.txt',
                'cluster.0/cluster_center_model.rmf3',
                'cluster.0/LPD_TestAll.mrc',
                'cluster.0/Sample_A/LPD_TestAll.mrc',
                'cluster.0/Sample_B/LPD_TestAll.mrc']
            if TEST_GNUPLOT:
                expected.extend([
                    'test.ChiSquare.pdf', 'test.Cluster_Population.pdf',
                    'test.Score_Dist.pdf', 'test.Top_Score_Conv.pdf'])

            for e in expected:
                os.unlink(os.path.join(tmpdir, e))
            os.rmdir(os.path.join(tmpdir, 'cluster.0', 'Sample_A'))
            os.rmdir(os.path.join(tmpdir, 'cluster.0', 'Sample_B'))
            os.rmdir(os.path.join(tmpdir, 'cluster.0'))

    def test_exhaust_selection_resolution(self):
        """Test the master sampling exhaustiveness script with
        selection and resolution options"""
        try:
            import pyRMSD  # noqa: F401
        except ImportError:
            self.skipTest("this test requires the pyRMSD Python module")
        with IMP.test.temporary_working_directory() as tmpdir:
            self.make_models(tmpdir, make_rmf=True)
            gsm_dir = os.path.join(tmpdir, 'modeling', 'good_scoring_models')
            gnuplot = ['-gp'] if TEST_GNUPLOT else []
            self.run_python_module(
                exhaust,
                ['-n', 'test', '-p', gsm_dir,
                 '-ra', 'A.rmf3', '-rb', 'B.rmf3',
                 '-sn', self.get_input_file_name('selection.txt'),
                 '-r', '1', '-d', self.get_input_file_name('selection.txt'),
                 '-m', 'cpu_omp', '-c', '8', '-g', '0.5'] + gnuplot)

            # Check for expected files
            expected = [
                'Distances_Matrix.data.npy', 'cluster.0.all.txt',
                'cluster.0.sample_A.txt', 'cluster.0.sample_B.txt',
                'test.ChiSquare_Grid_Stats.txt',
                'test.Cluster_Population.txt', 'test.Cluster_Precision.txt',
                'test.KS_Test.txt', 'test.Sampling_Precision_Stats.txt',
                'test.Score_Hist_A.txt', 'test.Score_Hist_B.txt',
                'test.Top_Score_Conv.txt',
                'cluster.0/cluster_center_model.rmf3',
                'cluster.0/LPD_TestAll.mrc',
                'cluster.0/Sample_A/LPD_TestAll.mrc',
                'cluster.0/Sample_B/LPD_TestAll.mrc']
            if TEST_GNUPLOT:
                expected.extend([
                    'test.ChiSquare.pdf', 'test.Cluster_Population.pdf',
                    'test.Score_Dist.pdf', 'test.Top_Score_Conv.pdf'])

            for e in expected:
                os.unlink(os.path.join(tmpdir, e))
            os.rmdir(os.path.join(tmpdir, 'cluster.0', 'Sample_A'))
            os.rmdir(os.path.join(tmpdir, 'cluster.0', 'Sample_B'))
            os.rmdir(os.path.join(tmpdir, 'cluster.0'))

    def test_exhaust_pdb(self):
        """Test the master sampling exhaustiveness script with PDBs"""
        try:
            import pyRMSD  # noqa: F401
        except ImportError:
            self.skipTest("this test requires the pyRMSD Python module")
        with IMP.test.temporary_working_directory() as tmpdir:
            self.make_models(tmpdir)
            make_pdbs_from_rmfs(tmpdir)
            gsm_dir = os.path.join(tmpdir, 'modeling', 'good_scoring_models')
            self.run_python_module(
                exhaust,
                ['-n', 'test', '-p', gsm_dir,
                 '-m', 'cpu_omp', '-c', '8', '-a', '-g', '0.5', '-e', 'pdb'])


if __name__ == '__main__':
    IMP.test.main()
