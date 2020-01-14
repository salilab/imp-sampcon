import unittest
import subprocess
import sys
import os
import shutil
import utils


TESTDIR = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def make_models(tmpdir):
    """Get a set of good-scoring models to use as input"""
    script = utils.get_script(TOPDIR, 'select_good_scoring_models.py')

    mod_dir = os.path.join(tmpdir, 'modeling')
    shutil.copytree(os.path.join(TESTDIR, 'modeling'), mod_dir)

    subprocess.check_call(
        [sys.executable, script,
         '-rd', mod_dir, '-rp', 'run',
         '-sl', 'CrossLinkingMassSpectrometryRestraint_Distance_',
         '-pl', 'ConnectivityRestraint_Rpb1',
         'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen',
         'ExcludedVolumeSphere_None', 'Total_Score',
         '-alt', '1.0', '-aut', '1.0', '-mlt', '0.0', '-mut', '15.0',
         '-e'])


class Tests(unittest.TestCase):
    def test_exhaust(self):
        """Test the master sampling exhaustiveness script"""
        with utils.temporary_directory() as tmpdir:
            make_models(tmpdir)
            gsm_dir = os.path.join(tmpdir, 'modeling', 'good_scoring_models')
            script = utils.get_script(
                        TOPDIR, 'Master_Sampling_Exhaustiveness_Analysis.py')
            subprocess.check_call(
                [sys.executable, script, '-n', 'test', '-p', gsm_dir,
                 '-d', os.path.join(TESTDIR, 'input', 'density_ranges.txt'),
                 '-m', 'cpu_omp', '-c', '8', '-a', '-g', '0.5', '-gp'],
                cwd=tmpdir)

            # Check for expected files
            expected = [
                'Distances_Matrix.data.npy', 'Identities_A.txt',
                'Identities_B.txt', 'cluster.0.all.txt',
                'cluster.0.sample_A.txt',
                'cluster.0.sample_B.txt', 'test.ChiSquare.pdf',
                'test.ChiSquare_Grid_Stats.txt', 'test.Cluster_Population.pdf',
                'test.Cluster_Population.txt', 'test.Cluster_Precision.txt',
                'test.KS_Test.txt', 'test.Sampling_Precision_Stats.txt',
                'test.Score_Dist.pdf', 'test.Score_Hist_A.txt',
                'test.Score_Hist_B.txt', 'test.Top_Score_Conv.pdf',
                'test.Top_Score_Conv.txt',
                'cluster.0/cluster_center_model.rmf3',
                'cluster.0/LPD_TestAll.mrc',
                'cluster.0/Sample_A/LPD_TestAll.mrc',
                'cluster.0/Sample_B/LPD_TestAll.mrc']

            for e in expected:
                os.unlink(os.path.join(tmpdir, e))
            os.rmdir(os.path.join(tmpdir, 'cluster.0', 'Sample_A'))
            os.rmdir(os.path.join(tmpdir, 'cluster.0', 'Sample_B'))
            os.rmdir(os.path.join(tmpdir, 'cluster.0'))


if __name__ == '__main__':
    unittest.main()
