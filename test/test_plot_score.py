import unittest
import subprocess
import sys
import os
import utils


TESTDIR = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


class Tests(unittest.TestCase):
    def test_plot_score_one(self):
        """Test plot_score.py with one score"""
        script = utils.get_script(TOPDIR, 'plot_score.py')
        inp_dir = os.path.join(TESTDIR, 'input')
        with utils.temporary_directory() as tmpdir:
            subprocess.check_call(
                [sys.executable, script,
                 os.path.join(inp_dir, 'model_ids_scores.txt'),
                 'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen'],
                cwd=tmpdir)
            os.unlink(os.path.join(
                tmpdir,
                'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen.png'))

    def test_plot_score_all(self):
        """Test plot_score.py with all scores"""
        script = utils.get_script(TOPDIR, 'plot_score.py')
        inp_dir = os.path.join(TESTDIR, 'input')
        expected = [
            'ConnectivityRestraint_Rpb1.png',
            'CrossLinkingMassSpectrometryRestraint_Distance_.png',
            'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen.png',
            'Total_Score.png',
            'ExcludedVolumeSphere_None.png']
        with utils.temporary_directory() as tmpdir:
            subprocess.check_call(
                [sys.executable, script,
                 os.path.join(inp_dir, 'model_ids_scores.txt'), 'all'],
                cwd=tmpdir)
            for e in expected:
                os.unlink(os.path.join(tmpdir, e))

    def test_plot_score_bad(self):
        """Test plot_score.py with bad score"""
        script = utils.get_script(TOPDIR, 'plot_score.py')
        inp_dir = os.path.join(TESTDIR, 'input')
        out = subprocess.check_output(
            [sys.executable, script,
             os.path.join(inp_dir, 'model_ids_scores.txt'), 'garbage'],
            universal_newlines=True)
        self.assertEqual(out.rstrip('\r\n'),
            "garbage is not a valid score parameter. Use 'all' or one of: "
            "Model_index, Run_id, Replica_id, Frame_id, "
            "CrossLinkingMassSpectrometryRestraint_Distance_, "
            "ConnectivityRestraint_Rpb1, "
            "CrossLinkingMassSpectrometryRestraint_Data_Score_Chen, "
            "ExcludedVolumeSphere_None, Total_Score")


if __name__ == '__main__':
    unittest.main()
