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
        subprocess.check_call(
            [sys.executable, script,
             os.path.join(inp_dir, 'model_ids_scores.txt'),
             'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen'])
        os.unlink('CrossLinkingMassSpectrometryRestraint_Data_Score_Chen.png')

    def test_plot_score_all(self):
        """Test plot_score.py with all scores"""
        script = utils.get_script(TOPDIR, 'plot_score.py')
        inp_dir = os.path.join(TESTDIR, 'input')
        subprocess.check_call(
            [sys.executable, script,
             os.path.join(inp_dir, 'model_ids_scores.txt'), 'all'])
        os.unlink('ConnectivityRestraint_Rpb1.png')
        os.unlink('CrossLinkingMassSpectrometryRestraint_Distance_.png')
        os.unlink('CrossLinkingMassSpectrometryRestraint_Data_Score_Chen.png')
        os.unlink('Total_Score.png')
        os.unlink('ExcludedVolumeSphere_None.png')


if __name__ == '__main__':
    unittest.main()
