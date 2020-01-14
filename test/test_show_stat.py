import unittest
import subprocess
import sys
import os
import utils

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


class Tests(unittest.TestCase):
    def test_show_stat(self):
        """Test show_stat.py"""
        script = utils.get_script(TOPDIR, 'show_stat.py')
        o = subprocess.check_output(
            [sys.executable, script,
             os.path.join(TESTDIR, 'modeling', 'run1', 'output', 'stat.0.out')],
            universal_newlines=True)
        fields = frozenset(o.split("\n"))
        expected_keys = frozenset([
            'ConnectivityRestraint_Rpb1',
            'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen',
            'ExcludedVolumeSphere_None',
            'MonteCarlo_Nframe',
            'MonteCarlo_Temperature',
            'ReplicaExchange_CurrentTemp',
            'Total_Score',
            'rmf_file',
            'rmf_frame_index'
            ])
        self.assertLess(expected_keys, fields)


if __name__ == '__main__':
    unittest.main()
