import unittest
import subprocess
import sys
import os
import IMP.test

class Tests(IMP.test.TestCase):
    def test_show_stat(self):
        """Test show_stat.py"""
        stat = self.get_input_file_name(os.path.join('modeling', 'run1',
            'output', 'stat.0.out'))
        o = subprocess.check_output(
            [sys.executable, '-m', 'IMP.sampcon.show_stat', stat],
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
    IMP.test.main()
