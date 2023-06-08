import sys
if sys.version_info[0] >= 3:
    from io import StringIO as StdoutIO
else:
    from io import BytesIO as StdoutIO
import os
import IMP.test
from IMP.sampcon import show_stat


class Tests(IMP.test.TestCase):
    def test_show_stat_help(self):
        """Test show_stat module help"""
        self.check_runnable_python_module("IMP.sampcon.show_stat")

    def test_show_stat(self):
        """Test show_stat"""
        stat = self.get_input_file_name(
            os.path.join('modeling', 'run1', 'output', 'stat.0.out'))
        old_stdout = sys.stdout
        try:
            sys.stdout = StdoutIO()
            _ = self.run_python_module(show_stat, [stat])
            o = sys.stdout.getvalue()
        finally:
            sys.stdout = old_stdout
        fields = frozenset(o.split("\n"))
        expected_keys = frozenset([
            'ConnectivityRestraint_Score_Rpb1',
            'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen',
            'ExcludedVolumeSphere_Score',
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
