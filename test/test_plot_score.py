import unittest
import subprocess
import sys
import os
import IMP.test


class Tests(IMP.test.TestCase):
    def test_plot_score_one(self):
        """Test plot_score.py with one score"""
        with IMP.test.temporary_directory() as tmpdir:
            subprocess.check_call(
                [sys.executable, '-m', 'IMP.sampcon.plot_score',
                 self.get_input_file_name('model_ids_scores.txt'),
                 'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen'],
                cwd=tmpdir)
            os.unlink(os.path.join(
                tmpdir,
                'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen.png'))

    def test_plot_score_all(self):
        """Test plot_score.py with all scores"""
        expected = [
            'ConnectivityRestraint_Rpb1.png',
            'CrossLinkingMassSpectrometryRestraint_Distance_.png',
            'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen.png',
            'Total_Score.png',
            'ExcludedVolumeSphere_None.png']
        with IMP.test.temporary_directory() as tmpdir:
            subprocess.check_call(
                [sys.executable, '-m', 'IMP.sampcon.plot_score',
                 self.get_input_file_name('model_ids_scores.txt'), 'all'],
                cwd=tmpdir)
            for e in expected:
                os.unlink(os.path.join(tmpdir, e))

    def test_plot_score_bad(self):
        """Test plot_score.py with bad score"""
        out = subprocess.check_output(
            [sys.executable, '-m', 'IMP.sampcon.plot_score',
             self.get_input_file_name('model_ids_scores.txt'), 'garbage'],
             universal_newlines=True)
        self.assertEqual(out.rstrip('\r\n'),
            "garbage is not a valid score parameter. Use 'all' or one of: "
            "Model_index, Run_id, Replica_id, Frame_id, "
            "CrossLinkingMassSpectrometryRestraint_Distance_, "
            "ConnectivityRestraint_Rpb1, "
            "CrossLinkingMassSpectrometryRestraint_Data_Score_Chen, "
            "ExcludedVolumeSphere_None, Total_Score")


if __name__ == '__main__':
    IMP.test.main()
