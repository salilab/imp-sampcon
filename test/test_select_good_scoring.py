import unittest
import subprocess
import sys
import os
import shutil
import RMF
import IMP.test
from IMP.sampcon import select_good


class Tests(IMP.test.TestCase):
    def test_select_good_help(self):
        """Test select_good module help"""
        self.check_runnable_python_module("IMP.sampcon.select_good")

    def test_select_good_scoring_models(self):
        """Test select_good"""
        with IMP.test.temporary_directory() as tmpdir:
            mod_dir = os.path.join(tmpdir, 'modeling')
            shutil.copytree(self.get_input_file_name('modeling'), mod_dir)
            self.run_python_module(select_good,
                ['-rd', mod_dir, '-rp', 'run',
                 '-sl', 'CrossLinkingMassSpectrometryRestraint_Distance_',
                 '-pl', 'ConnectivityRestraint_Rpb1',
                 'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen',
                 'ExcludedVolumeSphere_None', 'Total_Score',
                 '-alt', '1.0', '-aut', '1.0', '-mlt', '0.0', '-mut', '30.0'])
            model_ids = os.path.join(mod_dir, 'filter', 'model_ids_scores.txt')
            with open(model_ids) as fh:
                wc = len(fh.readlines())
            self.assertEqual(wc, 97)
            os.unlink(model_ids)
            os.rmdir(os.path.join(mod_dir, 'filter'))

    def test_select_good_scoring_models_one_run(self):
        """Test select_good with only one run"""
        with IMP.test.temporary_directory() as tmpdir:
            mod_dir = os.path.join(tmpdir, 'modeling')
            shutil.copytree(self.get_input_file_name('modeling'), mod_dir)
            # Keep only run1
            shutil.rmtree(os.path.join(mod_dir, 'run2'))
            self.run_python_module(select_good,
                ['-rd', mod_dir, '-rp', 'run',
                 '-sl', 'CrossLinkingMassSpectrometryRestraint_Distance_',
                 '-pl', 'ConnectivityRestraint_Rpb1',
                 'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen',
                 'ExcludedVolumeSphere_None', 'Total_Score',
                 '-alt', '1.0', '-aut', '1.0', '-mlt', '0.0', '-mut', '30.0'])
            model_ids = os.path.join(mod_dir, 'filter', 'model_ids_scores.txt')
            with open(model_ids) as fh:
                wc = len(fh.readlines())
            self.assertEqual(wc, 49)
            os.unlink(model_ids)
            os.rmdir(os.path.join(mod_dir, 'filter'))

    def test_select_good_scoring_models_extract(self):
        """Test select_good with extract"""
        with IMP.test.temporary_directory() as tmpdir:
            mod_dir = os.path.join(tmpdir, 'modeling')
            shutil.copytree(self.get_input_file_name('modeling'), mod_dir)
            self.run_python_module(select_good,
                ['-rd', mod_dir, '-rp', 'run',
                 '-sl', 'CrossLinkingMassSpectrometryRestraint_Distance_',
                 '-pl', 'ConnectivityRestraint_Rpb1',
                 'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen',
                 'ExcludedVolumeSphere_None', 'Total_Score',
                 '-alt', '1.0', '-aut', '1.0', '-mlt', '0.0', '-mut', '12.0',
                 '-e'])
            gsm_dir = os.path.join(mod_dir, 'good_scoring_models')
            for score, num in (('A', 5), ('B', 4)):
                score_file = os.path.join(gsm_dir, 'scores%s.txt' % score)
                with open(score_file) as fh:
                    wc = len(fh.readlines())
                self.assertEqual(wc, num)
                os.unlink(score_file)
            model_ids = os.path.join(gsm_dir, 'model_ids_scores.txt')
            with open(model_ids) as fh:
                wc = len(fh.readlines())
            self.assertEqual(wc, 10)
            os.unlink(model_ids)

            model_ids = os.path.join(gsm_dir, 'model_sample_ids.txt')
            with open(model_ids) as fh:
                lines = fh.readlines()
            self.assertEqual(len(lines), 9)
            for line in lines:
                num, sample = line.rstrip('\r\n').split()
                rmf = os.path.join(gsm_dir, 'sample_%s' % sample,
                                   '%s.rmf3' % num)
                if hasattr(RMF.NodeHandle, 'replace_child'):
                    r = RMF.open_rmf_file_read_only(rmf)
                    cpf = RMF.CombineProvenanceConstFactory(r)
                    fpf = RMF.FilterProvenanceConstFactory(r)
                    rn = r.get_root_node().get_children()[0]
                    # Should be one Provenance node
                    prov, = [n for n in rn.get_children()
                             if n.get_type() == RMF.PROVENANCE]
                    # Top-level provenance should be FilterProvenance
                    self.assertTrue(fpf.get_is(prov))
                    fp = fpf.get(prov)
                    self.assertEqual(fp.get_method(), "Best scoring")
                    self.assertEqual(fp.get_frames(), 9)
                    # Next provenance should be CombineProvenance
                    prov, = prov.get_children()
                    self.assertTrue(cpf.get_is(prov))
                    cp = cpf.get(prov)
                    self.assertEqual(cp.get_runs(), 2)
                    self.assertEqual(cp.get_frames(), 100)
                os.unlink(rmf)
            os.unlink(model_ids)
            os.rmdir(os.path.join(gsm_dir, 'sample_A'))
            os.rmdir(os.path.join(gsm_dir, 'sample_B'))

            os.rmdir(gsm_dir)

    def test_select_good_scoring_models_extract_one_run(self):
        """Test select_good with extract, one run"""
        with IMP.test.temporary_directory() as tmpdir:
            mod_dir = os.path.join(tmpdir, 'modeling')
            shutil.copytree(self.get_input_file_name('modeling'), mod_dir)
            # Keep only run1
            shutil.rmtree(os.path.join(mod_dir, 'run2'))
            self.run_python_module(select_good,
                ['-rd', mod_dir, '-rp', 'run',
                 '-sl', 'CrossLinkingMassSpectrometryRestraint_Distance_',
                 '-pl', 'ConnectivityRestraint_Rpb1',
                 'CrossLinkingMassSpectrometryRestraint_Data_Score_Chen',
                 'ExcludedVolumeSphere_None', 'Total_Score',
                 '-alt', '1.0', '-aut', '1.0', '-mlt', '0.0', '-mut', '12.0',
                 '-e'])
            gsm_dir = os.path.join(mod_dir, 'good_scoring_models')
            for score, num in (('A', 2), ('B', 3)):
                score_file = os.path.join(gsm_dir, 'scores%s.txt' % score)
                with open(score_file) as fh:
                    wc = len(fh.readlines())
                self.assertEqual(wc, num)
                os.unlink(score_file)
            model_ids = os.path.join(gsm_dir, 'model_ids_scores.txt')
            with open(model_ids) as fh:
                wc = len(fh.readlines())
            self.assertEqual(wc, 6)
            os.unlink(model_ids)

            model_ids = os.path.join(gsm_dir, 'model_sample_ids.txt')
            with open(model_ids) as fh:
                lines = fh.readlines()
            self.assertEqual(len(lines), 5)
            for line in lines:
                num, sample = line.rstrip('\r\n').split()
                os.unlink(os.path.join(gsm_dir, 'sample_%s' % sample,
                                       '%s.rmf3' % num))
            os.unlink(model_ids)
            os.rmdir(os.path.join(gsm_dir, 'sample_A'))
            os.rmdir(os.path.join(gsm_dir, 'sample_B'))

            os.rmdir(gsm_dir)


if __name__ == '__main__':
    IMP.test.main()
