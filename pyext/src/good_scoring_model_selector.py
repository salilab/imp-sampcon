from __future__ import print_function, division
import RMF
import subprocess
import os
import shutil
import random
import glob
import operator


# If we have new enough IMP/RMF, do our own RMF slicing with provenance
if hasattr(RMF.NodeHandle, 'replace_child'):
    def rmf_slice(infile, frameid, outfile, num_runs, total_num_frames,
                  num_good_scoring):
        inr = RMF.open_rmf_file_read_only(infile)
        outr = RMF.create_rmf_file(outfile)
        cpf = RMF.CombineProvenanceFactory(outr)
        fpf = RMF.FilterProvenanceFactory(outr)
        RMF.clone_file_info(inr, outr)
        RMF.clone_hierarchy(inr, outr)
        RMF.clone_static_frame(inr, outr)
        inr.set_current_frame(RMF.FrameID(frameid))
        outr.add_frame("f0")
        RMF.clone_loaded_frame(inr, outr)
        rn = outr.get_root_node()
        children = rn.get_children()
        if len(children) == 0:
            return
        rn = children[0]  # Should be the top-level IMP node
        prov = [c for c in rn.get_children() if c.get_type() == RMF.PROVENANCE]
        if not prov:
            return
        prov = prov[0]
        # Add combine-provenance info
        newp = rn.replace_child(prov, "combine", RMF.PROVENANCE)
        cp = cpf.get(newp)
        cp.set_frames(total_num_frames)
        cp.set_runs(num_runs)
        # Add filter-provenance info
        newp = rn.replace_child(newp, "filter", RMF.PROVENANCE)
        fp = fpf.get(newp)
        fp.set_frames(num_good_scoring)
        fp.set_method("Best scoring")
        # todo: put in a more appropriate value
        fp.set_threshold(0.)

# Otherwise, fall back to the rmf_slice command line tool
else:
    def rmf_slice(infile, frameid, outfile, num_runs, total_num_frames,
                  num_good_scoring):
        FNULL = open(os.devnull, 'w')
        subprocess.call(['rmf_slice', infile, "-f", str(frameid), outfile],
                        stdout=FNULL, stderr=subprocess.STDOUT)


class GoodScoringModelSelector(object):
    # Authors: Shruthi Viswanath

    ''' Select good-scoring models based on scores and/or data satisfaction.
    Exrtact the corresponding RMFs and put them in a separate directory
    '''

    def __init__(self, run_directory, run_prefix):
        """Constructor.
        @param run_directory the directory containing subdirectories of runs
        @param run_prefix the prefix for each run directory. For e.g. if the
               subdirectories are modeling_run1, modeling_run2, etc. the
               prefix is modeling_run
        """
        self.run_dir = run_directory
        self.run_prefix = run_prefix

        # list with each member as a tuple (run id,replica id,frame id)
        # corresponding to good-scoring models
        self.all_good_scoring_models = []

    def _all_run_dirs(self):
        """Yield (pathname, runid) for all run directories (unsorted)"""
        for x in os.listdir(self.run_dir):
            if x.startswith(self.run_prefix):
                runid = x[len(self.run_prefix):]
                fullpath = os.path.join(self.run_dir, x)
                if os.path.isdir(fullpath) and runid.isdigit():
                    yield (fullpath, runid)

    def _get_subfields_for_criteria(
            self, field_headers, selection_keywords_list,
            printing_keywords_list):
        '''Given the list of keywords, get all the stat file entries
           corresponding to each keyword.'''

        # list of dicts corresponding to field indices for each keyword
        selection_fields = [{} for kw in selection_keywords_list]
        # for ambiguous crosslink distances, it will store all restraints
        # corresponding to the ambivalence in one dict.

        # just a placeholder
        printing_fields = [-1 for j in range(len(printing_keywords_list))]

        for fh_index in field_headers:

            for ki, kw in enumerate(selection_keywords_list):
                # need exact name of field unless it is a xlink distance
                if kw == field_headers[fh_index]:
                    selection_fields[ki][kw] = fh_index

                elif kw in field_headers[fh_index] and \
                        field_headers[fh_index].startswith(
                                "CrossLinkingMassSpectrometry"
                                "Restraint_Distance_"):
                    # handle ambiguous restraints
                    (prot1, res1, prot2, res2) = \
                        field_headers[fh_index].split("|")[3:7]
                    prot1 = prot1.split('.')[0]
                    prot2 = prot2.split('.')[0]

                    if (prot1, res1, prot2, res2) in selection_fields[ki]:
                        selection_fields[ki][
                            (prot1, res1, prot2, res2)].append(fh_index)
                    else:
                        # list of indices corresponding to all
                        # combinations of protein copies
                        selection_fields[ki][(prot1, res1, prot2, res2)] = \
                            [fh_index]

            for ki, kw in enumerate(printing_keywords_list):
                if kw == field_headers[fh_index]:
                    printing_fields[ki] = fh_index

        return selection_fields, printing_fields

    def _get_crosslink_satisfaction(
            self, crosslink_distance_values,
            crosslink_percentage_lower_threshold,
            crosslink_percentage_upper_threshold,
            xlink_distance_lower_threshold, xlink_distance_upper_threshold):
        '''For crosslinks, we want models with atleast x% (e.g. 90%) or more
        crosslink satisfaction. A crosslink is satisfied if the distance
        is between the lower and upper distance thresholds
        @param crosslink_distance_values values of distances in the
               current model
        @param crosslink_percentage_lower_threshold atleast x% of crosslinks
               should be within the below distance thresholds
        @param crosslink_percentage_upper_threshold atmost x% of crosslinks
               should be within the below distance thresholds
               (usually 100%: dummy parameter)
        @param xlink_distance_lower_threshold a crosslink should be at least
               this distance apart (usually 0) to be considered satisfied
        @param xlink_distance_upper_threshold a crosslink should be at most
               this distance apart to be considered satisfied
        '''
        satisfied_xlinks = 0.0
        for d in crosslink_distance_values:
            if d >= xlink_distance_lower_threshold \
                    and d <= xlink_distance_upper_threshold:
                satisfied_xlinks += 1.0

        percent_satisfied = \
            satisfied_xlinks/float(len(crosslink_distance_values))

        if percent_satisfied >= crosslink_percentage_lower_threshold \
                and percent_satisfied <= crosslink_percentage_upper_threshold:
            return percent_satisfied, True
        else:
            return percent_satisfied, False

    def _get_score_satisfaction(self, score, lower_threshold, upper_threshold):
        ''' Check if the score is within the thresholds
        '''
        if score <= upper_threshold and score >= lower_threshold:
            return True
        return False

    def _extract_models_from_trajectories(self, output_dir, num_runs,
                                          total_num_frames, sampleA_set):
        '''Given the list of all good-scoring model indices, extract
           their frames and store them ordered by the list index.
           Store the models in two sample subdirectories.'''
        sampleA_dir = os.path.join(output_dir, "sample_A")
        sampleB_dir = os.path.join(output_dir, "sample_B")
        os.mkdir(sampleA_dir)
        os.mkdir(sampleB_dir)

        num_gsm = sum(1 for e in self.all_good_scoring_models)
        print("Extracting", num_gsm, "good scoring models.")
        model_num = 1

        for i, gsm in enumerate(self.all_good_scoring_models):
            if model_num % (num_gsm/10) == 0:
                print(str(model_num / (num_gsm/10)*10)+"% Complete")
            model_num += 1

            (runid, replicaid, frameid) = gsm

            trajfile = os.path.join(self.run_dir, self.run_prefix+runid,
                                    'output', 'rmfs', replicaid+'.rmf3')

            sample_dir = sampleA_dir if i in sampleA_set else sampleB_dir
            rmf_slice(trajfile, frameid,
                      os.path.join(sample_dir, str(i)+'.rmf3'),
                      num_runs, total_num_frames,
                      len(self.all_good_scoring_models))

    def get_good_scoring_models(
            self, selection_keywords_list=[], printing_keywords_list=[],
            aggregate_lower_thresholds=[], aggregate_upper_thresholds=[],
            member_lower_thresholds=[], member_upper_thresholds=[],
            extract=False):
        ''' Loops over all stat files in the run directory and populates
        the list of good-scoring models.
        @param selection_keywords_list is the list of keywords in the PMI
               stat file that need to be checked for each datatype/score
               in the criteria list
        @param printing_keywords_list is the list of keywords in the PMI
               stat file whose values needs to be printed for selected models
        @param aggregate_lower_thresholds The list of lower bounds on the
               values corresponding to fields in the criteria_list.
               Aggregates are used for terms like % of crosslink satisfaction
               and thresholds of score terms
        @param aggregate_upper_thresholds The list of upper bounds on the
               values corresponding to fields in the criteria_list.
               Aggregates are used for terms like % of crosslink satisfaction
               and thresholds of score terms
        @param member_lower_thresholds The list of lower bounds for values
               of subcomponents of an aggregate term. E.g. for crosslink
               satisfaction the thresholds are on distances for each
               individual crosslink. For score terms this can be ignored
               since thresholds are mentioned in the aggregate fields.
        @param member_upper_thresholds The list of upper bounds for values
               of subcomponents of an aggregate term. E.g. for crosslink
               satisfaction the thresholds are on distances for each
               individual crosslink. For score terms this can be ignored
               since thresholds are mentioned in the aggregate fields.
        '''
        if extract:
            output_dir = os.path.join(self.run_dir, "good_scoring_models")
        else:
            output_dir = os.path.join(self.run_dir, "filter")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir, ignore_errors=True)
        os.mkdir(output_dir)

        outf = open(os.path.join(output_dir, "model_ids_scores.txt"), 'w')
        # Header line first
        print(' '.join(["Model_index", "Run_id", "Replica_id", "Frame_id"]
                       + selection_keywords_list + printing_keywords_list),
              file=outf)

        num_runs = 0
        total_num_frames = 0

        for each_run_dir, runid in sorted(self._all_run_dirs(),
                                          key=operator.itemgetter(1)):
            num_runs += 1

            print("Analyzing", runid)

            for each_replica_stat_file in sorted(
                    glob.glob(os.path.join(each_run_dir, "output",
                                           "stat.*.out")),
                    key=lambda x: int(x.strip('.out').split('.')[-1])):

                replicaid = each_replica_stat_file.strip(".out").split(".")[-1]

                rsf = open(each_replica_stat_file, 'r')

                for line_index, each_model_line in enumerate(rsf.readlines()):
                    # for each model in the current replica

                    if line_index == 0:
                        field_headers = eval(each_model_line.strip())
                        fields_for_selection, fields_for_printing = \
                            self._get_subfields_for_criteria(
                                field_headers, selection_keywords_list,
                                printing_keywords_list)
                        continue

                    frameid = line_index-1

                    dat = eval(each_model_line.strip())

                    total_num_frames += 1
                    model_satisfies = False
                    selection_criteria_values = []

                    for si, score_type in enumerate(selection_keywords_list):
                        if "crosslink" in score_type.lower() \
                                and "distance" in score_type.lower():

                            crosslink_distance_values = []

                            for xltype in fields_for_selection[si]:
                                crosslink_distance_values.append(
                                    min([float(dat[j]) for j in
                                         fields_for_selection[si][xltype]]))

                            satisfied_percent, model_satisfies = \
                                self._get_crosslink_satisfaction(
                                    crosslink_distance_values,
                                    aggregate_lower_thresholds[si],
                                    aggregate_upper_thresholds[si],
                                    member_lower_thresholds[si],
                                    member_upper_thresholds[si])

                            selection_criteria_values.append(satisfied_percent)

                        else:
                            score_value = float(
                                dat[fields_for_selection[si][score_type]])

                            model_satisfies = self._get_score_satisfaction(
                                score_value, aggregate_lower_thresholds[si],
                                aggregate_upper_thresholds[si])
                            selection_criteria_values.append(score_value)

                        if not model_satisfies:
                            break

                    if model_satisfies:

                        # Now get the printing criteria
                        printing_criteria_values = []
                        for si, score_type in enumerate(
                                printing_keywords_list):
                            if fields_for_printing[si] < 0:
                                raise KeyError(
                                    "Bad stat file key '%s': use `imp_sampcon "
                                    "show_stat <stat file>` to see all "
                                    "acceptable keys" % score_type)
                            score_value = float(dat[fields_for_printing[si]])
                            printing_criteria_values.append(score_value)

                        self.all_good_scoring_models.append(
                            (runid, replicaid, frameid))

                        # Print out the scores finally

                        print(' '.join(
                            [str(x) for x in
                             [len(self.all_good_scoring_models) - 1,
                              runid, replicaid, frameid]] +
                            ["%.2f" % s for s in selection_criteria_values] +
                            ["%.2f" % s for s in printing_criteria_values]),
                            file=outf)

                rsf.close()
        outf.close()

        if extract:
            sampA, sampB = self._split_good_scoring_models_into_two_subsets(
                num_runs,
                split_type="divide_by_run_ids" if num_runs > 1 else "random")
            self._extract_models_from_trajectories(
                output_dir, num_runs, total_num_frames, frozenset(sampA))
            return sampA, sampB

    def _split_good_scoring_models_into_two_subsets(
            self, num_runs, split_type="divide_by_run_ids"):
        ''' Get the list of good scoring models and split them into two
            samples. Return the two subsets.

        @param split_type how to split good scoring models into two samples.
               Current options are:
               (a) divide_by_run_ids : where the list of runids is divided
                   into 2. e.g. if we have runs from 1-50, good scoring
                   models from runs 1-25 is sample A and those from
                   runs 26-50 is sample B.
               (b) random : split the set of good scoring models into
                   two subsets at random.
        '''
        sampleA_indices = []
        sampleB_indices = []

        if split_type == "divide_by_run_ids":  # split based on run ids

            half_num_runs = num_runs/2
            for i, gsm in enumerate(self.all_good_scoring_models):
                if int(gsm[0]) <= half_num_runs:
                    sampleA_indices.append(i)
                else:
                    sampleB_indices.append(i)

        elif split_type == "random":
            sampleA_indices = random.sample(
                range(len(self.all_good_scoring_models)),
                len(self.all_good_scoring_models)//2)
            sampleB_indices = [
                i for i in range(len(self.all_good_scoring_models))
                if i not in sampleA_indices]

        # write model and sample IDs to a file
        with open(os.path.join(self.run_dir, 'good_scoring_models',
                               'model_sample_ids.txt'), 'w') as f:
            for i in sampleA_indices:
                print(i, "A", file=f)
            for i in sampleB_indices:
                print(i, "B", file=f)
        return sampleA_indices, sampleB_indices
