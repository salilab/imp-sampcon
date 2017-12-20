import IMP
import IMP.atom
import IMP.rmf
import subprocess
from subprocess import Popen
import os,sys,string,math
import shutil
import random
import glob

class GoodScoringModelSelector(object):
    # Authors: Shruthi Viswanath
    
    ''' Select good-scoring models based on scores and/or data satisfaction.
    Exrtact the corresponding RMFs and put them in a separate directory
    '''
    
    def __init__(self,run_directory,run_prefix):
        """Constructor.
        @param run_directory the directory containing subdirectories of runs
        @param run_prefix the prefix for each run directory. For e.g. if the subdirectories are modeling_run1, modeling_run2, etc. the prefix is modeling_run
        """
        self.run_dir=run_directory
        self.run_prefix=run_prefix
        self.all_good_scoring_models=[]# list with each member as a tuple (run id,replica id,frame id) corresponding to good-scoring models

    def _get_subfields_for_criteria(self,field_headers,selection_keywords_list,printing_keywords_list):
        ''' Given the list of keywords, get all the stat file entries corresponding to each keyword.'''

        selection_fields=[[] for j in range(len(selection_keywords_list))] # list of lists corresponding to field indices for each keyword

        printing_fields = [-1 for j in range(len(printing_keywords_list))] #just a placeholder

        for fh_index in field_headers:
            
            for ki,kw in enumerate(selection_keywords_list):
                if kw in field_headers[fh_index]:
                    selection_fields[ki].append(fh_index)
                    
            for ki,kw in enumerate(printing_keywords_list):
                if kw in field_headers[fh_index]:
                    printing_fields[ki] = fh_index
            
        return selection_fields,printing_fields
    
    def _get_crosslink_satisfaction(self,crosslink_distance_values,crosslink_percentage_lower_threshold,
                                     crosslink_percentage_upper_threshold,xlink_distance_lower_threshold,xlink_distance_upper_threshold):
        ''' For crosslinks, we want models with atleast x% (e.g. 90%) or more crosslink satisfaction. A crosslink is satisfied if the distance is between the lower and upper distance thresholds 
        @param crosslink_distance_values values of distances in the current model
        @param crosslink_percentage_lower_threshold atleast x% of crosslinks should be within the below distance thresholds
        @param crosslink_percentage_upper_threshold atmost x% of crosslinks should be within the below distance thresholds (usually 100%: dummy parameter)
        @param xlink_distance_lower_threshold a crosslink should be atleast this distance apart (usually 0) to be considered satisfied
        @param xlink_distance_upper_threshold a crosslink should be atmost this distance apart to be considered satisfied 
        '''
        satisfied_xlinks=0.0
        for d in crosslink_distance_values:
            if d>=xlink_distance_lower_threshold and d<=xlink_distance_upper_threshold:
                satisfied_xlinks+=1.0

        percent_satisfied=satisfied_xlinks/float(len(crosslink_distance_values))

        if percent_satisfied>=crosslink_percentage_lower_threshold and percent_satisfied<=crosslink_percentage_upper_threshold:
            return percent_satisfied,True
        else:
            return percent_satisfied,False

    
    def _get_score_satisfaction(self,score,lower_threshold,upper_threshold):
        ''' Check if the score is within the thresholds
        '''
        if score<=upper_threshold and score>=lower_threshold:
            return True
        return False

    def _extract_models_from_trajectories(self,output_dir): 
        ''' Given the list of all good-scoring model indices, extract their frames and store them ordered by the list index.
        '''
        
        for i,gsm in enumerate(self.all_good_scoring_models):
            (runid,replicaid,frameid)=gsm 
            
            trajfile=os.path.join(self.run_dir,self.run_prefix+runid,'output','rmfs',replicaid+'.rmf3')

            slice_location=os.path.join(os.environ['IMP_BIN_DIR'],'rmf_slice')
            
            #rmf_slice=Popen([slice_location,trajfile,"-f",str(frameid),os.path.join(output_dir,str(i)+'.rmf3')])
            #out,err=rmf_slice.communicate()
            
            rmf_slice = subprocess.call([slice_location,trajfile,"-f",str(frameid),os.path.join(output_dir,str(i)+'.rmf3')])
            
            

    def get_good_scoring_models(self,selection_keywords_list=[],printing_keywords_list=[],aggregate_lower_thresholds=[],
                                        aggregate_upper_thresholds=[],member_lower_thresholds=[],member_upper_thresholds=[],extract=False): 
        ''' Loops over all stat files in the run directory and populates the list of good-scoring models.
        @param selection_keywords_list is the list of keywords in the PMI stat file that need to be checked for each datatype/score in the criteria list
        @param printing_keywords_list is the list of keywords in the PMI stat file whose values needs to be printed for selected models
        @param aggregate_lower_thresholds The list of lower bounds on the values corresponding to fields in the criteria_list. Aggregates are used for terms like % of crosslink satisfaction and thresholds of score terms 
        @param aggregate_upper_thresholds The list of upper bounds on the values corresponding to fields in the criteria_list. Aggregates are used for terms like % of crosslink satisfaction and thresholds of score terms
        @param member_lower_thresholds The list of lower bounds for values of subcomponents of an aggregate term. E.g. for crosslink satisfaction the thresholds are on distances for each individual crosslink. For score terms this can be ignored since thresholds are mentioned in the aggregate fields.
        @param member_upper_thresholds The list of upper bounds for values of subcomponents of an aggregate term. E.g. for crosslink satisfaction the thresholds are on distances for each individual crosslink. For score terms this can be ignored since thresholds are mentioned in the aggregate fields.
        '''
        if extract:
            output_dir=os.path.join(self.run_dir,"good_scoring_models")
        else:
            output_dir=os.path.join(self.run_dir,"filter")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir, ignore_errors=True)
        os.mkdir(output_dir)
        
        outf=open(os.path.join(output_dir,"model_ids_scores.txt"),'w')
        # Header line first
        print >>outf,"Model_index","Run_id", "Replica_id","Frame_id",
        for skw in selection_keywords_list:
            print >>outf,skw,
        for pkw in printing_keywords_list:
            print >>outf,pkw,
        
        print >>outf
     
        num_runs = 0 

        for each_run_dir in sorted(os.listdir(self.run_dir)):  
         
            if not each_run_dir.startswith(self.run_prefix):
                continue

            runid=each_run_dir.split(self.run_prefix)[1]
        
            num_runs+=1
            
            print "Analyzing",runid
           
            for each_replica_stat_file in sorted(glob.glob(os.path.join(self.run_dir,each_run_dir,"output")+"/stat.*.out")):
                            
                    replicaid=each_replica_stat_file.strip(".out").split(".")[-1]

                    rsf=open(each_replica_stat_file,'r')
                    
                    for line_index,each_model_line in enumerate(rsf.readlines()): # for each model in the current replica
                        
                        if line_index==0:
                            field_headers=eval(each_model_line.strip())
                            fields_for_selection,fields_for_printing=self._get_subfields_for_criteria(field_headers,selection_keywords_list,printing_keywords_list)
                            continue
                        
                        frameid=line_index-1

                        dat=eval(each_model_line.strip())
                        
                        model_satisfies=False
                        selection_criteria_values=[]
                                                
                        for si,score_type in enumerate(selection_keywords_list):
                            if "crosslink" in score_type.lower() and "distance" in score_type.lower():
                                crosslink_distance_values=[float(dat[j]) for j in fields_for_selection[si]] # TODO : consider ambiguity
                      
                                satisfied_percent,model_satisfies=self._get_crosslink_satisfaction(crosslink_distance_values,aggregate_lower_thresholds[si],aggregate_upper_thresholds[si],member_lower_thresholds[si],member_upper_thresholds[si])
                                
                                selection_criteria_values.append(satisfied_percent)

                            else:
                                score_value=float(dat[fields_for_selection[si][0]])

                                model_satisfies=self._get_score_satisfaction(score_value,aggregate_lower_thresholds[si],aggregate_upper_thresholds[si])
                                selection_criteria_values.append(score_value)
                      
                            if not model_satisfies:
                                break

                        if model_satisfies:
                            
                            # Now get the printing criteria
                            printing_criteria_values=[]
                            for si,score_type in enumerate(printing_keywords_list):
                                score_value=float(dat[fields_for_printing[si]])
                                printing_criteria_values.append(score_value)
                                
                            self.all_good_scoring_models.append((runid,replicaid,frameid))
                            
                            # Print out the scores finally
                            
                            print >>outf,len(self.all_good_scoring_models)-1,runid,replicaid,frameid,
                            
                            for scv in selection_criteria_values:
                                print >>outf,"%.2f" %(scv),
                            for pcv in printing_criteria_values:    
                                print >>outf,"%.2f" %(pcv),
                            print >>outf
        
                    rsf.close()

        if extract:
            self._extract_models_from_trajectories(output_dir) 
        
            self._split_good_scoring_models_into_two_subsets(output_dir,num_runs,split_type="divide_by_run_ids")

       
    def _split_good_scoring_models_into_two_subsets(self,output_dir,num_runs,split_type="divide_by_run_ids"):
        ''' Get the listof good scoring models and split them into two samples, keeping the models in separate directories. 
        @param split_type how to split good scoring models into two samples. Current options are:
        (a) divide_by_run_ids : where the list of runids is divided into 2. e.g. if we have runs from 1-50, good scoring models from runs
        1-25 is sample A and those from runs 26-50 is sample B. 
        (b) random : split the set of good scoring models into two subsets at random.
        '''
        sampleA_indices=[]
        sampleB_indices=[]
        
        if split_type=="divide_by_run_ids": # split based on run ids
            
            half_num_runs= num_runs/2
            for i,gsm in enumerate(self.all_good_scoring_models):
                if int(gsm[0])<=half_num_runs:   
                    sampleA_indices.append(i)
                else: 
                    sampleB_indices.append(i)

        elif split_type=="random":
            # not implemented!
            sampleA_indices=[i for i in range(random.sample(len(self.all_good_scoring_models),len(self.all_good_scoring_models)/2))]
            sampleB_indices=[i for i in range(self.all_good_scoring_models) if i not in sampleA_indices]

        # write model and sample IDs to a file
        f=open(os.path.join(self.run_dir,'good_scoring_models','model_sample_ids.txt'),'w')
        
        # move models to corresponding sample directory
        sampleA_dir = os.path.join(output_dir,"sample_A")
        sampleB_dir = os.path.join(output_dir,"sample_B")
        os.mkdir(sampleA_dir)      
        os.mkdir(sampleB_dir)
        
        for i in sampleA_indices:
            print >>f,i,"A"
            shutil.move(os.path.join(output_dir,str(i)+'.rmf3'),os.path.join(sampleA_dir,str(i)+'.rmf3'))
        for i in sampleB_indices:
            print >>f,i,"B"
            shutil.move(os.path.join(output_dir,str(i)+'.rmf3'),os.path.join(sampleB_dir,str(i)+'.rmf3'))
        f.close()
        
