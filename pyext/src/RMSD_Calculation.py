import pyRMSD.RMSDCalculator
from pyRMSD.matrixHandler import MatrixHandler
from pyRMSD.condensedMatrix import CondensedMatrix

import numpy as np
import sys, os, glob
import random

import IMP
import IMP.atom
import IMP.rmf
import RMF

def get_pdbs_coordinates(path, idfile_A, idfile_B):

    pts = []
    conform = []
    num = 0
    masses = []
    radii = []
    
    models_name = []
     
    f1=open(idfile_A, 'w+')
    f2=open(idfile_B, 'w+')

    for str_file in sorted(glob.glob("%s/sample_A/*.pdb" % path),key=lambda x:int(x.split('/')[-1].split('.')[0])):
        print >>f1, str_file, num
        models_name.append(str_file)
        
        m = IMP.Model()
        mh = IMP.atom.read_pdb(file, m,IMP.atom.NonWaterNonHydrogenPDBSelector())
        mps = IMP.core.get_leaves(mh) 
        pts = [IMP.core.XYZ(p).get_coordinates() for p in mps]
        if num == 0:
            masses = [IMP.atom.Mass(p).get_mass() for p in mps]
            radii  = [IMP.core.XYZR(p).get_radius() for p in mps]
        conform.append(pts)

        pts = []
        num = num + 1

        
    for str_file in sorted(glob.glob("%s/sample_B/*.pdb" % path),key=lambda x:int(x.split('/')[-1].split('.')[0])):
        print >>f2, str_file, num
        models_name.append(str_file)
        
        m = IMP.Model()
        mh = IMP.atom.read_pdb(file, m,IMP.atom.NonWaterNonHydrogenPDBSelector())
        mps = IMP.core.get_leaves(mh)
        pts = [IMP.core.XYZ(p).get_coordinates() for p in mps]
        conform.append(pts)
        pts = []   
        num = num + 1
        
    return np.array(conform), masses, radii, models_name

def get_rmfs_coordinates(path, idfile_A, idfile_B, subunit_name,  subsample, rmfs_lists = None):

    conform = []
    num = 0
    masses = []
    radii = []
    ps_names = []
    
    f1=open(idfile_A, 'w+')
    f2=open(idfile_B, 'w+')

    models_name = []

    rmfs = {}
    print('-----------------', rmfs_lists, subsample)
    if rmfs_lists:
        rmfs_A = []
        rmfs_B = []
        for line in open(rmfs_lists[0], 'r'):
            vals = line.split()
            #rmfs_A.append(vals[0].split('/')[-1].split('.')[0])
            rmfs_A.append(path+"/sample_A/"+vals[0])
        for line in open(rmfs_lists[1], 'r'):
            vals = line.split()
            #rmfs_B.append(vals[0].split('/')[-1].split('.')[0])
            rmfs_B.append(path+"/sample_B/"+vals[0])
        
        # Sort rmfs
        rmfs_A = sorted(rmfs_A, key=lambda x:x.split('/')[-1].split('.')[0])
        rmfs_B = sorted(rmfs_B, key=lambda x:x.split('/')[-1].split('.')[0])
        rmfs['A'] = rmfs_A
        rmfs['B'] = rmfs_B
            
    else:
        rmfs_A = sorted(glob.glob("%s/sample_A/*.rmf3" % path),key=lambda x:x.split('/')[-1].split('.')[0])
        rmfs_B = sorted(glob.glob("%s/sample_B/*.rmf3" % path),key=lambda x:x.split('/')[-1].split('.')[0])
        rmfs['A'] = rmfs_A
        rmfs['B'] = rmfs_B

    if subsample:
        
        tot = len(rmfs_A)+len(rmfs_B)
        p_A = float(len(rmfs_A))/tot
        p_B = float(len(rmfs_B))/tot
        s_A = int(subsample * p_A)
        s_B = int(subsample * p_B)
        rmfs_A = random.sample(rmfs_A,s_A)
        rmfs_B = random.sample(rmfs_B,s_B)
        rmfs['A'] = rmfs_A
        rmfs['B'] = rmfs_B
    
    for sample_name,sample_id_file in zip(['A','B'],[f1,f2]):
        
        for str_file in rmfs[sample_name]:
            print >>sample_id_file, str_file, num
            models_name.append(str_file)

            m = IMP.Model()
            inf = RMF.open_rmf_file_read_only(str_file)
            h = IMP.rmf.create_hierarchies(inf, m)[0]
            IMP.rmf.load_frame(inf, 0)
            
            pts = []

            if subunit_name:
                s0 = IMP.atom.Selection(h, resolution=1,molecule=subunit_name)
            else:
                s0 = IMP.atom.Selection(h, resolution=1)
                
         
            for leaf in s0.get_selected_particles():
                
                p=IMP.core.XYZR(leaf)
                pts.append(p.get_coordinates())
                
                if num == 0 and sample_name=='A':
                    masses.append(IMP.atom.Mass(leaf).get_mass())
                    radii.append(p.get_radius())
                    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(leaf))
                    
                    
                    if IMP.atom.Fragment.get_is_setup(leaf): #TODO not tested on non-fragment systems
                        residues_in_bead = IMP.atom.Fragment(leaf).get_residue_indexes()
                        
                        ps_names.append(mol_name+"_"+str(min(residues_in_bead))+"_"+str(max(residues_in_bead)))
                            
                    else:
                        residue_in_bead = str(IMP.atom.Residue(leaf).get_index())
                        
                        ps_names.append(mol_name+"_"+residue_in_bead+"_"+residue_in_bead)
            
   
            conform.append(pts)
            pts = []
            num = num + 1
        
    return ps_names, masses, radii, np.array(conform), models_name

def get_rmsds_matrix(conforms, mode, sup, cores):
    print "Mode:",mode,"Superposition:",sup,"Number of cores:",cores

    if(mode=="cpu_serial" and not sup):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("NOSUP_OMP_CALCULATOR", conforms)

    elif(mode=="cpu_omp" and not sup):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("NOSUP_OMP_CALCULATOR", conforms)
        calculator.setNumberOfOpenMPThreads(int(cores))

    elif(mode=="cpu_omp" and sup):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", conforms)
        calculator.setNumberOfOpenMPThreads(int(cores))

    elif(mode=="cuda" and sup):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_CUDA_MEM_CALCULATOR", conforms)  

    else:
        print "Wrong values to pyRMSD ! Please Fix"
        exit()

    rmsd = calculator.pairwiseRMSDMatrix()
    rmsd_matrix=CondensedMatrix(rmsd)
    inner_data = rmsd_matrix.get_data()
    np.save("Distances_Matrix.data", inner_data)

    return inner_data
