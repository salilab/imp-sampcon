from __future__ import print_function
import pyRMSD.RMSDCalculator
from pyRMSD.matrixHandler import MatrixHandler
from pyRMSD.condensedMatrix import CondensedMatrix

import numpy as np
import sys, os, glob

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
        print(str_file, num, file=f1)
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
        print(str_file, num, file=f2)
        models_name.append(str_file)
        
        m = IMP.Model()
        mh = IMP.atom.read_pdb(file, m,IMP.atom.NonWaterNonHydrogenPDBSelector())
        mps = IMP.core.get_leaves(mh)
        pts = [IMP.core.XYZ(p).get_coordinates() for p in mps]
        conform.append(pts)
        pts = []   
        num = num + 1
        
    return np.array(conform), masses, radii, models_name

def get_rmfs_coordinates(path, idfile_A, idfile_B, subunit_name):

    conform = []
    num = 0
    masses = []
    radii = []
    ps_names = []
    
    f1=open(idfile_A, 'w+')
    f2=open(idfile_B, 'w+')

    models_name = []
    
    for sample_name,sample_id_file in zip(['A','B'],[f1,f2]):
        
        for str_file in sorted(glob.glob("%s/sample_%s/*.rmf3" % (path,sample_name)),key=lambda x:int(x.split('/')[-1].split('.')[0])):
            print(str_file, num, file=sample_id_file)
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
                    copy_number = "X"
                    # Need to find the copy number from the molecule
                    # In PMI, this is three levels above the individual residues/beads
                    mol_p = IMP.atom.Hierarchy(p).get_parent().get_parent().get_parent()
                    if IMP.atom.Copy().get_is_setup(mol_p):
                        copy_number = str(IMP.atom.Copy(mol_p).get_copy_index())
                    
                    if IMP.atom.Fragment.get_is_setup(leaf): #TODO not tested on non-fragment systems
                        residues_in_bead = IMP.atom.Fragment(leaf).get_residue_indexes()
                        
                        ps_names.append(mol_name+"_"+str(min(residues_in_bead))+"_"+str(max(residues_in_bead))+"_"+copy_number)
                            
                    else:
                        residue_in_bead = str(IMP.atom.Residue(leaf).get_index())
                        
                        ps_names.append(mol_name+"_"+residue_in_bead+"_"+residue_in_bead+"_"+copy_number)
            
            conform.append(pts)
            pts = []
            num = num + 1
        
    return ps_names, masses, radii, np.array(conform), models_name

def get_rmsds_matrix(conforms, mode, sup, cores):
    print("Mode:",mode,"Superposition:",sup,"Number of cores:",cores)

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
        print("Wrong values to pyRMSD ! Please Fix")
        exit()

    rmsd = calculator.pairwiseRMSDMatrix()
    rmsd_matrix=CondensedMatrix(rmsd)
    inner_data = rmsd_matrix.get_data()
    np.save("Distances_Matrix.data", inner_data)

    return inner_data
