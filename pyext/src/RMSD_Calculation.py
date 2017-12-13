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
        
    f1=open(idfile_A, 'w+')
    f2=open(idfile_B, 'w+')

    for str_file in sorted(glob.glob("%s/sample_A/*.pdb" % path)):
        print >>f1, str_file, num
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

        
    for str_file in sorted(glob.glob("%s/sample_B/*.pdb" % path)):
        print >>f2, str_file, num
        m = IMP.Model()
        mh = IMP.atom.read_pdb(file, m,IMP.atom.NonWaterNonHydrogenPDBSelector())
        mps = IMP.core.get_leaves(mh)
        pts = [IMP.core.XYZ(p).get_coordinates() for p in mps]
        conform.append(pts)
        pts = []   
        num = num + 1
        
    return np.array(conform), masses 

def get_rmfs_coordinates(path, idfile_A, idfile_B):

    conform = []
    num = 0
    masses = []
    radii = []
    ps_names = []
    
    f1=open(idfile_A, 'w+')
    f2=open(idfile_B, 'w+')

    models_name = []
    
    for str_file in sorted(glob.glob("%s/sample_A/*.rmf3" % path)):
        print >>f1, str_file, num
        models_name.append(str_file)

        m = IMP.Model()
        inf = RMF.open_rmf_file_read_only(str_file)
        h = IMP.rmf.create_hierarchies(inf, m)[0]
        IMP.rmf.load_frame(inf, 0)
        
        pts = []

        s0 = IMP.atom.Selection(h, resolution=1)
        for leaf in s0.get_selected_particles():
            
            p=IMP.core.XYZR(leaf)
            pts.append(p.get_coordinates())
            if num == 0:
                masses.append(IMP.atom.Mass(leaf).get_mass())
                radii.append(p.get_radius())
                ps_names.append(IMP.atom.get_molecule_name(IMP.atom.Hierarchy(leaf)))
        conform.append(pts)
        pts = []
        num = num + 1
        
    for str_file in sorted(glob.glob("%s/sample_B/*.rmf3" % path)):
        print >>f2, str_file, num
        models_name.append(str_file)
        m = IMP.Model()
        inf = RMF.open_rmf_file_read_only(str_file)
        h = IMP.rmf.create_hierarchies(inf, m)[0]
        IMP.rmf.load_frame(inf, 0)

        pts = []
        s0 = IMP.atom.Selection(h, resolution=1)
        for leaf in s0.get_selected_particles():
            p=IMP.core.XYZ(leaf)
            pts.append(p.get_coordinates())

        conform.append(pts)                                                                                                                                                                                                     
        pts = []
        num = num + 1

    return ps_names, masses, radii, np.array(conform), models_name

def get_rmfs_subunit_coordinates(path, subunit_name):
    pts = []
    conform = []
    num = 0
    
    for str_file in glob.glob("%s/*.rmf3" % path):
        print >>f1, str_file, num
        m = IMP.Model()
        inf = RMF.open_rmf_file_read_only(str_file)
        h = IMP.rmf.create_hierarchies(inf, m)[0]
        particle2s = IMP.core.get_leaves(h)
        IMP.rmf.load_frame(inf, 0)

        s0 = IMP.atom.Selection(h, molecule=subunit_name, resolution=1)
        for leaf in s0.get_selected_particles():
            if "bead" not in leaf.get_name():
                p=IMP.core.XYZ(leaf)
                pts.append(p.get_coordinates())

        conform.append(pts)
        pts = []
        num = num + 1
                
    for str_file in glob.glob("%s/*.rmf3" % path):
        print >>f2, str_file, num
        m = IMP.Model()
        inf = RMF.open_rmf_file_read_only(str_file)
        h = IMP.rmf.create_hierarchies(inf, m)[0]
        particle2s = IMP.core.get_leaves(h)
        IMP.rmf.load_frame(inf, 0)

        s0 = IMP.atom.Selection(h, molecule=subunit_name, resolution=1)
        for leaf in s0.get_selected_particles():
            if "bead" not in leaf.get_name():
                p=IMP.core.XYZ(leaf)
                pts.append(p.get_coordinates())

        conform.append(pts)                                                                                                                                                                                                     
        pts = []
        num = num + 1

    return np.array(conform)


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
