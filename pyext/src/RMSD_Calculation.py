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
        mh = IMP.atom.read_pdb(str_file, m,
                               IMP.atom.NonWaterNonHydrogenPDBSelector())
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
        mh = IMP.atom.read_pdb(str_file, m,
                               IMP.atom.NonWaterNonHydrogenPDBSelector())
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
                pts.append([p.get_coordinates()[i] for i in range(3)])
                
                if num == 0 and sample_name=='A':
                    masses.append(IMP.atom.Mass(leaf).get_mass())
                    radii.append(p.get_radius())
                    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(leaf))
                    # traverse up the Hierarchy to get copy number of the molecule that the bead belongs to
                    copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(leaf)) 
                    
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


def get_rmfs_coordinates_one_rmf(path, rmf_A, rmf_B, subunit_name):

    # Open RMFs and get total number of models
    rmf_fh = RMF.open_rmf_file_read_only(rmf_A)
    n_models = [rmf_fh.get_number_of_frames()]
    rmf_fh = RMF.open_rmf_file_read_only(rmf_B)
    n_models.append(rmf_fh.get_number_of_frames())

    masses = []
    radii = []
    ps_names = []

    models_name = []

    # Build hierarchy from the RMF file
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
    IMP.rmf.load_frame(rmf_fh, 0)
    m.update()

    ######
    # Initialize output array

    pts = 0 # number of particles in each model

    #   Get selection
    if subunit_name:
        s0 = IMP.atom.Selection(h, resolution=1, molecule=subunit_name)
    else:
        s0 = IMP.atom.Selection(h, resolution=1)
    #   Count particles
    for leaf in s0.get_selected_particles():
        p=IMP.core.XYZR(leaf)
        pts+=1
    #   Initialize array
    conform = np.empty([n_models[0]+n_models[1], pts, 3])

    mod_id = 0 # index for each model in conform.
    for rmf_file in [rmf_A, rmf_B]:

        models_name.append(rmf_file)

        rmf_fh = RMF.open_rmf_file_read_only(rmf_file)
        h = IMP.rmf.create_hierarchies(rmf_fh, m)[0]
        
        print("Opening RMF file:", rmf_file, "with", rmf_fh.get_number_of_frames(), "frames")
        for f in range(rmf_fh.get_number_of_frames()):
            if f%100==0:
                #pass
                print("  -- Opening frame", f, "of", rmf_fh.get_number_of_frames())
            IMP.rmf.load_frame(rmf_fh, f)

            m.update()
            if subunit_name:
                s0 = IMP.atom.Selection(h, resolution=1, molecule=subunit_name)
            else:
                s0 = IMP.atom.Selection(h, resolution=1)
            particles = s0.get_selected_particles()
            # Copy particle coordinates
            for i in range(len(particles)):
                leaf=particles[i]
                p=IMP.core.XYZR(leaf)
                pxyz = p.get_coordinates()
                conform[mod_id][i][0] = pxyz[0]
                conform[mod_id][i][1] = pxyz[1]
                conform[mod_id][i][2] = pxyz[2]

                # Just for the first model, update the masses and radii and log the particle name in ps_names
                if mod_id == 0 and rmf_file==rmf_A:
                    masses.append(IMP.atom.Mass(leaf).get_mass())
                    radii.append(p.get_radius())
                    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(leaf))
                    # traverse up the Hierarchy to get copy number of the molecule that the bead belongs to
                    copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(leaf))
                    if IMP.atom.Fragment.get_is_setup(leaf): #TODO not tested on non-fragment systems
                        residues_in_bead = IMP.atom.Fragment(leaf).get_residue_indexes()
                        ps_names.append(mol_name+"_"+str(min(residues_in_bead))+"_"+str(max(residues_in_bead))+"_"+copy_number)
                    else:
                        residue_in_bead = str(IMP.atom.Residue(leaf).get_index())
                        ps_names.append(mol_name+"_"+residue_in_bead+"_"+residue_in_bead+"_"+copy_number)
            mod_id+=1

    return ps_names, masses, radii, conform, models_name, n_models

def get_rmsds_matrix(conforms, mode, sup, cores):
    print("Mode:",mode,"Superposition:",sup,"Number of cores:",cores)

    if(mode=="cpu_serial" and not sup):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("NOSUP_OMP_CALCULATOR", conforms)

    elif(mode=="cpu_omp" and not sup):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("NOSUP_OMP_CALCULATOR", conforms)
        calculator.setNumberOfOpenMPThreads(cores)

    elif(mode=="cpu_omp" and sup):
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", conforms)
        calculator.setNumberOfOpenMPThreads(cores)

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
