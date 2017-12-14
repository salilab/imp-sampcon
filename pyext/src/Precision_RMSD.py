import os,sys
import numpy
import string
import IMP
import IMP.em
import pyRMSD.RMSDCalculator

def parse_custom_ranges(ranges_file):
    fl = open(ranges_file, 'r')
    density_custom_ranges = fl.readlines()[0].strip()
    exec(density_custom_ranges)
       
    fl.close()
    
    return density_custom_ranges
    
    
def get_particles_from_superposed(cluster_conform_i, cluster_conform_0, masses, radii, align): 

    m=IMP.Model()
    ps = []
    if align:
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_SERIAL_CALCULATOR", numpy.array([cluster_conform_0, cluster_conform_i]))
    else:
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("NOSUP_SERIAL_CALCULATOR", numpy.array([cluster_conform_0, cluster_conform_i]))

    rmsd, superposed_fit = calculator.pairwise(0, 1, get_superposed_coordinates = True)

    for particle_index in range(len(superposed_fit[1])): 
        p = IMP.Particle(m, "%s" % str(particle_index))
    
        IMP.core.XYZ.setup_particle(p, superposed_fit[1][particle_index])
        IMP.core.XYZR.setup_particle(p, float(radii[particle_index]))
        IMP.atom.Mass.setup_particle(p, float(masses[particle_index]))
        m.update()
        ps.append(p)

    return rmsd, m, ps

class GetModelDensity(object):
    """Compute mean density maps from structures.
    Keeps a dictionary of density maps,
    keys are in the custom ranges. When you call add_subunits_density, it adds
    particle coordinates to the existing density maps.
    """

    def __init__(self, custom_ranges=None, resolution=20.0, voxel=5.0, bead_names = None):
        """Constructor.
        @param list of particles decorated with mass, radius, and XYZ
           @param resolution The MRC resolution of the output map (in Angstrom unit)
           @param voxel The voxel size for the output map (lower is slower)
        """

        self.MRCresolution = resolution
        self.voxel = voxel
        self.count_models = 0.0
        self.densities={}
        self.bead_names = bead_names
        self.custom_ranges=custom_ranges
        
        # for each custom range get the particle indices that will be added to the density for that custom range
        self.particle_indices_in_custom_ranges={}
        
        for density_name in self.custom_ranges:
            self.particle_indices_in_custom_ranges[density_name]=[]
        
        # go through each bead, put it in the appropriate custom range(s)
        for index,beadname in enumerate(self.bead_names):
            for density_name in self.custom_ranges:
                for domain in self.custom_ranges[density_name]: # each domain in the list custom_ranges[density_name]
                     if self._is_contained(beadname,domain):
                        self.particle_indices_in_custom_ranges[density_name].append(index)
                        #print beadname,"is in",domain
                        break # already added particle to this custom range
    
    def normalize_density(self):
        pass
   
    def _create_density_from_particles(self, ps, name,
                                      kernel_type='GAUSSIAN'):
        '''Internal function for adding to densities.
        pass XYZR particles with mass and create a density from them.
        kernel type options are GAUSSIAN, BINARIZED_SPHERE, and SPHERE.'''
        kd = {
            'GAUSSIAN': IMP.em.GAUSSIAN,
            'BINARIZED_SPHERE': IMP.em.BINARIZED_SPHERE,
            'SPHERE': IMP.em.SPHERE}
        dmap = IMP.em.SampledDensityMap(ps, self.MRCresolution, self.voxel)
        dmap.calcRMS()
        dmap.set_was_used(True)

        if name not in self.densities:
            self.densities[name] = dmap
        else:
            bbox1 = IMP.em.get_bounding_box(self.densities[name])
            bbox2 = IMP.em.get_bounding_box(dmap)
            bbox1 += bbox2
            dmap3 = IMP.em.create_density_map(bbox1,self.voxel)
            dmap3.set_was_used(True)
            dmap3.add(dmap)
            dmap3.add(self.densities[name])
            self.densities[name] = dmap3
    
    def _is_contained(self,bead_name,domain):
        """ domain can be the name of a single protein or a tuple (protein_name,start_residue,end_residue)
        bead is a string of type moleculeName_startResidue_endResidue
        """
        if type(domain)==tuple:
            bead_res_start,bead_res_end,bead_protein = bead_name.split("_")
            bead_residues = set(range(int(bead_res_start),int(bead_res_end)+1))
            domain_protein = domain[0]
            domain_residues = set(range(int(domain[1]),int(domain[2])+1))
            
            if bead_protein == domain_protein and not domain_residues.isdisjoint(bead_residues):
                return True
            
        else:
            if domain in bead_name:
                return True
        
        return False
                
    def add_subunits_density(self, ps):
        """Add a frame to the densities.
        @param ps List of particles decorated with XYZR and Mass.
        """
        self.count_models += 1.0
        # initialize custom list of particles
        particles_custom_ranges={}
        for density_name in self.custom_ranges:
            particles_custom_ranges[density_name]=[]
        
        # add each particle to the relevant custom list
        for density_name in self.custom_ranges:
            for particle_index in self.particle_indices_in_custom_ranges[density_name]:
                particles_custom_ranges[density_name].append(ps[particle_index])
                                    
        #finally, add each custom particle list to the density
        for density_name in self.custom_ranges:
            self._create_density_from_particles(particles_custom_ranges[density_name],density_name)
     
    def get_density_keys(self):
        return list(self.densities.keys())

    def get_density(self,name):
        """Get the current density for some component name"""
        if name not in self.densities:
            return None
        else:
            return self.densities[name]
        
    def write_mrc(self, path="./",file_prefix=""):
        for density_name in self.densities:
            self.densities[density_name].multiply(1. / self.count_models)
            IMP.em.write_map(
                self.densities[density_name],
                path + "/" + file_prefix + "_"+ density_name + ".mrc",
                IMP.em.MRCReaderWriter())
