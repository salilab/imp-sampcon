import os,sys
import numpy
import string
import IMP
import IMP.em
import pyRMSD.RMSDCalculator

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

    def __init__(self, custom_ranges=None, resolution=20.0, voxel=5.0, molnames = None):
        """Constructor.
        @param list of particles decorated with mass, radius, and XYZ
           @param resolution The MRC resolution of the output map (in Angstrom unit)
           @param voxel The voxel size for the output map (lower is slower)
        """

        self.MRCresolution = resolution
        self.voxel = voxel
        self.count_models = 0.0
        self.densities={}
        self.molnames = molnames
        self.custom_ranges=custom_ranges
    
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
        for particle, molname in zip(ps, self.molnames):
            for density_name in self.custom_ranges:
                if molname in self.custom_ranges[density_name]: 
                    particles_custom_ranges[density_name].append(particle)
                    
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
