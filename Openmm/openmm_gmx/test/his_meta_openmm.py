from gmx_openmm import *

class HisMetadynamics(GMXSimulation):
    def __init__(self, gro_file, top_file, subset=None,
             temperature=300*u.kelvin, pressure=1.0*u.atmosphere,
             timestep=2.0*u.femtosecond,
             time_per_frame=50*u.picosecond,
             residue=168, biasFactor=1.5, height=1.0*u.kilojoule_per_mole, frequency=100,
             output_basename=None, output_index=None):
        
    
        GMXSimulation.__init__(self,gro_file=gro_file,
        top_file=top_file,
        subset=subset,
        temperature=temperature,
        pressure=pressure,
        timestep=timestep,
        time_per_frame=time_per_frame,
        output_index=output_index,
        output_basename=output_basename)
        
        
        self.residue = residue       
        self.biasFactor = biasFactor
        self.height = height
        self.frequency = frequency
        self.variables = self.get_chi1_chi2()
        

    def run_meta(self, time):
    
        # simulation time
        n_steps_per_frame = 2500 # 2fs/step
        n_steps = int(time / (2*u.femtoseconds)) + 1

        biasDir = self.output_folder + f'/His{self.residue}_biasDir'
        create_dir(biasDir)
        
        self.system_meta = build_system(self.top)
        variables = self.get_chi1_chi2()
        self.meta = metadynamics.Metadynamics(self.system_meta, self.variables, self.temperature, self.biasFactor, 
                                            self.height, self.frequency, saveFrequency=100, biasDir=biasDir)

        integrator = build_integrator(self.temperature, self.timestep)
        simu = Simulation(self.top.topology, self.system_meta, integrator)
        
        prog_reporter = parmed.openmm.ProgressReporter(
           f=f"{self.output_name}_meta_His{self.residue}.progress",
           reportInterval=n_steps_per_frame,
           totalSteps=n_steps,
           potentialEnergy=True,
           kineticEnergy=True,
           totalEnergy=False,
           temperature=True,
           volume=True)
               
        full_reporter = md.reporters.DCDReporter(
            f"{self.output_name}_meta_His{self.residue}.dcd",
            reportInterval=n_steps_per_frame)

        simu.reporters.append(prog_reporter)
        simu.reporters.append(full_reporter)

        if hasattr(self, 'simulation_eq'):
            simu = copy_simu_state(self.simulation_eq, simu)
            self.meta.step(simu, n_steps)
            simu.saveCheckpoint(f'{self.output_name}_meta.chk')
        else:
            raise ValueError('Run simulation_eq first!')

        self.free_energy = self.meta.getFreeEnergy()
        np.save(f'{self.output_name}_free_energy_His{self.residue}.npy', self.free_energy)
        print(self.free_energy)


    def get_chi1(self, residue, md_top):
        N = md_top.select(f'chainid 0 and residue {residue} and name N ')
        CA = md_top.select(f'chainid 0 and residue {residue} and name CA')
        CB = md_top.select(f'chainid 0 and residue {residue} and name CB')
        CG = md_top.select(f'chainid 0 and residue {residue} and name CG')
        return np.array([N, CA, CB, CG]).transpose()[0]


    def get_chi2(self, residue, md_top):
        CA = md_top.select(f'chainid 0 and residue {residue} and name CA ')
        CB = md_top.select(f'chainid 0 and residue {residue} and name CB')
        CG = md_top.select(f'chainid 0 and residue {residue} and name CG')
        ND1 = md_top.select(f'chainid 0 and residue {residue} and name ND1')
        return np.array([CA, CB, CG, ND1]).transpose()[0]


    def get_chi1_chi2(self):
        'return chi1 and chi2 collective varibles of a His'
        
        pdbfile = self.gro_file[:-4]+'.pdb'
        traj = md.load(pdbfile)
        md_top = traj.topology

        cv1 = CustomTorsionForce('theta')
        cv1.addTorsion(*self.get_chi1(self.residue, md_top)[:])
        chi1 = BiasVariable(cv1, -np.pi, np.pi, 0.5, True)
        cv2 = CustomTorsionForce('theta')
        cv2.addTorsion(*self.get_chi2(self.residue, md_top)[:])
        chi2 = BiasVariable(cv2, -np.pi, np.pi, 0.5, True)

        return [chi1, chi2]


def meta_make_parser():
    from argparse import ArgumentParser
    parser = ArgumentParser()

    parser.add_argument('--residue', default=168, type=int,
                        help='His as collective varibles')
    
    
    parser.add_argument('--biasFactor', default=1.5, type=float,
                        help='biasFactor for metadynamics')

    parser.add_argument('--height', default=1.0, type=float,
                        help='height for metadynamics')
        
    parser.add_argument('--frequency', default=100, type=float,
                        help='frequency for metadynamics')
    return parser
