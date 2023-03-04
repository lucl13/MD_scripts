from __future__ import print_function

from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np 
import openmmtools.integrators
import mdtraj as md
import parmed
import re
from simtk import unit as u
from sys import stdout
import copy

import logging
logger = logging.getLogger(__name__)

class GMXSimulation(object):
    """
    Create a basic OpenMM simulation from Gromacs files.

    The resulting simulation has sane defaults for biomolecular simulation.
    There's only limited customizability -- if you want more
    customizability, just use OpenMM directly!


    """
    def __init__(self, gro_file, top_file, subset=None,
                 temperature=300*u.kelvin, pressure=1.0*u.atmosphere,
                 timestep=2.0*u.femtosecond,
                 time_per_frame=50*u.picosecond,
                 output_basename=None, output_index=None):
        self.timestep = timestep
        self.temperature = temperature
        self.pressure = pressure
        # load the stuff from the files
        self.gro = GromacsGroFile(gro_file)
        self.gro_file = gro_file
        self.subset = subset
        box_vectors = Vec3(*self.gro.getPeriodicBoxVectors())
        self.top = GromacsTopFile(top_file, periodicBoxVectors=box_vectors)


        # set up the reporters
        if output_basename is None:
            output_basename = re.sub('\.gro$', '', gro_file)
        if output_index is not None:
            output_basename += str(output_index)  # make 0-padded?
        self.output_basename = output_basename
        self.n_steps_per_frame = int(time_per_frame / timestep)

        # create a output_folder

        self.output_folder = 'output'
        create_dir(self.output_folder)

        # output_name

        self.output_name = self.output_folder + '/' + self.output_basename
        
        
    
    def read_chk(self,simulation_type):
         
        # build three simulations
        self.build_simulation_cst()
        self.build_simulation_eq()
        self.build_simulation_pro()
        
        if simulation_type == 'cst':
            self.simulation_cst.loadCheckpoint(f'{self.output_name}_cst.chk')
            print(f'{simulation_type} chk loaded!')
        
        if simulation_type == 'eq':

            self.simulation_eq.loadCheckpoint(f'{self.output_name}_eq.chk')
            print(f'{simulation_type} chk loaded!')

        
        if simulation_type == 'pro':

            self.simulation_pro.loadCheckpoint(f'{self.output_name}_pro.chk')
            print(f'{simulation_type} chk loaded!')

        
        
    def run_cst(self, gro_file, time, output_file=None):
        
        
        self.build_simulation_cst()
        
        n_steps = int(time / self.timestep) + 1
        prog_reporter = parmed.openmm.ProgressReporter(
            f=self.output_name + "_cst.progress",
            reportInterval=self.n_steps_per_frame,
            totalSteps=n_steps,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=False,
            temperature=True,
            volume=True
        )
        
   
        # full reporter, including velocities
        full_reporter = md.reporters.DCDReporter(
            self.output_name + "_cst.dcd",
            reportInterval=self.n_steps_per_frame
        )
    
        

        self.simulation_cst.reporters.append(prog_reporter)
        self.simulation_cst.reporters.append(full_reporter)
        
        ##; Harmonic restraind
        traj = md.load(gro_file[:-4]+'.pdb')
        mdtraj_topol = traj.topology
        force = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
        k = 1000 * u.kilojoule_per_mole / u.nanometer**2
        force.addGlobalParameter("k", k)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        for i in mdtraj_topol.select('protein'):
            r0 = traj.xyz[0][i]*u.nanometers
            force.addParticle(i, [r0[0], r0[1], r0[2]])

        self.simulation_cst.system.addForce(force)
        
        self.initialize_random_setup(self.gro_file, self.temperature * u.kelvin)

        self.simulation_cst.minimizeEnergy()
        self.simulation_cst.step(n_steps)
        self.simulation_cst.saveCheckpoint(f'{self.output_name}_cst.chk')
        
    def run_eq(self, time):
                
        self.build_simulation_eq()
        
        n_steps = int(time / self.timestep) + 1
        prog_reporter = parmed.openmm.ProgressReporter(
            f=self.output_name + "_eq.progress",
            reportInterval=self.n_steps_per_frame,
            totalSteps=n_steps,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=False,
            temperature=True,
            volume=True
        )
        
        
        
        # full reporter, including velocities
        full_reporter = md.reporters.DCDReporter(
            self.output_name + "_eq.dcd",
            reportInterval=self.n_steps_per_frame
        )
        
        self.simulation_eq.reporters.append(prog_reporter)
        self.simulation_eq.reporters.append(full_reporter)     
        
        if hasattr(self, 'simulation_cst'):
    
            self.simulation_eq = copy_simu_state(self.simulation_cst, self.simulation_eq)
            self.simulation_eq.step(n_steps)
            self.simulation_eq.saveCheckpoint(f'{self.output_name}_eq.chk')
        else:
            raise ValueError('Run simulation_cst first!')

        
    def run_pro(self, time):
        
        
        self.build_simulation_pro()
        
        n_steps = int(time / self.timestep) + 1
        
        # log reporter (stdout)
        log_reporter = StateDataReporter(
            stdout,
            reportInterval=self.n_steps_per_frame,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=False,
            temperature=True,
            volume=True
        )

        prog_reporter = parmed.openmm.ProgressReporter(
            f=self.output_name + "_pro.progress",
            reportInterval=self.n_steps_per_frame,
            totalSteps=n_steps,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=False,
            temperature=True,
            volume=True
        )
        
        
        # full reporter, including velocities
        full_reporter = md.reporters.DCDReporter(
            self.output_name + "_pro.dcd",
            reportInterval=self.n_steps_per_frame
        )
        
        self.simulation_pro.reporters.append(log_reporter)                         
        self.simulation_pro.reporters.append(prog_reporter)
        self.simulation_pro.reporters.append(full_reporter)
        
        
        # subset reporter (no H2O)
        if self.subset is not None:
            subset_reporter = md.reporters.DCDReporter(
                self.output_name + "_dry.dcd",
                reportInterval=self.n_steps_per_frame,
                atomSubset=self.subset
            )

        if self.subset is not None:
            self.simulation_pro.reporters.append(subset_reporter)
     

        if hasattr(self, 'simulation_eq'):
            self.simulation_pro = copy_simu_state(self.simulation_eq, self.simulation_pro)
            self.simulation_pro.step(n_steps)
            self.simulation_eq.saveCheckpoint(f'{self.output_name}_eq.chk')
        
        else:
            raise ValueError('Run simulation_eq first!')
        
        
    def _debug_state(self):
        state = self.simulation.context.getState(getPositions=True,
                                                 getVelocities=True,
                                                 getEnergy=True)
        coords = state.getPositions(asNumpy=True)
        box = state.getPeriodicBoxVectors(asNumpy=True)
        velocities = state.getVelocities(asNumpy=True)
        potential = state.getPotentialEnergy()
        kinetic = state.getKineticEnergy()
        print(coords)
        print(box)
        print(velocities)
        print(potential, kinetic)

    def initialize_random_setup(self, gro_file, temperature):
        gro = app.GromacsGroFile(gro_file)
        self.simulation_cst.context.setPositions(gro.positions)
        self.simulation_cst.context.setVelocitiesToTemperature(temperature)
        
    def build_simulation_cst(self):
        
        # set up the integrator
        integrator = build_integrator(self.temperature, self.timestep)
        system = build_system(self.top)
        system.addForce(MonteCarloBarostat(self.pressure, self.temperature, 25))
        self.simulation_cst = build_simulation(self.top, system, integrator)
    
    def build_simulation_eq(self):
        
        # set up the integrator
        integrator = build_integrator(self.temperature, self.timestep)
        system = build_system(self.top)
        system.addForce(MonteCarloBarostat(self.pressure, self.temperature, 25))
        self.simulation_eq = build_simulation(self.top, system, integrator)
        
    def build_simulation_pro(self):
        
        # set up the integrator
        integrator = build_integrator(self.temperature, self.timestep)
        system = build_system(self.top)
        system.addForce(MonteCarloBarostat(self.pressure, self.temperature, 25))
        self.simulation_pro = build_simulation(self.top, system, integrator)
        
    
def build_integrator(temperature, timestep):
    integrator = openmmtools.integrators.VVVRIntegrator(
        temperature,
        1.0 / u.picoseconds,
        timestep
    )
    integrator.setConstraintTolerance(0.00001)
    return integrator

def build_system(topology):
    system = topology.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.2*u.nanometers,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005
    )

    return system

def build_simulation(topology, system, integrator):
    simulation = app.Simulation(topology.topology, system, integrator)
    return simulation

def copy_simu_state(simu_from, simu_to):
    
    state = simu_from.context.getState(getPositions=True,getVelocities=True)
    simu_to.context.setState(state)
    return simu_to   

def make_parser():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-t', '--topology')
    parser.add_argument('--subset', default="not water")
    parser.add_argument('-T', '--temperature', default=300.0, type=float,
                        help='Simulation temperature (K)')
    parser.add_argument('-P', '--pressure', default=1.0, type=float,
                        help='Simulation pressure (atm)')
    parser.add_argument('--time-per-frame', default=10, type=float,
                        help='Time between saved frames (ps)')
    parser.add_argument('--timestep', default=2.0, type=float,
                        help='Timestep for integration (fs)')
    parser.add_argument('--device-index', default=0)
    parser.add_argument('-n', '--nruns', default=1, type=int)
    parser.add_argument("--start-run", default=0, type=int)
    #parser.add_argument('--posre', default=0.0,
                        #help="Position-restrained run time (ns)")
    #parser.add_argument('--equil', default=0.0,
                        #help="Equilibration time (ns)")
    parser.add_argument('--run', default=0.0,
                        help="Production run time (ns)")
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--log', default=None)
    parser.add_argument('gro_file')
    return parser

def prepare_logging(logfile, debug):
    logger.setLevel(logging.INFO)
    if debug:
        logger.setLever(logging.DEBUG)

    if logfile is None:
        handler = logging.StreamHandler()
    else:
        handler = logging.FileHandler(logfile)

    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

def write_pdb(gmx_simulation):
    positions = gmx_simulation.simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(gmx_simulation.simulation.topology, positions, open('output.pdb', 'w'))

    
def create_dir(dir):
    import os
    import shutil
    if os.path.exists(dir):
        #shutil.rmtree(dir)
        #print('output folder exist!')
        pass 
    else:
        os.makedirs(dir)


