#!/usr/bin/python

from __future__ import print_function

import socket
import sys
import os
import subprocess
import mdtraj as md 

use_multiprocessing = True
if use_multiprocessing:
    import multiprocessing
    max_cpus = 24 # We might want to not run on the full number of cores, as Rosetta take about 2 Gb of memory per instance

###################################################################################################################################################################
# Important: The variables below are set to values that will make the run complete faster (as a tutorial example), but will not give scientifically valid results.
#            Please change them to the "normal" default values before a real run.
###################################################################################################################################################################

chenlin_mutation = int(sys.argv[1])
rosetta_scripts_path = os.path.expanduser("/usr/local/rosetta_mpi/main/source/bin/cartesian_ddg.mpi.macosclangrelease")
residue_to_mutate = ('A', chenlin_mutation, '') # Residue position to perfrom saturation mutatagenesis. Format: (Chain, PDB residue number, insertion code).

if not os.path.isfile(rosetta_scripts_path):
    print('ERROR: "rosetta_scripts_path" variable must be set to the location of the "rosetta_scripts" binary executable')
    print('This file might look something like: "cartesian_ddg.mpi.macosclangreleas"')
    raise Exception('Rosetta cartesian_ddg missing')

def run_Cartesian_ddg_saturation( name, input_path, input_pdb_path, mut_aa ):
    output_directory = os.path.join( 'output_saturation', os.path.join( '%s_%s' % (name, mut_aa) ) )
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    mutation_chain, mutation_resi, mutation_icode = residue_to_mutate
    mutfile_path = os.path.join( output_directory, 'mutate_%s%d%s_to_%s.mutfile' % (mutation_chain, mutation_resi, mutation_icode, mut_aa) )

    traj = md.load(input_pdb_path)
    original_residue = traj.topology.residue(mutation_resi-1).code

    with open( mutfile_path, 'w') as f:
        f.write( f'total 1 \n1\n{original_residue} {mutation_resi} {mut_aa}' )

    Cartesian_ddG_args = [
        os.path.abspath(rosetta_scripts_path),
        "-s %s" % os.path.abspath(input_pdb_path),
        '-ddg:iterations 3',
        '-ddg::cartesian',
        '-ddg::dump_pdbs False',
        '-bbnbrs 1',
        '-fa_max_dis 9.0',
        '-score:weights ref2015_cart ',
        '-relax:cartesian',
        '-relax:min_type lbfgs_armijo_nonmonotone',
        f'-ddg::mut_file {os.path.abspath( mutfile_path )} ',
        '-use_input_sc',
        '-flip_HNQ',
        '-optimization:default_max_cycles 200',
        '-crystal_refine ',
        '-ignore_zero_occupancy false',
        '-ex1',
        '-ex2',
    ]

    log_path = os.path.join(output_directory, 'rosetta.out')

    print( 'Running Rosetta with args:' )
    print( ' '.join(Cartesian_ddG_args) )
    print( 'Output logged to:', os.path.abspath(log_path) )
    print()

    outfile = open(log_path, 'w')
    process = subprocess.Popen(Cartesian_ddG_args, stdout=outfile, stderr=subprocess.STDOUT, close_fds = True, cwd = output_directory)
    returncode = process.wait()
    outfile.close()

if __name__ == '__main__':
    mutation_chain, mutation_resi, mutation_icode = residue_to_mutate
    cases = []
    # for nstruct_i in range(1, nstruct + 1 ):
    for case_name in os.listdir('inputs'):
        case_path = os.path.join( 'inputs', case_name )
        for f in os.listdir(case_path):
            if f.endswith('.pdb'):
                input_pdb_path = os.path.join( case_path, f )
                break


        for mut_aa in 'ACDEFGHIKLMNPQRSTVWY':
            cases.append( ('%s_%s%d%s' % (case_name, mutation_chain, mutation_resi, mutation_icode), case_path, input_pdb_path, mut_aa) )

    if use_multiprocessing:
        pool = multiprocessing.Pool( processes = min(max_cpus, multiprocessing.cpu_count()) )

    for args in cases:
        if use_multiprocessing:
            pool.apply_async( run_Cartesian_ddg_saturation, args = args )
        else:
            run_Cartesian_ddg_saturation( *args )

    if use_multiprocessing:
        pool.close()
        pool.join()
