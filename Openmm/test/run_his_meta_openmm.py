from his_meta_openmm import *

#from chenlin_openmm import * 

arg_md = ['--topology','1_168h169c_4mer_topol.top','1_168h169c_4mer_em.gro']
arg_meta = ['--residue','168']

cst_runtime = u.picoseconds*2
eq_runtime = u.picoseconds*2
pro_runtime = u.picoseconds*1
meta_runtime = u.picoseconds*2

parser = make_parser()
meta_parser = meta_make_parser()
opts = parser.parse_args(args=arg_md)
meta_opts = meta_parser.parse_args(args=arg_meta)

prepare_logging(opts.log, opts.debug)
mdtraj_topol = md.load(opts.gro_file[:-4]+'.pdb').topology
if opts.subset.lower() == "none":
    subset = None
else:
    subset = mdtraj_topol.select(opts.subset)

for run_num in range(opts.start_run, opts.start_run + opts.nruns):
    if opts.start_run == 0 and opts.nruns == 1:
        output_index = None
    else:
        output_index = run_num

    if output_index is not None:
        logger.info("Starting run " + str(output_index))
        
    name = re.sub(r'\.gro$', '', opts.gro_file)


    hismeta = HisMetadynamics(
        gro_file = opts.gro_file,
        top_file = opts.topology,
        subset = subset,
        temperature = opts.temperature*u.kelvin,
        pressure = opts.pressure*u.atmosphere,
        timestep = opts.timestep*u.femtoseconds,
        time_per_frame = opts.time_per_frame*u.picoseconds,
        residue = meta_opts.residue, 
        biasFactor = meta_opts.biasFactor,
        height = meta_opts.height*u.kilojoule_per_mole, 
        frequency = meta_opts.frequency,
        output_index = output_index,
        output_basename = name)

    #hismeta.run_cst(opts.gro_file, cst_runtime)
    #hismeta.run_eq(eq_runtime)
    #hismeta.run_pro(pro_runtime)
    hismeta.read_chk('eq')
    
    hismeta.run_meta(meta_runtime)
    print('Finished!')

