from gmx_openmm import * 

arg_md = ['--topology','1_168h169c_4mer_topol.top','1_168h169c_4mer_em.gro']
#arg_meta = ['--residue','168']

cst_runtime = u.picoseconds*1
eq_runtime = u.picoseconds*1
pro_runtime = u.picoseconds*1
#meta_runtime = u.picoseconds*1

parser = make_parser()
#meta_parser = meta_make_parser()
opts = parser.parse_args(args=arg_md)
#meta_opts = meta_parser.parse_args(args=arg_meta)

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
    gmx_simu = GMXSimulation(
        gro_file=opts.gro_file,
        top_file=opts.topology,
        subset=subset,
        temperature=opts.temperature*u.kelvin,
        pressure=opts.pressure*u.atmosphere,
        timestep=opts.timestep*u.femtoseconds,
        time_per_frame=opts.time_per_frame*u.picoseconds,
        output_index=output_index,
        output_basename=name
    )
    
    #print('cst running')
    #gmx_simu.run_cst(opts.gro_file, cst_runtime)

    #print('eq running')
    #gmx_simu.run_eq(eq_runtime)

    #print('pro running')
    #gmx_simu.run_pro(pro_runtime)

    gmx_simu.read_chk('eq')
    gmx_simu.run_pro(pro_runtime)