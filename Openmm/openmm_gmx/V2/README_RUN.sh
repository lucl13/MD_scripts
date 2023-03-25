#!/bin/csh
#

# Generated by CHARMM-GUI (http://www.charmm-gui.org) v3.7
#
# This folder contains a pre-optimized PDB structure and OpenMM inputs.
# All input files were optimized for OpenMM v6.2 or above, so lower version of OpenMM can cause some errors.
# You can get the latest development version of OpenMM at the git repository:
# https://github.com/pandegroup/openmm

set init = step3_input
set equi_prefix = step4_equilibration
set prod_prefix = step5_production
set prod_step   = step5
set num_rep = 3

# Equilibration
#set input_param = "-t 4_apo_4XYJ_17-784_topol.top -g 4_apo_4XYJ_17-784_md.gro"
#python -u openmm_run.py -i ${equi_prefix}.inp ${input_param} -orst ${equi_prefix}.rst -odcd ${equi_prefix}.dcd > ${equi_prefix}.out


# For REUS, please prepare initial configuration for each replica and use instead
set num_rep_m1 = `expr ${num_rep} - 1`
foreach x (`(seq 0 1 ${num_rep_m1})`)
    cp ${equi_prefix}.rst ${equi_prefix}_${x}.rst
end

# Production
# The OpenMM check point file (.chk) cannot be used in a different machine environment.
# So please make sure if you are using the same GPU and CUDA version of machine while doing additional
# production steps with the check point file.
set cnt = 1
set cntmax = 3

while ( ${cnt} <= ${cntmax} )
    @ pcnt = ${cnt} - 1
    set istep = ${prod_step}_${cnt}
    set pstep = ${prod_step}_${pcnt}
    if ( ${cnt} == 1 ) then
        set pstep = ${equi_prefix}
    endif
    set input_param = "-t 4_apo_4XYJ_17-784_topol.top -g 4_apo_4XYJ_17-784_md.gro -irst ${pstep}"
    python -u openmm_tremd.py --round 100 --steps 5000 -i ${prod_prefix}.inp ${input_param} -orst ${istep} -odcd ${istep} > ${istep}.out
    @ cnt += 1
end