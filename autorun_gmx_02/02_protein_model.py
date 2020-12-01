#!/usr/bin/python3             
                                    ###   使用说明     ###
                   ##   用这个脚本之前，需要注意以下几点：                                ##
                   ##   1.不同体系修改name,n,num,T，name_p；                            ##
                   ##   2.修改下面代码行带## ！！的5个部分；                             ##
                   ##   3.手动修改mdp中的define,P，steps,energygrp（md.mdp）。          ##



import subprocess    
name=['PETase-MHETase']        


for j in name:         
    for i in range(1,51):    
        subprocess.run(f'cp topol.top {i}_topol.top', shell=True) 
        subprocess.run('gmx editconf -f {}_{}.gro  -o {}_box.gro    -c   -d 1.2'.format(i,j,i), shell=True)                             ##. !!修改02 盒子大小
        subprocess.run('gmx solvate -cp  {}_box.gro  -cs spc216.gro -p  {}_topol.top  -o  {}_sol.gro '.format(i,i,i), shell=True)
        subprocess.run('gmx grompp  -f ions.mdp    -c {}_sol.gro  -p  {}_topol.top  -o {}_ions.tpr  -maxwarn  1  '.format(i,i,i), shell=True) 
        subprocess.run('echo "sol\n" | gmx genion -s {}_ions.tpr  -o {}_ions.gro  -neutral  -p  {}_topol.top'.format(i,i,i), shell=True) 
        subprocess.run('gmx grompp  -f mini.mdp -c {}_ions.gro  -p {}_topol.top  -o {}_em.tpr   -maxwarn  1'.format(i,i,i), shell=True) 
        subprocess.run('gmx mdrun -v -deffnm {}_em'.format(i), shell=True) 
#        subprocess.run("echo  '\"{}\" |  13  \n q'   | gmx make_ndx -f {}_em.gro -o {}_index.ndx".format(n[0],i,i),shell=True)  ## 修改04  名称
#       subprocess.run('sed -i  ""  /^tc-grps/s/"Protein_lig Water_and_ions"/"{}_{} Water_and_ions"/   {}_md.mdp'.format(n[0],n[1],i),shell=True)     ## 修改05  名称
#      subprocess.run('sed -i  ""  /^tc-grps/s/"Protein_lig Water_and_ions"/"{}_{} Water_and_ions"/   npt.mdp'.format(n[0],n[1]),shell=True)     ## 修改05  名称
#        subprocess.run('gmx grompp -f npt.mdp -c {}_em.gro -r {}_em.gro  -n {}_index.ndx -p {}_topol.top -o {}_npt.tpr  -maxwarn  1'.format(i,i,i,i,i),shell=True)


exit()


