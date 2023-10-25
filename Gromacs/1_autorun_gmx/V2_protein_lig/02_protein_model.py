#!/usr/bin/python3             

import subprocess    

names = ['4XYK_17-784-FLC']         


for name in names:         
    for i in range(1,2):    
        subprocess.run(f'cp topol.top {i}_{name}_topol.top', shell=True) 
        subprocess.run(f'echo "protein" | gmx editconf -f {i}_{name}.pdb -o {i}_{name}_box.pdb -c -d 1.2', shell=True)                             ##. !!修改02 盒子大小
        subprocess.run(f'gmx solvate -cp {i}_{name}_box.pdb -cs spc216.gro -p {i}_{name}_topol.top -o {i}_{name}_sol.pdb', shell=True)
        subprocess.run(f'gmx grompp -f ions.mdp -c {i}_{name}_sol.pdb -p {i}_{name}_topol.top -o {i}_{name}_ions.tpr -maxwarn 1', shell=True) 
        subprocess.run(f'echo "sol\n" | gmx genion -s {i}_{name}_ions.tpr -o {i}_{name}_ions.pdb -neutral -conc 0.15 -p {i}_{name}_topol.top', shell=True) 
       # subprocess.run(f"sed -ie 's/CL/CLA/g' {i}_{name}_topol.top".format(i), shell=True)
        subprocess.run(f'gmx grompp -f mini.mdp -c {i}_{name}_ions.pdb  -p {i}_{name}_topol.top -o {i}_{name}_em.tpr -maxwarn 1', shell=True) 
        subprocess.run(f'gmx mdrun -v -deffnm {i}_{name}_em'.format(i), shell=True) 

exit()


