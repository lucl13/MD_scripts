#!/usr/bin/python3

import subprocess
name="PETase-MHETase-md"

for i in range(1, 51):
	subprocess.run(f'echo "protein system"| gmx trjconv -f {i}_{name}.xtc -s {i}_{name}.tpr -o cluster_{i}.xtc -pbc cluster', shell=True)
	subprocess.run(f'echo "protein system" | gmx trjconv -f cluster_{i}.xtc -s {i}_{name}.tpr -o fit_{i}.xtc -fit rot+trans' , shell=True)
	subprocess.run(f'echo "non-water"| gmx trjconv -f fit_{i}.xtc -s {i}_{name}.tpr   -o {i}_{name}_nowat.xtc -pbc mol', shell=True)
	subprocess.run(f'echo "non-water" | gmx trjconv -f {i}_em.gro   -s {i}_{name}.tpr  -o {i}_{name}_nowat.gro', shell=True) 
	