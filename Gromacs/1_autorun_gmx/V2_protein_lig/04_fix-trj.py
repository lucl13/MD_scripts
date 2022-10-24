#!/usr/bin/python3

import subprocess

names = ['4XYK_17-784-FLC']   

for name in names:
	for i in range(1,2):
		subprocess.run(f'echo "protein_lig system"| gmx trjconv -f {i}_{name}_md.xtc -s {i}_{name}_md.tpr -o {i}_{name}_cluster.xtc -pbc cluster -n index.ndx', shell=True)
		subprocess.run(f'echo "protein_lig system" | gmx trjconv -f {i}_{name}_cluster.xtc -s {i}_{name}_md.tpr -o {i}_{name}_fit.xtc -fit rot+trans -n index.ndx' , shell=True)
		subprocess.run(f'echo "non-water"| gmx trjconv -f {i}_{name}_fit.xtc -s {i}_{name}_md.tpr   -o {i}_{name}_nowat.xtc -pbc mol', shell=True)
		subprocess.run(f'echo "non-water" | gmx trjconv -f {i}_{name}_em.gro   -s {i}_{name}_md.tpr  -o {i}_{name}_nowat.gro', shell=True) 
		
