#!/usr/bin/python3

import subprocess

names = ['apo_2PFK']   

for name in names:
	for i in range(1,2):
		subprocess.run(f'echo "protein system"| gmx trjconv -f {i}_{name}_md.xtc -s {i}_{name}_md.tpr -o {i}_{name}_cluster.xtc -pbc cluster', shell=True)
		subprocess.run(f'echo "protein system" | gmx trjconv -f {i}_{name}_cluster.xtc -s {i}_{name}_md.tpr -o {i}_{name}_fit.xtc -fit rot+trans' , shell=True)
		subprocess.run(f'echo "non-water"| gmx trjconv -f {i}_{name}_fit.xtc -s {i}_{name}_md.tpr   -o {i}_{name}_nowat.xtc -pbc mol', shell=True)
		subprocess.run(f'echo "non-water" | gmx trjconv -f {i}_{name}_em.gro   -s {i}_{name}_md.tpr  -o {i}_{name}_nowat.gro', shell=True) 
		
		
