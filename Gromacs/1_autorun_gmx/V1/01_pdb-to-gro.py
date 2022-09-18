#!/usr/bin/python3             

import subprocess    

for i in range(1,51):          
        #subprocess.run(' sed -i  -e  "/^HETATM/d"  {}.pdb '.format(i), shell=True)    
        #subprocess.run('grep   -v  HOH   {}.pdb > {}_clean.pdb'.format(i,i), shell=True)                           
        subprocess.run('echo 9 |gmx pdb2gmx  -f {}_PETase-MHETase.pdb  -o {}_PETase-MHETase.gro  -water tip3p  -p {}_topol.top'.format(i,i,i), shell=True)


exit()


