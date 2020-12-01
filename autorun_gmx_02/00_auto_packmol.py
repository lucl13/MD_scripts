import subprocess
path='/Volumes/px-PETase/5.models/01.protein'
for i in range(1, 51):
    input_file = ['tolerance 2.0\n',
     'filetype pdb\n',
     'seed -1\n',
     f'output {i}_PETase-MHETase.pdb\n',
     f'structure {path}/MHETase-6jtu_C.pdb\n',
     '  number 1 \n',
     'center \n',
     f'  fixed 50 50 50 {i*7.2} {120+i*7.2} {240+i*7.2}\n',
     'chain C \n',
     'add_amber_ter \n',
     'end structure\n',
     f'structure {path}/PETase-5xjh.pdb\n',
     '  number 1 \n',
     'center \n',
     f'  fixed 50 50 120 {240+i*7.2} {120+i*7.2} {i*7.2}\n',
     'chain A \n',
     'add_amber_ter \n',
     'end structure\n']
    with open('rand_50.inp','w') as f:
        f.writelines(input_file)
        
    subprocess.run('packmol < rand_50.inp',shell=True)