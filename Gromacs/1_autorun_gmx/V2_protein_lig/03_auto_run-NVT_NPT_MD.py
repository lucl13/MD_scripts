#!/usr/bin/python3   

import subprocess 

md_time = 60 # ns
md_steps = int(md_time*1e6/2)  
   
names = ['4XYK_17-784-FLC']   

def if_job_finished(job):
    filename = job + '.log'
    try:
        tail_output = subprocess.check_output(['tail', filename])
        last_line = tail_output.strip().decode("utf-8").split('\n')[-1]
        return last_line.startswith('Finished')
    except:
        return False

def run_NVT(i, name):
    subprocess.run(f'gmx grompp -f nvt.mdp -c {i}_{name}_em.gro -r {i}_{name}_em.gro -p {i}_{name}_topol.top -o {i}_{name}_nvt.tpr -n {i}_{name}_index.ndx', shell=True) 
    p = subprocess.run(f'gmx mdrun -v -deffnm {i}_{name}_nvt', shell=True) 
    if p.returncode == 0:
        print(f'{i}_{name} NVT job Done!')
    else: 
        raise Exception(f'{i}_{name} NVT job Exception !!!') 

def run_NPT(i, name):
    subprocess.run(f'gmx grompp -f npt.mdp -c {i}_{name}_nvt.gro -r {i}_{name}_nvt.gro -t {i}_{name}_nvt.cpt -p {i}_{name}_topol.top -o {i}_{name}_npt.tpr -maxwarn 1 -n {i}_{name}_index.ndx', shell=True) 
    p = subprocess.run(f'gmx mdrun -v -deffnm {i}_{name}_npt', shell=True) 
    if p.returncode == 0:
        print(f'{i}_{name} NPT job Done!')
    else: 
        raise Exception(f'{i}_{name} NPT job Exception !!!') 

def run_MD(i, name, md_steps):
    subprocess.run(f"sed -ie 's/nsteps     = .*/nsteps     = {md_steps}/g' md.mdp", shell=True)
    subprocess.run(f'gmx grompp -f md.mdp -c {i}_{name}_npt.gro -r {i}_{name}_npt.gro -t {i}_{name}_npt.cpt -p {i}_{name}_topol.top -o {i}_{name}_md.tpr -n {i}_{name}_index.ndx', shell=True) 
    p = subprocess.run(f'gmx mdrun -v -deffnm {i}_{name}_md', shell=True) 
    if p.returncode == 0:
        print(f'{i}_{name} MD job Done!')
    else: 
        raise Exception(f'{i}_{name} MD job Exception !!!') 


for name in names:         
    for i in range(1,2):  

        if if_job_finished(f'{i}_{name}_md'):
            print('MD finished! Noting to do')

        else:
            if if_job_finished(f'{i}_{name}_npt'):

                print('Running MD .....')
                run_MD(i, name, md_steps)
            else:

                if if_job_finished(f'{i}_{name}_nvt'): 

                    print('Running NPT .....')
                    run_NPT(i, name)

                    print('Running MD .....')
                    run_MD(i, name, md_steps)

                else:

                    print('Running NVT .....')
                    run_NVT(i, name)

                    print('Running NPT .....')
                    run_NPT(i, name)

                    print('Running MD .....')
                    run_MD(i, name, md_steps)

exit()




