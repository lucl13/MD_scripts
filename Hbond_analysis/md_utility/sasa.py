import os

import MDAnalysis
import mdtraj as md
import numpy as np
import pandas as pd

from md_utility.utils import get_res_info
from md_utility.utils import is_protein


def get_sasa_res(traj, case_name):
    # try to load the npy file
    data_dir = './sasa/data/'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    data_file = data_dir + case_name + '_sasa_res.npy'

    try:
        sasa_res = np.load(data_file)
        print(f'{data_file} loaded')
    except:
        sasa_res = md.shrake_rupley(traj[1000:], mode='residue')

        # write to a numpy file
        np.save(data_file, sasa_res)
        print(f'{data_file} saved')

    return sasa_res


def parse_sasa(case_name, top_file, sasa_res):
    topology = md.load(top_file).topology
    u = MDAnalysis.Universe(top_file)
    sasa_avg_res = np.average(sasa_res, axis=0)
    n_protein_res = 0
    sasa_protein_res = []
    protein_res_list = []
    for res in topology.residues:
        if is_protein(res):
            n_protein_res += 1
            protein_res_list.append(res.index)
            sasa_protein_res.append(sasa_avg_res[res.index])
    one_chain_length = int(n_protein_res / 4)

    #sasa_df = pd.DataFrame({'sasa': np.array(sasa_protein_res),
    #                        'residue': np.array(
    #                            [get_res_info(case_name, u, u.residues[i].atoms[0].index, with_atom_name=False)
    #                             for i in protein_res_list])})

    sasa_df = pd.DataFrame({'sasa': np.array(sasa_protein_res),
                            'residue': np.array(
                                [u.residues[i].resnum  for i in protein_res_list])})

    from md_utility.utils import sort_res_info_str_df
    #sasa_mean_df = sort_res_info_str_df(sasa_df.groupby('residue').mean(), 'residue')
    sasa_mean_df = sasa_df.groupby('residue').mean()
    return sasa_df, sasa_mean_df
