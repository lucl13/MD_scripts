import MDAnalysis
import mdtraj
import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns

from md_utility.utils import is_protein


def get_rmsf_all_atom(traj: mdtraj.Trajectory, case_name: str):
    # try to load the npy file
    try:
        rmsf_all_atom = np.load('./data/' + case_name + '_rmsf_all_atom.npy')
        print(f'{case_name}_rmsf_all_atom.npy file loaded')
    except:
        rmsf_all_atom = md.rmsf(traj[1000:], traj[0])

        # write to a numpy file
        np.save('./data/' + case_name + '_rmsf_all_atom.npy', rmsf_all_atom)
        print(f'{case_name}_rmsf_all_atom.npy file saved')

    return rmsf_all_atom


def get_backbone_rmsf_byres(topology, rmsf_all_atom):
    n_protein_res = 0
    rmsf_by_res = []
    for res in topology.residues:
        if is_protein(res):
            n_protein_res += 1

            atom_index = topology.select(f'resid {res.index} and (name CA or name N or name C or name O)')
            rmsf_by_res.append(np.average(rmsf_all_atom[atom_index]))

    one_chain_legnth = int(n_protein_res / 4)
    rmsf_df = pd.DataFrame({'rmsf': np.array(rmsf_by_res),
                            'residue': np.array([j for i in range(4) for j in range(topology.residue(0).resSeq,
                                                                                    topology.residue(
                                                                                        0).resSeq + one_chain_legnth)])})

    rmsf_mean_df = rmsf_df.groupby('residue').mean()
    return rmsf_df, rmsf_mean_df


def plot_rmsf(rmsf_df, rmsf_mean_df, ax):
    sns.lineplot(x='residue', y='rmsf', data=rmsf_df, ax=ax)
    # sns.scatterplot(x='residue', y='rmsf', data=rmsf_mean_df[rmsf_mean_df['rmsf'] > 0.3], ax=ax)

    # text the residue number of the residues with high RMSF
    # for i in rmsf_mean_df[rmsf_mean_df['rmsf'] > 0.3].index:
    # ax.text(i, rmsf_mean_df.loc[i, 'rmsf'], i, fontsize=12)


def rmsd_diff(rmsf_mean_1, rmsf_mean_2, pdb_file):
    diff_mean = rmsf_mean_1['rmsf'] - rmsf_mean_2['rmsf']
    diff_mean.apply(abs).sort_values(ascending=False).head(20)

    # save_pdb_with_bfactor

    u = MDAnalysis.Universe(pdb_file)
    u.add_TopologyAttr('tempfactors')  # add empty attribute for all atoms
    protein = u.select_atoms('protein')  # select protein atoms

    for residue, r_value in zip(protein.residues, pd.concat([diff_mean] * 4)):
        residue.atoms.tempfactors = r_value
    u.atoms.write('./rmsf/2pfk_rmsf-diff_complex-apo.pdb')

    return diff_mean
