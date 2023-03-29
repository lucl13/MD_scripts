import MDAnalysis
import numpy as np
import pandas as pd


def is_protein(res):
    """
    replace the original mdtraj is_protein function,
    since it gives wrong return for special residues like HISP/ASPP
    :param res: residue
    :return: bool
    """
    atom_names = [i.name for i in res.atoms]

    if {'CA', 'C', 'N'}.issubset(set(atom_names)):
        return True
    else:
        return False


def get_res_info(case_name, u: MDAnalysis.Universe, atom_index, with_atom_name=True):
    """
    get the residue info give an atom index

    :param with_atom_name: bool
    :param case_name: str
    :param u: MDAnalysis.Universe
    :param atom_index:
    :return: residue info: format->ALA180-N
    """
    # atom_index: res_info dict
    data_file = f'./hbond/data/{case_name}_res-info_dict.npy'
    try:
        atom_index_to_res_info = np.load(data_file, allow_pickle=True).item()
    except:
        atom_index_to_res_info = {}
        for atom in u.atoms:
            resname = str(atom.resname)
            resnum = str(atom.resnum)
            atom_name = str(atom.name)
            atom_index_to_res_info[atom.index] = ''.join([resname, resnum, '-', atom_name])
        np.save(data_file, atom_index_to_res_info)

    if with_atom_name:
        if isinstance(atom_index, np.int64):
            res_info = atom_index_to_res_info[atom_index]
            return res_info
        else:
            res_info = pd.Series([atom_index_to_res_info[i] for i in atom_index])
            return res_info
    else:
        if isinstance(atom_index, np.int64):
            res_info = atom_index_to_res_info[atom_index].split('-')[0]
            return res_info
        else:
            res_info = pd.Series([atom_index_to_res_info[i].split('-')[0] for i in atom_index])
            return res_info


def plot_on_pdb(case_name: str, pdb: str, df_per_res_1: pd.DataFrame, column_name: str,
                df_per_res_2: pd.DataFrame = None,
                with_atom_name=True,
                save_pdb=False):
    # case_name = pdb_file.split('/')[-1][:-4]
    u = MDAnalysis.Universe(pdb)
    u.add_TopologyAttr('tempfactors')  # add empty attribute for all atoms
    protein = u.select_atoms('protein')  # select protein atoms

    if df_per_res_2 is None:
        for residue in protein.residues:
            try:
                per_res_1 = df_per_res_1.loc[get_res_info(case_name, u, residue.atoms[0].index, with_atom_name)][
                    column_name]
                residue.atoms.tempfactors = per_res_1
            except:
                residue.atoms.tempfactors = 0

    else:
        for residue in protein.residues:
            try:

                per_res_1 = df_per_res_1.loc[get_res_info(case_name, u, residue.atoms[0].index, with_atom_name)][
                    column_name]
                per_res_2 = df_per_res_2.loc[get_res_info(case_name, u, residue.atoms[0].index, with_atom_name)][
                    column_name]
                residue.atoms.tempfactors = per_res_1 - per_res_2
            except:
                residue.atoms.tempfactors = 0
    if save_pdb:
        u.atoms.write(f'./hbond/{case_name}_hbond.pdb')
    return u


def sort_res_info_str_df(df, index_name):
    import re
    """
    Sort the residue column by the residue num in residue info ('ALA33')
    :param df: Dataframe
    :column_name: str
    :return: Dataframe
    """
    df = df.reset_index()
    df['resnum'] = df[index_name].apply(lambda x: re.search('[A-Z]*(\d*).*', x).group(1))
    df = df.sort_values(by='resnum').drop(columns='resnum')
    df.reset_index(inplace=True, drop=True)
    return df
