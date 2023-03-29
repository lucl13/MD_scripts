import numpy as np
import pandas as pd
import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from md_utility.utils import get_res_info


# parse the raw_hbonds

def parse_hbond_results(case_name, u: MDAnalysis.Universe, raw_hbonds_reslus):
    columns = ['frame', 'donor_index', 'hydrogen_index', 'acceptor_index', 'distance', 'angle']
    df = pd.DataFrame(raw_hbonds_reslus, columns=columns)
    df['donor'] = get_res_info(case_name, u, df['donor_index'])
    df['hydrogen'] = get_res_info(case_name, u, df['hydrogen_index'])
    df['acceptor'] = get_res_info(case_name, u, df['acceptor_index'])

    new_columns = ['frame', 'donor_index', 'hydrogen_index', 'acceptor_index', 'donor', 'hydrogen', 'acceptor',
                   'distance', 'angle']
    return df[new_columns]


# group the hbonds

# df_hbond.groupby(['donor', 'hydrogen', 'acceptor']).size()
# repeat or greater than hbonds.n_fames from the tetramer
# df_hbond[(df_hbond['donor'] == 'VAL85-N') & (df_hbond['hydrogen'] == 'VAL85-HN' )]

def group_hbonds(df_hbond: pd.DataFrame):
    df_hb_grouped = df_hbond.groupby(['donor', 'hydrogen', 'acceptor']).size().reset_index()
    n_frames = len(set(df_hbond.iloc[:, 0]))
    df_hb_grouped[0] = df_hb_grouped[0] / n_frames / 4 * 100  # average of the tetramer
    df_hb_grouped.columns = ['donor', 'hydrogen', 'acceptor', 'hbond_ratio']
    max_per_res = df_hb_grouped.groupby('donor').max()
    return df_hb_grouped, max_per_res


def run_hbond_analysis(case_name, u: MDAnalysis.Universe):
    data_file = './hbond/data/' + case_name + '_hbond.npy'
    try:
        hbonds_results = np.load(data_file)
        print(f'{case_name}_hbond.npy file loaded')
    except:
        hbonds = HBA(universe=u, d_a_cutoff=3.5)
        hbonds.donors_sel = 'protein and name N'
        hbonds.hydrogens_sel = 'protein and name HN'
        hbonds.acceptors_sel = 'protein and backbone and not name H*'
        hbonds.run(start=0)

        # save file
        np.save(data_file, hbonds.results.hbonds)
        print(f'{case_name}_hbond.npy file saved')

        hbonds_results = hbonds.results.hbonds
    df_hbond = parse_hbond_results(case_name, u, hbonds_results)
    df_hb_grouped, max_per_res = group_hbonds(df_hbond)


    return max_per_res
