#!/usr/bin/env python
# coding: utf-8


import os

import pandas as pd
import pyrosetta
import multiprocessing
import time
import sys
import pyrosetta_help as ph
import re



def runBackrub(file, job_id):

    job_name = f'{file[:-4]}_{int(job_id):05}'
    logger = ph.configure_logger()
    pyrosetta.init('-backrub:ntrials 100 -out:overwrite true -extra_res_fa RIP.params -ignore_zero_occupancy false')

    pose = pyrosetta.pose_from_pdb(file)
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())

    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    mm.set_jump(True)

    backrub = pyrosetta.rosetta.protocols.backrub.BackrubProtocol()
    backrub.set_taskfactory(tf)
    backrub.set_movemap(mm)
    backrub.apply(pose)

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    pose.dump_pdb(f'{OUTPUT_DIR}/{job_name}.pdb')

    score_dict = dict(pose.scores)
    score_dict['description'] = f'backrub_{int(job_id):05}'
    score_df = score_df = pd.DataFrame(score_dict, index=[0])
    score_scv = f'{OUTPUT_DIR}/score_backub.csv'
    if os.path.isfile(score_scv):
        score_df.to_csv(score_scv, mode='a', header=False, index=False)
    else:
        score_df.to_csv(score_scv, index=False)

    entries = ph.get_all_log_entries()
    log_text = [str(e['text']) for e in entries]

    with open(f'{OUTPUT_DIR}/{job_name}.log', 'w') as f:
        for line in log_text:
            ansi_escape = re.compile(r'\x1b(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
            line_cleaned = ansi_escape.sub('', line)
            f.write(f"{line_cleaned}\n")


def runBackrubs(file, job_ids):
    for job_id in job_ids:
        #ph.get_all_log_entries()

        runBackrub(file, job_id)




def make_chunks(data: list, thread_count) -> dict:
    """
    Takes a list and splits it into parts based on the thread count
    :param data: a list that needs to be split up amongst the threads
    :param thread_count: the number of threads to use
    :return: None
    """
    threads = {}
    for x in range(0, thread_count):
        threads[x] = []
    thread = 0
    for x in range(0, len(data)):
        threads[thread].append(data[x])
        thread += 1
        if thread == thread_count:
            thread = 0
    return threads




if __name__ == '__main__':


    file ='2dri-RIP_relaxed.pdb'
    n_sctrut = 22
    
    global OUTPUT_DIR 
    OUTPUT_DIR = './output_2/'

    job_ids = [i for i in range(1, n_sctrut+1)]
    chunks = make_chunks(job_ids, 22)
    threads = []

    for thread_num in range(22):
        current_thread = multiprocessing.Process(target=runBackrubs, args=(file, chunks[thread_num]))
        threads.append(current_thread)
        current_thread.start()
    for t in threads:
        t.join()

