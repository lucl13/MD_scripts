#!/usr/bin/env python
# coding: utf-8


import os

import pandas as pd
import pyrosetta
import multiprocessing
import time

def runBackrub(pose, job_id):

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())

    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    mm.set_jump(True)

    OUTPUT_DIR = './output_3/'

    backrub = pyrosetta.rosetta.protocols.backrub.BackrubProtocol()
    backrub.set_taskfactory(tf)
    backrub.set_movemap(mm)
    backrub.apply(pose)

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    pose.dump_pdb(f'{OUTPUT_DIR}/1_apo_1urp_backrubed_{int(job_id):05}.pdb')

    score_df = pd.DataFrame()
    score_dict = dict(pose.scores)
    score_dict['description'] = f'backrub_{int(job_id):05}'
    score_df = score_df.append(score_dict, ignore_index=True)
    score_scv = f'{OUTPUT_DIR}/score_backub.csv'
    if os.path.isfile(score_scv):
        score_df.to_csv(score_scv, mode='a', header=False, index=False)
    else:
        score_df.to_csv(score_scv, index=False)


def runBackrubs(file, job_ids):
    pyrosetta.init('-backrub:ntrials 10000', '-out:overwrite')
    pose = pyrosetta.pose_from_pdb(file)
    for job_id in job_ids:
        runBackrub(pose, job_id)


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


    file ='1_apo_1urp_relaxed.pdb'
    n_sctrut = 10000

    job_ids = [i for i in range(1, n_sctrut+1)]
    chunks = make_chunks(job_ids, 22)
    threads = []

    for thread_num in range(22):
        current_thread = multiprocessing.Process(target=runBackrubs, args=(file, chunks[thread_num]))
        threads.append(current_thread)
        current_thread.start()
    for t in threads:
        t.join()

