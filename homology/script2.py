from modeller import *

env = environ()
aln = alignment(env)
for (pdb, chain) in (('tem1', 'A'), ('tem2', 'A'), ('tem3', 'M'),
                     ('tem4', 'A'), ('tem5', 'A'), ('tem6', 'E'), ('tem7', 'A')):
    m = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
aln.malign()
aln.malign3d()
aln.compare_structures()
aln.id_table(matrix_file='family.mat')
env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)
