from modeller import *

env = environ()
aln = alignment(env)
mdl = model(env, file='6q2f_temp', model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='6q2f_temp', atom_files='6q2f_temp.pdb')
aln.append(file='target.ali', align_codes='target')
aln.align2d()
aln.write(file='temp-target.ali', alignment_format='PIR')
aln.write(file='temp-target.pap', alignment_format='PAP')
