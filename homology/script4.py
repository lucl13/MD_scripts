from modeller import *
from modeller.automodel import *
#from modeller import soap_protein_od

env = environ()
a = automodel(env, alnfile='temp-target.ali',
              knowns=('6q2f_temp'), sequence='target',
              assess_methods=(assess.DOPE,assess.GA341))
a.starting_model = 1
a.ending_model = 10
a.make()
