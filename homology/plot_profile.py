import  pylab
import modeller
import matplotlib.pyplot as plt
import numpy as np 
#matplotlib.rcParams["xtick.direction"]="inout"
#matplotlib.rcParams["ytick.direction"]="inout"
#matplotlib.rcParams["lines.color"]="red"
#plt.style.use("chenlin")
#matplotlib.rcParams['font.sans-serif'] = ['Arial']




def r_enumerate(seq):
    """Enumerate a sequence in reverse order"""
    # Note that we don't use reversed() since Python 2.3 doesn't have it
    num = len(seq) - 1
    while num >= 0:
        yield num, seq[num]
        num -= 1

def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to
       the alignment sequence `seq`."""
    # Read all non-comment and non-blank lines from the file:
    f = open(profile_file)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in r_enumerate(seq.residues):
        for gap in range(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals

e = modeller.environ()
a = modeller.alignment(e, file='temp-target.ali')

template = get_profile('6q2f_temp.profile', a['6q2f_temp'])
model = get_profile('target.B99990005.profile', a['target'])

# Plot the template and model profiles in the same plot for comparison:
plt.figure(1, figsize=(10,8))
ax=plt.gca()
ax.spines['left'].set_linewidth(3)
ax.spines['right'].set_linewidth(3)
ax.spines['top'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3) 
plt.xlabel('Alignment position',fontsize=24)
plt.ylabel('DOPE per-residue score',fontsize=24)
plt.ylim(-0.1,0.08)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.plot(model, color='red', linewidth=2, label='Model')
plt.plot(template, color='green', linewidth=2, label='Template_1')
plt.legend()
plt.savefig('dope_profile.png', dpi=65)
