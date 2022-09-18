import numpy as np 

import matplotlib.pyplot as plt 

def parse_line(raw_line, delimiter=' '):
    line = []
    for i in raw_line.strip().split(delimiter):
        if i != '':
            line.append(i)
    return line

def loadxvg(name):
    f = open(name)
    raw_data = f.readlines()
    f.close()
    data = []
    for i in raw_data:
        line = parse_line(i)
        if (line[0] == "#" or line[0]  == "@" or line[0]== "@TYPE") and len(line)>=2:
            if line[1] == "xaxis":
                x_label = line[3].strip('"')
            elif line[1] == "yaxis":
                y_label = line[3].strip('"()"')
            else:
                continue
        elif len(line)>=2 :
            data.append([float(line[0]), float(line[1])])
    return np.array(data), x_label, y_label


plt.style.use('chenlin')

data, x_label, y_label = loadxvg('rmsf_1.xvg')
plt.plot(data[:,0], data[:,1])
plt.xlabel(x_label)
plt.ylabel(y_label)

#plt.savefig('rmsd.png')
#plt.show()
