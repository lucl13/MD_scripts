"""
Generated by CHARMM-GUI (http://www.charmm-gui.org)

omm_readparams.py

This module is for reading coordinates and parameters in OpenMM.

Correspondance: jul316@lehigh.edu or wonpil@lehigh.edu
Last update: June 18, 2021
"""

import os, json

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

def read_top(filename, grofile, fftype='GROMACS'):
    gro = GromacsGroFile(grofile)
    box_vectors = Vec3(*gro.getPeriodicBoxVectors())
    top = GromacsTopFile(filename, periodicBoxVectors=box_vectors)
    return top

def read_crd(filename, fftype='GROMACS'):
    crd = GromacsGroFile(filename)
    return crd

def read_charmm_rst(filename):
    charmm_rst = CharmmRstFile(filename)

    for i, line in enumerate(charmm_rst.header):
        line = line.strip()
        words = line.split()
        if len(line) != 0:
            if words[0] == 'CRYSTAL' or words[0] == '!CRYSTAL':
                line1 = charmm_rst.header[i+1]
                line2 = charmm_rst.header[i+2]
                boxlx = Vec3(float(line1.split()[0].replace("D", "E")), 0.0, 0.0)
                boxly = Vec3(0.0, float(line1.split()[2].replace("D", "E")), 0.0)
                boxlz = Vec3(0.0, 0.0, float(line2.split()[2].replace("D", "E")))
                box = (boxlx, boxly, boxlz)
                break

    positions = charmm_rst.positionsold
    new_positions = []

    for position in positions:
        oldx = position[0]/angstrom
        oldy = position[1]/angstrom
        oldz = position[2]/angstrom

        newx = oldx + boxlx[0]/2.0
        newy = oldy + boxly[1]/2.0
        newz = oldz + boxlz[2]/2.0

        new_positions.append(Vec3(newx, newy, newz))

    charmm_rst.box = Quantity(box, angstroms)
    charmm_rst.positions = Quantity(new_positions, angstroms)

    return charmm_rst

def read_params(filename):
    parFiles = ()
    for line in open(filename, 'r'):
        if '!' in line: line = line.split('!')[0]
        parfile = line.strip()
        if len(parfile) != 0: parFiles += ( parfile, )

    params = CharmmParameterSet( *parFiles )
    return params

