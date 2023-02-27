#!/usr/bin/env python3

#=========================================================================================
# Peacemaker -- A Quantum Cluster Equilibrium Code.
#
# Copyright 2004-2006 Barbara Kirchner, University of Bonn
# Copyright 2007-2012 Barbara Kirchner, University of Leipzig
# Copyright 2013-2023 Barbara Kirchner, University of Bonn
#
# This file is part of Peacemaker.
#
# Peacemaker is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Peacemaker is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Peacemaker.  If not, see <http://www.gnu.org/licenses/>
#=========================================================================================

import sys
import os.path
import getopt
import numpy as np
from numpy import linalg as LA
import math


# VdW volumes taken from Bondi's compilation.
# Masses taken from http://www.csudh.edu/oliver/chemdata/atmass.htm.
elements = {
        "Ac"  :  [2.47,  227.02800],
        "Al"  :  [1.84,   26.98154],
        "Am"  :  [2.44,  243.00000],
        "Sb"  :  [2.06,  121.75700],
        "Ar"  :  [1.88,   39.94800],
        "As"  :  [1.85,   74.92159],
        "At"  :  [2.02,  210.00000],
        "Ba"  :  [2.68,  137.32700],
        "Bk"  :  [2.44,  247.00000],
        "Be"  :  [1.53,    9.01218],
        "Bi"  :  [2.07,  208.98037],
        "Bh"  :  [0.00,  262.00000],
        "B"   :  [1.92,   10.81100],
        "Br"  :  [1.85,   79.90400],
        "Cd"  :  [2.18,  112.41100],
        "Ca"  :  [2.31,   40.07800],
        "Cf"  :  [2.45,  251.00000],
        "C"   :  [1.70,   12.01100],
        "Ce"  :  [2.42,  140.11500],
        "Cs"  :  [3.43,  132.90543],
        "Cl"  :  [1.75,   35.45270],
        "Cr"  :  [2.06,   51.99610],
        "Co"  :  [2.00,   58.93320],
        "Cu"  :  [1.96,   63.54600],
        "Cm"  :  [2.45,  247.00000],
        "Db"  :  [0.00,  262.00000],
        "Dy"  :  [2.31,  162.50000],
        "Es"  :  [2.45,  252.00000],
        "Er"  :  [2.29,  167.26000],
        "Eu"  :  [2.35,  151.96500],
        "Fm"  :  [2.45,  257.00000],
        "F"   :  [1.47,   18.99840],
        "Fr"  :  [3.48,  223.00000],
        "Gd"  :  [2.34,  157.25000],
        "Ga"  :  [1.87,   69.72300],
        "Ge"  :  [2.11,   72.61000],
        "Au"  :  [2.14,  196.96654],
        "Hf"  :  [2.23,  178.49000],
        "Hs"  :  [0.00,  265.00000],
        "He"  :  [1.40,    4.00260],
        "Ho"  :  [2.30,  164.93032],
        "H"   :  [1.10,    1.00794],
        "In"  :  [1.93,  114.82000],
        "I"   :  [1.98,  126.90447],
        "Ir"  :  [2.13,  192.22000],
        "Fe"  :  [2.04,   55.84700],
        "Kr"  :  [2.02,   83.80000],
        "La"  :  [2.43,  138.90550],
        "Lr"  :  [2.46,  262.00000],
        "Pb"  :  [2.02,  207.20000],
        "Li"  :  [1.82,    6.94100],
        "Lu"  :  [2.24,  174.96700],
        "Mg"  :  [1.73,   24.30500],
        "Mn"  :  [2.05,   54.93805],
        "Mt"  :  [0.00,  266.00000],
        "Md"  :  [2.46,  258.00000],
        "Hg"  :  [2.23,  200.59000],
        "Mo"  :  [2.17,   95.94000],
        "Nd"  :  [2.39,  144.24000],
        "Ne"  :  [1.54,   20.17970],
        "Np"  :  [2.39,  237.04800],
        "Ni"  :  [1.97,   58.69340],
        "Nb"  :  [2.18,   92.90638],
        "N"   :  [1.55,   14.00674],
        "No"  :  [2.46,  259.00000],
        "Os"  :  [2.16,  190.20000],
        "O"   :  [1.52,   15.99940],
        "Pd"  :  [2.10,  106.42000],
        "P"   :  [1.80,   30.97376],
        "Pt"  :  [2.13,  195.08000],
        "Pu"  :  [2.43,  244.00000],
        "Po"  :  [1.97,  209.00000],
        "K"   :  [2.75,   39.09830],
        "Pr"  :  [2.40,  140.90765],
        "Pm"  :  [2.38,  145.00000],
        "Pa"  :  [2.43,  231.03590],
        "Ra"  :  [2.83,  226.02500],
        "Rn"  :  [2.20,  222.00000],
        "Re"  :  [2.16,  186.20700],
        "Rh"  :  [2.10,  102.90550],
        "Rb"  :  [3.03,   85.46780],
        "Ru"  :  [2.13,  101.07000],
        "Rf"  :  [0.00,  261.00000],
        "Sm"  :  [2.36,  150.36000],
        "Sc"  :  [2.15,   44.95591],
        "Sg"  :  [0.00,  263.00000],
        "Se"  :  [1.90,   78.96000],
        "Si"  :  [2.10,   28.08550],
        "Ag"  :  [2.11,  107.86820],
        "Na"  :  [2.27,   22.98977],
        "Sr"  :  [2.49,   87.62000],
        "S"   :  [1.80,   32.06600],
        "Ta"  :  [2.22,  180.94790],
        "Tc"  :  [2.16,   98.00000],
        "Te"  :  [2.06,  127.60000],
        "Tb"  :  [2.33,  158.92534],
        "Tl"  :  [1.96,  204.38330],
        "Th"  :  [2.45,  232.03810],
        "Tm"  :  [2.27,  168.93421],
        "Sn"  :  [2.17,  118.71000],
        "Ti"  :  [2.11,   47.88000],
        "W"   :  [2.18,  183.85000],
        "U"   :  [2.41,  238.02890],
        "V"   :  [2.07,   50.94150],
        "Xe"  :  [2.16,  131.29000],
        "Yb"  :  [2.26,  173.04000],
        "Y"   :  [2.32,   88.90585],
        "Zn"  :  [2.01,   65.39000],
        "Zr"  :  [2.23,   91.22400],
        "Xx"  :  [0.00,    0.00000]
        }            



options = {'-r': False, '-f': False, '-x': False,
           '-t': False, '-v': False, '-m': False}

try:
    opts, args = getopt.getopt(sys.argv[1:], "frxtvm")
except getopt.GetoptError:
    print("Usage:\t{} [option] <files>\n\t[-r] \t rotates molecule with inertia tensor\n\t[-x] \t creates a xyz format output\n\t[-t] \t creates a tex format output".format(sys.argv[0]))
    sys.exit(2)
for opt, arg in opts:
  if opt in options:
    options[opt] = True

if len(args) == 0:
    print("Usage:\t{} [option] <files>\n\t[-r] \t rotates molecule with inertia tensor\n\t[-x] \t creates a xyz format output\n\t[-t] \t creates a tex format output".format(sys.argv[0]))
    sys.exit(1)

nX = 0
W = np.zeros((len(args), 3))

for j, filename in enumerate(args):
  filename = args[j]

  with open(filename, "r") as fin:
    # read xyz file
    line = fin.readline()
    nAtoms = int(line)
    pos = np.zeros((nAtoms, 3))
    elSymbols = []
    m = np.zeros(nAtoms)
    vol = np.zeros(nAtoms)
    line = fin.readline()
    comment = line[:-1]
    for i in range(nAtoms):
      line = fin.readline()
      symbol, x, y, z = line.split()
      symbol = symbol.lower().capitalize()
      if symbol == "Xx": nX = nX + 1
      m[i] = elements[symbol][1]
      vol[i] = 4.0/3.0*math.pi*elements[symbol][0]**3
      elSymbols.append(symbol)
      pos[i,:] = float(x), float(y), float(z)

    formula = ""
    for e in elements:
      count = elSymbols.count(e)
      if (count > 0):
        formula = formula + e
        if (count > 1):
          formula = formula + str(count)

    # shift COM to origin
    com = np.zeros(3)
    mtot = 0.0
    vtot = 0.0
    for i in range(nAtoms):
      com += m[i]*pos[i,:]
      mtot += m[i]
      vtot += vol[i]

    if (options['-v']):
      print(vtot)
      continue
    elif (options['-m']):
      print(mtot)
      continue
    elif (options['-f']):
      print(formula)
      continue

    com /= mtot
    for i in range(nAtoms):
      pos[i,:] -= com

    # calculate inertia tensor
    inertia = np.zeros((3, 3))
    for i in range(nAtoms):
      inertia[0, 0] += m[i]*(pos[i,1]**2 + pos[i,2]**2)
      inertia[1, 1] += m[i]*(pos[i,0]**2 + pos[i,2]**2)
      inertia[2, 2] += m[i]*(pos[i,0]**2 + pos[i,1]**2)
      inertia[0, 1] -= m[i]*pos[i,0]*pos[i,1]
      inertia[0, 2] -= m[i]*pos[i,0]*pos[i,2]
      inertia[1, 2] -= m[i]*pos[i,1]*pos[i,2]
    inertia[1,0] = inertia[0,1]
    inertia[2,0] = inertia[0,2]
    inertia[2,1] = inertia[1,2]

    # diagonalize inertia tensor and print
    w, v = LA.eig(inertia)

    p = w.argsort()
    w = w[p]
    v = v[:,p]

#    print(v)

    # if -r option is active: rotate coordinates with inertia tensor
    if (options['-r']):
      for i in range(nAtoms):
        pos[i,:] = v.T.dot(pos[i,:])

    # if -x option is active: print out in .xyz format
    if (options['-x']):
      print("%u" % (nAtoms-nX))
      print("{:s}".format(comment))
      for i in range(nAtoms):
        if elSymbols[i]=="CL": elSymbols[i]="Cl"
        if elSymbols[i]=="Xx": continue
        print("{:3s} {:11.8f} {:13.8f} {:13.8f}".format(elSymbols[i], pos[i,0], pos[i,1], pos[i,2]))
      print("")

    # if -t option is active: print out in tex format
    elif (options['-t']):
      print("\\begin{longtable}[l]{lrrr}")
      print("\multicolumn{{4}}{{l}}{{\\textbf{{{:s}}}}}\\\\".format(filename))
      print("%u&&&\\\\" % (nAtoms-nX))
      print("\multicolumn{{4}}{{l}}{{{:s}}}\\\\".format(comment))

      for i in range(nAtoms):
        if elSymbols[i]=="XX": continue
        print("{:3s} & {:11.8f} & {:13.8f} & {:13.8f} \\\\".format(elSymbols[i], pos[i,0], pos[i,1], pos[i,2]))
      print("\\end{longtable}\n")


    # default, if no option was specified
    else:
      asym = (2 * 1/w[1] - 1/w[0] - 1/w[2])/(1/w[0] - 1/w[2])

      print(formula)
      print("volume {:.4f}".format(vtot))
      print("mass {:.4f}".format(mtot))
      print("inertia {:.4f} {:.4f} {:.4f}".format(w[0], w[1], w[2]))
      print("asymmetry parameter {:.4f}\n".format(asym))
    W[j][0] = w[0]
    W[j][1] = w[1]
    W[j][2] = w[2]
  fin.close()

# if two input files where specified, calculate geometrical distance
if len(args) == 2 and not any(options):
  gdist = pow(pow((1/W[0][0]-1/W[1][0])/(1/W[0][0]), 2) + pow((1/W[0][1]-1/W[1][1])/(1/W[0][1]), 2) + pow((1/W[0][2]-1/W[1][2])/(1/W[0][2]), 2), 0.5)
  print("geometrical distance: {:.4f}".format(gdist))
