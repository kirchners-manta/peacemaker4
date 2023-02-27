#!/usr/bin/env python

import sys
import os
import re
import numpy as np
import math
import getopt
import itertools
import shutil
from numpy import linalg as LA


# Masses taken from http://www.csudh.edu/oliver/chemdata/atmass.htm.
masses = {
        "Ac" : 227.028,
        "Al" : 26.981539,
        "Am" : 243,
        "Sb" : 121.757,
        "Ar" : 39.948,
        "As" : 74.92159,
        "At" : 210,
        "Ba" : 137.327,
        "Bk" : 247,
        "Be" : 9.012182,
        "Bi" : 208.98037,
        "Bh" : 262,
        "B"  : 10.811,
        "Br" : 79.904,
        "Cd" : 112.411,
        "Ca" : 40.078,
        "Cf" : 251,
        "C"  : 12.011,
        "Ce" : 140.115,
        "Cs" : 132.90543,
        "Cl" : 35.4527,
        "Cr" : 51.9961,
        "Co" : 58.93320,
        "Cu" : 63.546,
        "Cm" : 247,
        "Db" : 262,
        "Dy" : 162.50,
        "Es" : 252,
        "Er" : 167.26,
        "Eu" : 151.965,
        "Fm" : 257,
        "F"  : 18.9984032,
        "Fr" : 223,
        "Gd" : 157.25,
        "Ga" : 69.723,
        "Ge" : 72.61,
        "Au" : 196.96654,
        "Hf" : 178.49,
        "Hs" : 265,
        "He" : 4.002602,
        "Ho" : 164.93032,
        "H"  : 1.00794,
        "In" : 114.82,
        "I"  : 126.90447,
        "Ir" : 192.22,
        "Fe" : 55.847,
        "Kr" : 83.80,
        "La" : 138.9055,
        "Lr" : 262,
        "Pb" : 207.2,
        "Li" : 6.941,
        "Lu" : 174.967,
        "Mg" : 24.3050,
        "Mn" : 54.93805,
        "Mt" : 266,
        "Md" : 258,
        "Hg" : 200.59,
        "Mo" : 95.94,
        "Nd" : 144.24,
        "Ne" : 20.1797,
        "Np" : 237.048,
        "Ni" : 58.6934,
        "Nb" : 92.90638,
        "N"  : 14.00674,
        "No" : 259,
        "Os" : 190.2,
        "O"  : 15.9994,
        "Pd" : 106.42,
        "P"  : 30.973762,
        "Pt" : 195.08,
        "Pu" : 244,
        "Po" : 209,
        "K"  : 39.0983,
        "Pr" : 140.90765,
        "Pm" : 145,
        "Pa" : 231.0359,
        "Ra" : 226.025,
        "Rn" : 222,
        "Re" : 186.207,
        "Rh" : 102.90550,
        "Rb" : 85.4678,
        "Ru" : 101.07,
        "Rf" : 261,
        "Sm" : 150.36,
        "Sc" : 44.955910,
        "Sg" : 263,
        "Se" : 78.96,
        "Si" : 28.0855,
        "Ag" : 107.8682,
        "Na" : 22.989768,
        "Sr" : 87.62,
        "S"  : 32.066,
        "Ta" : 180.9479,
        "Tc" : 98,
        "Te" : 127.60,
        "Tb" : 158.92534,
        "Tl" : 204.3833,
        "Th" : 232.0381,
        "Tm" : 168.93421,
        "Sn" : 118.710,
        "Ti" : 47.88,
        "W"  : 183.85,
        "U"  : 238.0289,
        "V"  : 50.9415,
        "Xe" : 131.29,
        "Yb" : 173.04,
        "Y"  : 88.90585,
        "Zn" : 65.39,
        "Zr" : 91.224,
        "Xx" : 0.0
}

def read_xyz(filename):
  with open(filename, "r") as fin:
    line = fin.readline()
    natoms = int(line)
    xyz = np.zeros((natoms, 3))
    m = np.zeros(natoms)
    line = fin.readline()
    for i in range(natoms):
      line = fin.readline()
      label, x, y, z = line.split()
      label = label.lower()
      label = label.capitalize()
      m[i] = masses[label]
      xyz[i,:] = float(x), float(y), float(z)
  fin.close()
  return xyz, m


def inertia(m, xyz):
  # calculate inertia tensor
  inertia = np.zeros((3, 3))
  for i, v in enumerate(xyz):
    inertia[0, 0] += m[i]*(v[1]**2 + v[2]**2)
    inertia[1, 1] += m[i]*(v[0]**2 + v[2]**2)
    inertia[2, 2] += m[i]*(v[0]**2 + v[1]**2)
    inertia[0, 1] -= m[i]*v[0]*v[1]
    inertia[0, 2] -= m[i]*v[0]*v[2]
    inertia[1, 2] -= m[i]*v[1]*v[2]
  inertia[1,0] = inertia[0,1]
  inertia[2,0] = inertia[0,2]
  inertia[2,1] = inertia[1,2]

  # diagonalize inertia tensor
  w, v = LA.eig(inertia)

  v = v[:, w.argsort()]
  w.sort()
  return w, v


class Cluster:
  def __init__(self, root, file):
    self.root = root
    self.file = file
    self.filename = os.path.join(root, file)
    self.xyz, self.m = read_xyz(self.filename)
    self.mtot = sum(self.m)

    # shift COM to origin
    com = np.zeros(3)
    for i, v in enumerate(self.xyz):
      com += self.m[i]*v

    com /= self.mtot
    self.xyz -= com

    self.w, self.v = inertia(self.m, self.xyz)



filenames = ["input.xyz"]

dflag = False

try:
    opts, args = getopt.getopt(sys.argv[1:], "hd")
except getopt.GetoptError:
    print("Usage:\t{} [options] <filename>\n\t[-d] \t delete duplicates".format(sys.argv[0]))
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-d"):
      dflag = True
    if opt in ("-h"):
      print("Usage:\t{} [options] <filename>\n\t[-d] \t delete duplicates".format(sys.argv[0]))
      sys.exit(2)

if len(args) != 0:
  filenames = []
  for a in args:
    filenames.append(a)

for filename in filenames:

  filelist = []
  clusters = []

  for root, dirs, files in os.walk("."):
    for file in files:
      if file.endswith(filename):
        filelist.append([root, file])

  filelist.sort()

  for f in filelist:
#    print(f)
    clusters.append(Cluster(f[0], f[1]))

  duplicates = []

  for i, c1 in enumerate(clusters):
    for j, c2 in enumerate(clusters):
      if (j > i):
        if (c2.mtot != c1.mtot):
          break
        if (c2.mtot == c1.mtot):
          gdist = pow(pow((1/c1.w[0]-1/c2.w[0])/(1/c1.w[0]), 2)
                + pow((1/c1.w[1]-1/c2.w[1])/(1/c1.w[1]), 2)
                + pow((1/c1.w[2]-1/c2.w[2])/(1/c1.w[2]), 2), 0.5)
          if (gdist < 0.01):
            duplicates.append(c2.root)
#            print(c1.filename, c2.filename, gdist)

  duplicates = list(dict.fromkeys(duplicates))
  duplicates.sort()

  for f in duplicates:
    print(f)

  if (dflag):
    for f in duplicates:
      shutil.rmtree(f)

























