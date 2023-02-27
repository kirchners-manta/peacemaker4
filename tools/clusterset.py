#!/usr/bin/env python3

import sys
import os
import re
import numpy as np
import math
import getopt
import itertools


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


class Cluster:
  def __init__(self, filename):
    self.filename = filename
    self.name = filename[:-4]
    self.compNames = ""
    self.comp = np.zeros(len(base), dtype=int)
    self.volume = 0
    self.energy = 0
    self.ebind = 0
    self.sigma = 1
    self.isMonomer = False

    self.getComposition()
    self.getEnergy()
    self.getSigma()

  def getComposition(self):
    for i, b in enumerate(base):
      if (b in self.name):
        rest = self.name.split(b, 1)[1]
        rest = re.sub(r'[-{}]', ' ', rest)
        rest = rest.split()[0]
        if rest and rest[0].isdigit():
          self.comp[i] = int(rest[0])

  # TODO: Implement function to calculate sigma.
  def getSigma(self):
    self.sigma = 1

  def getVolume(self):
    with open(self.filename, "r") as fin:
      natoms = int(fin.readline())
      comment = fin.readline()
      for i in range(natoms):
        label = fin.readline().split()[0]
        self.volume += 4.0/3.0*math.pi*elements[label][0]**3

  def getEnergy(self):
    with open(self.filename, "r") as f:
      line = f.readline()
      self.energy = float(f.readline())


base = []
nmax = 9999

options = {}

try:
  opts, args = getopt.getopt(sys.argv[1:], "")
except getopt.GetoptError:
  print("Usage:\t{} [options] <basename> <basename2> ... [n_max]".format(sys.argv[0]))
  sys.exit(2)
for opt, arg in opts:
  if opt in options:
    options[opt] = True

if len(args) == 0:
  print("Usage:\t{} [options] <basename> <basename2> ... [n_max]".format(sys.argv[0]))
  sys.exit(1)


for i, arg in enumerate(args):
  if (arg.isdigit()):
    nmax = int(arg)
    continue
  base.append(arg)

clusters = []

for filename in os.listdir("."):
  if filename.endswith(".xyz"):
    if any(filename[0] == b for b in base):
      name = filename[:-4]
      name = name.split("-")
      name = re.split("[0-9]", name[0])
      name = "".join(name)
      if all(n in base for n in name):
        clusters.append(Cluster(filename))
        clusters[-1].compNames=name


clist = [[] for _ in range(len(base)+1)]

for c in clusters:
  for i in range(len(base)):
    if (sum(c.comp) <= nmax):
      if (sum(c.comp) == c.comp[i]):
        clist[i].append(c)
        break
      elif (c.comp[i] != 0):
        clist[-1].append(c)
        break

monomerEnergies = np.zeros(len(base))

for i, l in enumerate(clist):
  l.sort(key=lambda x: (len(x.name), x.compNames, x.name, sum(x.comp)))
  for c in l:
    if (sum(c.comp) == 1):
      c.isMonomer = True
      c.getVolume()
      monomerEnergies[i] = c.energy
      break

for l in clist:
  l.sort(key=lambda x: (len(x.name), x.compNames, x.name, sum(x.comp)))
  for c in l:
    c.ebind = c.energy
    for i in range(len(base)):
      c.ebind = c.ebind - c.comp[i] * monomerEnergies[i]
    c.ebind = c.ebind * 2625.5


lst = list(itertools.product([0, 1], repeat=len(base)))

for l in lst[1:]:
    with open(''.join(base[i]*l[i] for i in range(len(base)))+".clusterset", "w") as f:
      for cl in clist:
        for c in cl:
          if (all(c.comp[k] * l[k] == c.comp[k] for k in range(len(c.comp))) and (sum(c.comp) <= nmax) and (sum(c.comp) != 0)):
            f.write("[" + c.name + "]\n")
            f.write("  composition " + ''.join((str(c.comp[j])+" ")*l[j] for j in range(len(c.comp))) + "\n")
            if (c.isMonomer):
              f.write("  monomer\n")
              f.write("  volume " + str(c.volume) + "\n")
              monomer = False
            f.write("  sigma " + str(c.sigma) + "\n")
            f.write("  energy " + str(c.ebind) + "\n")
            f.write("  coordinates ../clusters/" + c.filename + "\n")
            f.write("  frequencies ../clusters/" + c.name + ".flist\n\n")


    f.close()
