#!/usr/bin/env python3

import subprocess
import sys
import os.path
import getopt
import numpy as np
from numpy import linalg as LA


iflag = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], "")
except getopt.GetoptError:
    print("Usage:\t{} <datafile> <number> <datafile2> <number2> ...".format(sys.argv[0]))
    sys.exit(2)
if len(args) % 2 != 0:
  print("Usage:\t{} <datafile> <number> <datafile2> <number2> ...".format(sys.argv[0]))
  sys.exit(2)

mols = int(len(args)/2)
files = []
numbers = [[] for _ in range(mols)]

for i in range(len(args)):
    if i % 2 == 0:
       files.append(args[i])
    else:
      numbers[int((i-1)/2)] = int(args[i])


lines = [[] for _ in range(mols)]
comments = [[] for _ in range(mols)]

natoms = [[] for _ in range(mols)]
atypes = [[] for _ in range(mols)]
nbonds = [[] for _ in range(mols)]
btypes = [[] for _ in range(mols)]
nangles = [[] for _ in range(mols)]
antypes = [[] for _ in range(mols)]
ndihedrals = [[] for _ in range(mols)]
dtypes = [[] for _ in range(mols)]
nimpropers = [[] for _ in range(mols)]
itypes = [[] for _ in range(mols)]

masses = [[] for _ in range(mols)]
pairc = [[] for _ in range(mols)]
bondc = [[] for _ in range(mols)]
anglec = [[] for _ in range(mols)]
dihedc = [[] for _ in range(mols)]
impropc = [[] for _ in range(mols)]

atoms = [[] for _ in range(mols)]
bonds = [[] for _ in range(mols)]
angles = [[] for _ in range(mols)]
dihedrals = [[] for _ in range(mols)]
impropers = [[] for _ in range(mols)]


for i in range(mols):
  with open(files[i], "r") as fdata:
    # open general LAMMPS data file and save content to memory
    tmplines = fdata.readlines()

    for j in range(len(tmplines)):
      sline, scomment = (tmplines[j].split("#") + ["\n"]*99)[:2]
      lines[i].append(sline)
      comments[i].append(scomment)

    for j in range(len(lines[i])):
      if "atoms" in lines[i][j]:
        natoms[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        natoms[i] = int(tmplist[0])
      elif "atom types" in lines[i][j]:
        atypes[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        atypes[i] = int(tmplist[0])
      elif "bonds" in lines[i][j]:
        nbonds[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        nbonds[i] = int(tmplist[0])
      elif "bond types" in lines[i][j]:
        btypes[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        btypes[i] = int(tmplist[0])
      elif "angles" in lines[i][j]:
        nangles[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        nangles[i] = int(tmplist[0])
      elif "angle types" in lines[i][j]:
        antypes[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        antypes[i] = int(tmplist[0])
      elif "dihedrals" in lines[i][j]:
        ndihedrals[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        ndihedrals[i] = int(tmplist[0])
      elif "dihedral types" in lines[i][j]:
        dtypes[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        dtypes[i] = int(tmplist[0])
      elif "impropers" in lines[i][j]:
        nimpropers[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        nimpropers[i] = int(tmplist[0])
      elif "improper types" in lines[i][j]:
        itypes[i].append([0])
        tmplist = (lines[i][j].split())[:2]
        itypes[i] = int(tmplist[0])



      elif "Masses" in lines[i][j]:
        for k in range(int(atypes[i])):
          masses[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          masses[i][k] = list(map(float, tmplist))
          masses[i][k].append(comments[i][j+k+2].rstrip())

      elif "Pair Coeffs" in lines[i][j]:
        for k in range(int(atypes[i])):
          pairc[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          pairc[i][k] = list(map(float, tmplist))
          pairc[i][k].append(comments[i][j+k+2].rstrip())

      elif "Bond Coeffs" in lines[i][j]:
        for k in range(int(btypes[i])):
          bondc[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          bondc[i][k] = list(map(float, tmplist))
          bondc[i][k].append(comments[i][j+k+2].rstrip())

      elif "Angle Coeffs" in lines[i][j]:
        for k in range(int(antypes[i])):
          anglec[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          anglec[i][k] = list(map(float, tmplist))
          anglec[i][k].append(comments[i][j+k+2].rstrip())

      elif "Dihedral Coeffs" in lines[i][j]:
        for k in range(int(dtypes[i])):
          dihedc[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          dihedc[i][k] = list(map(float, tmplist))
          dihedc[i][k].append(comments[i][j+k+2].rstrip())

      elif "Improper Coeffs" in lines[i][j]:
        for k in range(int(itypes[i])):
          impropc[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          impropc[i][k] = list(map(float, tmplist))
          impropc[i][k].append(comments[i][j+k+2].rstrip())

      elif "Atoms" in lines[i][j]:
        for k in range(int(natoms[i])):
          atoms[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          atoms[i][k] = list(map(float, tmplist))
          atoms[i][k].append(comments[i][j+k+2].rstrip())

      elif "Bonds" in lines[i][j]:
        for k in range(int(nbonds[i])):
          bonds[i].append([0,0,0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          bonds[i][k] = list(map(float, tmplist))

      elif "Angles" in lines[i][j]:
        for k in range(int(nangles[i])):
          angles[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          angles[i][k] = list(map(float, tmplist))

      elif "Dihedrals" in lines[i][j]:
        for k in range(int(ndihedrals[i])):
          dihedrals[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          dihedrals[i][k] = list(map(float, tmplist))

      elif "Impropers" in lines[i][j]:
        for k in range(int(nimpropers[i])):
          impropers[i].append([0])
          tmplist = lines[i][j+k+2].split()
          del tmplist[0]
          impropers[i][k] = list(map(float, tmplist))

  fdata.close


for i in range(mols):

  n_max = max([item[0] for item in atoms[i]])
  for j in range(1, int(numbers[i])):
    for k in range(natoms[i]):
      atoms[i].append(atoms[i][k][:])
      atoms[i][k+j*natoms[i]][0] += n_max*j

  if bonds[i]:
    for j in range(1, int(numbers[i])):
      for k in range(nbonds[i]):
        bonds[i].append(bonds[i][k][:])
        bonds[i][k+j*nbonds[i]][1] += natoms[i]*j
        bonds[i][k+j*nbonds[i]][2] += natoms[i]*j

  if angles[i]:
    for j in range(1, int(numbers[i])):
      for k in range(nangles[i]):
        angles[i].append(angles[i][k][:])
        angles[i][k+j*nangles[i]][1] += natoms[i]*j
        angles[i][k+j*nangles[i]][2] += natoms[i]*j
        angles[i][k+j*nangles[i]][3] += natoms[i]*j

  if dihedrals[i]:
    for j in range(1, int(numbers[i])):
      for k in range(ndihedrals[i]):
        dihedrals[i].append(dihedrals[i][k][:])
        dihedrals[i][k+j*ndihedrals[i]][1] += natoms[i]*j
        dihedrals[i][k+j*ndihedrals[i]][2] += natoms[i]*j
        dihedrals[i][k+j*ndihedrals[i]][3] += natoms[i]*j
        dihedrals[i][k+j*ndihedrals[i]][4] += natoms[i]*j

  if impropers[i]:
    for j in range(1, int(numbers[i])):
      for k in range(nimpropers[i]):
        impropers[i].append(impropers[i][k][:])
        impropers[i][k+j*nimpropers[i]][1] += natoms[i]*j
        impropers[i][k+j*nimpropers[i]][2] += natoms[i]*j
        impropers[i][k+j*nimpropers[i]][3] += natoms[i]*j
        impropers[i][k+j*nimpropers[i]][4] += natoms[i]*j

  natoms[i] = natoms[i] * numbers[i]
  nbonds[i] = nbonds[i] * numbers[i]
  nangles[i] = nangles[i] * numbers[i]
  ndihedrals[i] = ndihedrals[i] * numbers[i]
  nimpropers[i] = nimpropers[i] * numbers[i]

atypestotal = 0
btypestotal = 0
antypestotal = 0
dtypestotal = 0
itypestotal = 0
natomstotal = 0

for i in range(1, mols):

  natomstotal += natoms[i-1]
  atypestotal += atypes[i-1]
  if bonds[i-1]:
    btypestotal += btypes[i-1]
  if angles[i-1]:
    antypestotal += antypes[i-1]
  if dihedrals[i-1]:
    dtypestotal += dtypes[i-1]
  if impropers[i-1]:
    itypestotal += itypes[i-1]
  n_max = max([item[0] for item in atoms[i-1]])

  for k in range(natoms[i]):
      atoms[i][k][0] += n_max
      atoms[i][k][1] += atypestotal


  if bonds[i]:
    for k in range(nbonds[i]):
      bonds[i][k][0] += btypestotal
      bonds[i][k][1] += natomstotal
      bonds[i][k][2] += natomstotal

  if angles[i]:
    for k in range(nangles[i]):
      angles[i][k][0] += antypestotal
      angles[i][k][1] += natomstotal
      angles[i][k][2] += natomstotal
      angles[i][k][3] += natomstotal

  if dihedrals[i]:
    for k in range(ndihedrals[i]):
      dihedrals[i][k][0] += dtypestotal
      dihedrals[i][k][1] += natomstotal
      dihedrals[i][k][2] += natomstotal
      dihedrals[i][k][3] += natomstotal
      dihedrals[i][k][4] += natomstotal

  if impropers[i]:
    for k in range(nimpropers[i]):
      impropers[i][k][0] += itypestotal
      impropers[i][k][1] += natomstotal
      impropers[i][k][2] += natomstotal
      impropers[i][k][3] += natomstotal
      impropers[i][k][4] += natomstotal


print("LAMMPS data file for use with OGOLEM\n")
print("{} atoms".format(sum(natoms)))
print("{} atom types".format(sum(atypes)))
for i in range(mols):
  if bonds[i]:
    print("{} bonds".format(sum(filter(None, nbonds))))
    print("{} bond types".format(sum(filter(None, btypes))))
    break
for i in range(mols):
  if angles[i]:
    print("{} angles".format(sum(filter(None, nangles))))
    print("{} angle types".format(sum(filter(None, antypes))))
    break
for i in range(mols):
  if dihedrals[i]:
    print("{} dihedrals".format(sum(filter(None, ndihedrals))))
    print("{} dihedral types".format(sum(filter(None, dtypes))))
    break
for i in range(mols):
  if impropers[i]:
    print("{} impropers".format(sum(filter(None, nimpropers))))
    print("{} improper types".format(sum(filter(None, itypes))))
    break

print("\n-100.0 100.0 xlo xhi")
print("-100.0 100.0 ylo yhi")
print("-100.0 100.0 zlo zhi")


print("\nMasses\n")
k=0
for i in range(mols):
  for j in range(atypes[i]):
    k += 1
    print("%4d   %16.4f   # %s" % (k, masses[i][j][0], masses[i][j][1]) )

print("\nPair Coeffs\n")
k=0
for i in range(mols):
  if pairc[i]:
    for j in range(atypes[i]):
      k += 1
      print("%4d   %16.6f   %16.6f   # %s" % (k, pairc[i][j][0], pairc[i][j][1], pairc[i][j][2]) )

print("\nAtoms\n")
k=0
for i in range(mols):
  if atoms[i]:
    for j in range(natoms[i]):
      k += 1
      print("%4d %4d %4d  %16.6f   %16.6f  %16.6f   %16.6f  # %s" % (k, atoms[i][j][0], atoms[i][j][1], atoms[i][j][2], atoms[i][j][3], atoms[i][j][4], atoms[i][j][5], atoms[i][j][6]) )


for l in range(mols):
  if bonds[l]:
    print("\nBond Coeffs\n")
    k=0
    for i in range(mols):
      if bondc[i]:
        for j in range(btypes[i]):
          k += 1
          print("%4d   %16.6f   %16.6f   # %s" % (k, bondc[i][j][0], bondc[i][j][1], bondc[i][j][2]) )

    print("\nBonds\n")
    k=0
    for i in range(mols):
      if bonds[i]:
        for j in range(nbonds[i]):
          k += 1
          print("%4d %4d %4d %4d" % (k, bonds[i][j][0], bonds[i][j][1], bonds[i][j][2]) )
    break


for l in range(mols):
  if angles[l]:
    print("\nAngle Coeffs\n")
    k=0
    for i in range(mols):
      if anglec[i]:
        for j in range(antypes[i]):
          k += 1
          print("%4d   %16.6f   %16.6f   # %s" % (k, anglec[i][j][0], anglec[i][j][1], anglec[i][j][2]) )

    print("\nAngles\n")
    k=0
    for i in range(mols):
      if angles[i]:
        for j in range(nangles[i]):
          k += 1
          print("%4d %4d %4d %4d %4d" % (k, angles[i][j][0], angles[i][j][1], angles[i][j][2], angles[i][j][3]) )
    break


for l in range(mols):
  if dihedrals[l]:
#    print("\nDihedral Coeffs\n")
#    k=0
#    for i in range(mols):
#      if dihedc[i]:
#        for j in range(dtypes[i]):
#          k += 1
#          print("%4d   %16.6f   %16i   %16i   # %s" % (k, dihedc[i][j][0], dihedc[i][j][1], dihedc[i][j][2], dihedc[i][j][3]) )
    print("\nDihedral Coeffs\n")
    k=0
    for i in range(mols):
      if dihedc[i]:
        for j in range(dtypes[i]):
          k += 1
          print("%4d   %16.6f   %16.6f   %16.6f   %16.6f   # %s" % (k, dihedc[i][j][0], dihedc[i][j][1], dihedc[i][j][2], dihedc[i][j][3], dihedc[i][j][4]) )

    print("\nDihedrals\n")
    k=0
    for i in range(mols):
      if dihedrals[i]:
        for j in range(ndihedrals[i]):
          k += 1
          print("%4d %4d %4d %4d %4d %4d" % (k, dihedrals[i][j][0], dihedrals[i][j][1], dihedrals[i][j][2], dihedrals[i][j][3], dihedrals[i][j][4]) )
    break


for l in range(mols):
  if impropers[l]:
    print("\nImproper Coeffs\n")
    k=0
    for i in range(mols):
      if impropc[i]:
        for j in range(itypes[i]):
          k += 1
          print("%4d   %16.6f   %16.6f   # %s" % (k, impropc[i][j][0], impropc[i][j][1], impropc[i][j][2]) )

    print("\nImpropers\n")
    k=0
    for i in range(mols):
      if impropers[i]:
        for j in range(nimpropers[i]):
          k += 1
          print("%4d %4d %4d %4d %4d %4d" % (k, impropers[i][j][0], impropers[i][j][1], impropers[i][j][2], impropers[i][j][3], impropers[i][j][4]) )
    break


print("\n")















































