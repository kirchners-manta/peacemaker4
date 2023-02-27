#!/usr/bin/python3.8

import sys
import os, errno
import subprocess
import shutil
import fileinput
import getopt
import numpy as np
from scipy.optimize import differential_evolution

# Usage: pymaker2.py qce.input w.clusterset
#
# Run pymaker2.py in a directory with the necessary input files prepared.
# Use a QCE input file similar to the example shown below and add additional
# settings if needed. Set the parameters to single values in the QCE input
# file. Do not perform a parameter sampling, its results will be ignored in
# the interactive mode. Use the interactive key word and set the optimizer
# option to 1111, which allows the optimization of all four parameters amf,
# bxv, amf_temp, and bxv_temp. Set the boundaries and starting values at the
# bottom of this script. If any parameter should not be optimized, set both
# its boundaries to the same value, e.g. (0.0, 0.0). It is recommended to
# not optimize amf_temp and constrain it to 0.0.
# Do not optimize amf_temp or bxv_temp if you don't provide an isobar as
# experimental reference.
#
# The qce.input file should look like this:
#
# [system]
#     components 1
#
# [qce]
#     amf 0.0
#     bxv 1.0
#     amf_temp 0.0
#     bxv_temp 0.0
#     interactive
#     optimizer 1111
#
# [ensemble]
#     temperature 273.15 399.15 127 # K
#     pressure 1.01325 # bar
#     monomer_amounts 1.0 # mol
#
# [reference]
#     isobar isob.dat
#
# [output]
#     noprogress
#


counter = 0

# Defines the function to be optimized. Takes a set of QCE parameters and runs a
# single point calculation, then returns the error as function value.
def fx(x):

    # Set default error/function value very high, in case that Peacemaker crashes
    # and doesn't return an error.
    fx=1.0E10

    global counter

    # Write the file out again
    with open("tmp.input", 'w') as file:
      file.write("{:.12e} {:.12e} {:.12e} {:.12e}".format(x[0], x[1], x[2], x[3]))
      file.close()
    shutil.move("tmp.input", "imode.input")

    while (not os.path.isfile("sp.out")):
        pass

    error = ''
    while (not error):
      with open("sp.out", "r") as file:
        error = file.readlines()
        file.close()

    shutil.move("sp.out", "old.out")

    fx = float(error[0])

    counter = counter + 1
    print("%i. amf = %.6E, bxv = %.6E, amf_temp = %.6E, bxv_temp = %.6E, error = %.6E" % (counter, x[0], x[1], x[2], x[3], fx))

    return fx



try:
    opts, args = getopt.getopt(sys.argv[1:], "")
except getopt.GetoptError:
    print("Usage:\t{} <qce input> <clusterset>".format(sys.argv[0]))
    sys.exit(2)
for opt, arg in opts:
    pass

if len(args) != 2:
    print("Usage:\t{} <qce input> <clusterset>".format(sys.argv[0]))
    sys.exit(1)

qce_input = args[0]
clusterset = args[1]



# Define start values and boundaries of amf, bxv, amf_temp, and bxv_temp in that order.
x0=[1.0, 1.0, 0.0, 0.005]
bounds=[(0.0, 2.0), (0.5, 1.5), (0.0, 0.0), (0.0, 0.01)]

# Start optimization algorithm to optimize function fx, which
# returns an error value for a given set of QCE parameters.
with open("peace.out", "w") as outfile:
    subprocess.Popen(["/home/ingenmey/peacemaker2-svn/branches/temperature/peacemaker", qce_input, clusterset], stdout=outfile)

    ret = differential_evolution(fx, bounds=bounds, strategy='best1bin')

    # When done, print result.
    print("global minimum: amf = %.10E, bxv = %.10E, amf_temp = %.10E, bxv_temp = %.10E, error = %.6E" % (ret.x[0], ret.x[1], ret.x[2], ret.x[3], ret.fun))

    with open("stop_imode", 'w') as file:
        file.write("global minimum: amf = %.10E, bxv = %.10E, amf_temp = %.10E, bxv_temp = %.10E, error = %.6E" % (ret.x[0], ret.x[1], ret.x[2], ret.x[3], ret.fun))
        file.close()

    outfile.close()
