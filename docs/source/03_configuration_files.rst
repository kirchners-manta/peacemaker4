The toml format
============

The toml format is a human-readable configuration file format that is easy to read and write.
It is used to define the parameters for the QCE input file and the clusterset file in Peacemaker.
More information about the toml format can be found in the `TOML documentation <https://toml.io/en/>`_.

A typical toml file is structured with sections, each defined by square brackets, and key-value pairs 
within those sections.
The key-value pairs are separated by an equal sign, and values can be strings, numbers or arrays.
An example of a toml file, which is actually a QCE input file, is shown below:

.. code-block:: toml

   [system]
        components = 1

    [qce]
        amf = [0.0, 0.5, 101]
        bxv = [1.0, 2.0, 101]

   [ensemble]
        temperature = 298.15
        pressure = 101325.0
        monomer_amounts = 1.0

Since the values can have different types, they are presented in different ways:

* Strings are enclosed in double quotes, e.g., `"example"`.
* Numbers are written as is, e.g., `3.14`.
* Arrays are enclosed in square brackets, e.g., `[1, 2, 3]`.

Input files from earlier versions of Peacemaker
------------------------------

If you want to use QCE-input files and clusterset files from earlier versions of Peacemaker,
you can convert them to the toml format using the the following python scripts, located in the
``tools/convert2toml`` directory of the peacemaker4 repository:

* QCEinput2toml.py
* clusterset2toml.py

These scripts can be used from the command line as follows:

.. code-block:: bash

   python QCEinput2toml.py <QCE-input file> <output toml file>
   python clusterset2toml.py <clusterset file> <output toml file>


QCE input file
================
The QCE input file is a toml file that contains all necessary information about the system to be 
investigated and the parameters to be sampled.
In the following, all sections and key-value pairs of the QCE input file are explained in detail.
Sections are defined by square brackets, and key-value pairs are separated by an equal sign as 
shown in the example above.

[system]
------------------------------
**components = N**

    The number of components in the system.
    :math:`N = 1` for a pure system, :math:`N = 1` for a binary mixture, :math:`N = 1` for a ternary mixture, etc.
    Note that although it is possible to run a pure system as binary system, where the amount of 
    one of the species is set to zero, we strongly encourage you to run such calculations as a pure system.
    Results will be the same in either case, but slow convergence may arise for some temperatures 
    if the amount of monomers of one component is sufficiently small.
    Optional. Default: 1

[qce]
------------------------------
**amf = [A]** <br>
**amf = [min, max, steps]**

    The mean field parameter ``amf`` in units of :math:`\mathrm{J\cdot m^3 \cdot mol^{-2}}`.
    Can be specified either as a single value A, or as a range min, max, steps, where min is the start, 
    max the end, and steps the number of data points (including both boundaries).
    Optional. Default: 0.0