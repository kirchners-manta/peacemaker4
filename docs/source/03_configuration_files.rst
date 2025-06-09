The toml format
============

The toml format is a human-readable configuration file format that is easy to read and write.
It is used to define the parameters for the QCE input file and the clusterset file in Peacemaker.
More information about the toml format can be found in the `TOML documentation <https://toml.io/en/>`_.

A typical toml file is structured with sections, each defined by square brackets, and key-value pairs 
within those sections.
The key-value pairs are separated by an equal sign, and values can be strings, numbers, arrays or boleans.
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
* Booleans are written as `true` or `false`.

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
.. line-block::
    **components = N**
    The number of components in the system. :math:`N = 1` for a pure system, :math:`N = 1` for a binary mixture, :math:`N = 1` for a ternary mixture, etc. Note that although it is possible to run a pure system as binary system, where the amount of one of the species is set to zero, we strongly encourage you to run such calculations as a pure system. Results will be the same in either case, but slow convergence may arise for some temperatures if the amount of monomers of one component is sufficiently small.
    *Optional. Default: 1*

[qce]
------------------------------
.. line-block::
    **amf = [A]** 
    **amf = [min, max, steps]**
       The mean field parameter :math:`a_{mf}` in units of :math:`\mathrm{J m^3 mol^{-2}}`. Can be specified either as a single value **A**, or as a range of three values, where **min** is the start, **max** the end, and **steps** the number of data points (including both boundaries).
       *Optional. Default: 0.0*

    **bxv = [A]**
    **bxv = [min, max, steps]**
       The exclusion volume scaling parameter :math:`b_{xv}`. Can be specified either as a single value **A**, or as a range of three values, where **min** is the start, **max** the end, and **steps** the number of data points (including both boundaries).
       *Optional. Default: 1.0*

    **amf_temp = [A]**
    **amf_temp = [min, max, steps]**
       The linear temperature dependence parameter :math:`a_{mf,temp}` of the mean field. The specification is similar to the one for :math:`a_{mf}`. This is an experimental feature and should only be used with care.
       *Optional. Default: 0.0*

    **bxv_temp = [A]**
    **bxv_temp = [min, max, steps]**
       The linear temperature dependence parameter :math:`b_{xv,temp}` of the exclusion volume. The specification is similar to the one for :math:`b_{xv}`. This is an experimental feature and should only be used with care.
       *Optional. Default: 0.0*

    **grid_iterations = N**
       The number of iterations for the parameter sampling if a sampling grid is specified. With each iteration, the grid center is moved to the best parameter pair and the grid size is decreased with a factor of 0.2.
       *Optional. Default: 1*

    **rotor_cutoff = A**
       The cutoff frequency in :math:`cm^{-1}` at which the RRHO-correction for low frequencies will be used. To limit their influence on the entropy, vibrational modes with a frequency below A will be treated as hindered rotations, employing a switching function to smooth the transition between harmonic oscillator and rigid rotator. If set to 0, no correction will be applied.
       *Optional. Default: 0*

    **optimizer = ["amf",  "bxv",  "amf_temp",  "bxv_temp"]**
       Enables the Nelder-Mead algorithm for parameter optimization. Possible values are:
       **"amf"**: Optimize the mean field parameter :math:`a_{mf}`
       **"bxv"**: Optimize the exclusion volume scaling parameter :math:`b_{xv}`
       **"amf_temp"**: Optimize the linear temperature dependence of the mean field parameter :math:`a_{mf,temp}`
       **"bxv_temp"**: Optimize the linear temperature dependence of the exclusion volume parameter :math:`b_{xv,temp}`
       Parameters can be given in any combination and order. By default, no optimization is performed.

    **max_deviation = A**
       The maximum relative deviation of the Gibbs energy. Used to check convergence of the QCE iteration. A QCE cycle has converged, if 

.. math::
 
    |\frac{G_{i} - G_{i-1}}{G_{i-1}}| < A .

.. line-block::
       where :math:`G_i` is the Gibbs energy of the i-th iteration.
       *Optional. Default: 1.0e-9*

    **volume_damping_factor = A**
       The volume damping factor used to damp the initial volume guess if one of the polynomials did not converge. Shall be between 0 and 1. Damping is performed by :math:`\gamma_V = 1 \pm A`, depending on the mode of the temperature loop.
       *Optional. Default: 0.01*

    **qce_iterations = N**
        The maximum number of iterations in a QCE cycle.
        *Optional. Default: 100*

    **newton_iterations = N**
        The maximum number of iterations in the Newton-Raphson cycle used to solve the n d-dimensional population polynomial equations.
        *Optional. Default: 100*

[ensemble]
------------------------------
.. line-block::
    **temperature = [A]**
    **temperature = [min, max, steps]**
        The temperature in units of :math:`K`. Can be specified either as a single value **A**, or as a range of three values, where **min** is the start, **max** the end, and **steps** the number of data points (including both boundaries).
        *Optional. Default: 298.15*

    **pressure = A**
        The pressure in units of :math:`bar`. 
        *Optional. Default: 1.01325.0*

    **monomer_amounts = [N, M, ...]**
        The molar amounts of the components in the system. The number of values must match the number of components specified in the **system** section. The values are given in units of :math:`mol` and must sum up to 1.0.
        *Required.*

[reference]
------------------------------
This section is optioanl. It enables comparison to experimental reference data.
It is disabled by default.
Further details on parameter sampling are given in the last section of the documentation.

.. line-block::
    **density = [A, B]**
    **density = [A, B, C]**
        Reference density **B** in units of :math:`g cm^{-3}` at reference temperature **A** in :math:`K` and an optional error weight **C**.
        *Optional.*

    **isobar_file = "path/to/isobar/file"**
    **isobar_weight = A**
        Path to an isobar file and an optional error weight **A**. Isobar files contain two columns representing the temperature in :math:`K` and volume in :math:`L`.
        *Optional.*

    **phase_transition = [A]**
    **phase_transition = [A, B]**
        Reference temperature of phase transition **A** in units of K and an optional error weight **B**. 
        *Optional.*

[output]
------------------------------
This section is optional. It enables the output of additional files and is disabled by default.

.. line-block::
    **contribuions = true/false**
        Enables the output of contributions of each degree of freedom to the thermodynamic quantities. If set to true, contribution output is enabled for all possible thermodynamic quantities, which are helmholtz-contributions, internal-contributions, entropy-contributions and cv-contributions.
        *Optional.*

    **helmholtz_contributions = true/false**
        Enables the output of contributions to the Helmholtz free energy. If set to true, contributions are written to a file called **helmholtz_contrib.dat**.
        *Optional.*

    **internal_contributions = true/false**
        Enables the output of contributions to the internal energy. If set to true, contributions are written to a file called **internal_contrib.dat**.
        *Optional.*

    **entropy_contributions = true/false**
        Enables the output of contributions to the entropy. If set to true, contributions are written to a file called **entropy_contrib.dat**.
        *Optional.*

    **cv_contributions = true/false**
        Enables the output of contributions to the heat capacity at constant volume. If set to true, contributions are written to a file called **cv_contrib.dat**.
        *Optional.*

    **progress_bar = true/false**
        Enables the output of a progress bar during the calculation. If set to true, a progress bar is displayed in the terminal.
        *Optional. Default: true*

Example QCE input files
------------------------------

**Single Point Calculation**
The input file shown bellow will run a QCE "single point" calculation for a one-component system using the clusterset 
specified in the command line and explained in the following section. Default options are used in most cases.

.. code-block:: toml

    [qce]
        amf = 0.1
        bxv = 1.3

    [ensemble]
        temperature = [200.0, 400.0, 201]
        pressure = 1.01325
        monomer_amounts = 1.0

**Parameter Sampling**
This input will perform an :math:`a_{mf}`, :math:`b_{xv}` parameter sampling for a pure substance.
Reference data are provided by an isobar file.

.. code-block:: toml

    [system]
        components = 1

    [qce]
        amf = [0.0, 0.5, 101]
        bxv = [1.0, 2.0, 101]

    [ensemble]
        temperature = [200.0, 400.0, 201]
        pressure = 1.01325
        monomer_amounts = 1.0

    [reference]
        isobar_file = "isobar.dat"

**Parameter Optimization**
The following input will perform an :math:`a_{mf}`, :math:`b_{xv}` parameter optimization for a ternary mixture, following a 
rough sampling on a small grid.
Reference data are provided by a density at :math:`298.15 K` and a temperature of phase transition.

.. code-block:: toml

    [system]
        components = 3

    [qce]
        amf = [0.0, 2.0, 11]
        bxv = [0.5, 1.5, 11]
        optimizer = ["amf", "bxv"]
        grid_iterations = 2

    [ensemble]
        temperature = [273.15, 400.15, 128]
        pressure = 1.01325
        monomer_amounts = [0.6, 0.1, 0.3]

    [reference]
        density = [298.15, 0.9248]
        phase_transition = 332.61


Clusterset file
========================
The clusterset file is a toml file that contains the paths to the structure and frequency files 
of each cluster, as well as information about the clusters, such as their composition or energy.

Here, the sections are the clusters, which are defined by square brackets, and the key-value pairs
within those sections are equivalent for each cluster.
The clusterset file is structured as follows:

.. code-block:: toml

    [cluster1]
        isMonomer = true/false
        composition = [N, M, ...]
        sigma = N
        coordinates = "path/to/structure/file1.xyz"
        frequencies = "path/to/frequency/file1.dat"
        energy = A
        volume = A
        frequency_scale = A
        anharmonicity = A

    [cluster2]
        isMonomer = true/false
        composition = [N, M, ...]
        sigma = N
        coordinates = "path/to/structure/file2.xyz"
        frequencies = "path/to/frequency/file2.dat"
        energy = A
        volume = A
        frequency_scale = A
        anharmonicity = A


The Keywords are explained in detail below:

.. line-block::
    **isMonomer = true/false**
        Specifies whether the cluster is a monomer or not. If set to true, the cluster is treated as a monomer.
        *Optional but must be present once for each component. Default: false*

    **composition = [N, M, ...]**
        The composition of the cluster as an array of integers, where each integer represents the number of molecules of a certain type in the cluster.
        *Required.*

    **sigma = N**
        The rotational symmetry number of the cluster.
        *Optional. Default: 1*

    **coordinates = "path/to/structure/file1.xyz"**
        The path to the structure file of the cluster in XYZ format. Units are Angstrom.
        *Required.*

    **frequencies = "path/to/frequency/file1.dat"**
        Path to a frequency file. It contains the number of frequencies in line 1, followed by a comment line, followed by one frequency per line. Units are :math:`cm^{-1}`.
        *Required.*

    **energy = A**
        The adiabatic interaction energy of the cluster in units of :math:`\mathrm{kJ mol^{-1}}` (negative energies represent stable clusters).
        *Required.*

    **volume = A**
        The volume of the cluster in units of :math:`\mathrm{A^3}`.
        Must only be specified for monomers.

    **frequency_scale = a**
        The frequency scaling factor for the cluster. This is used to scale the frequencies of the cluster.
        *Optional. Default: 1.0*

    **anharmonicity = A**
        Anharmonicity constant for the cluster.
        *Optional. Default: 0.0*
