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
| **amf = [A]** 
| **amf = [min, max, steps]**

    The mean field parameter :math:`a_{mf}` in units of :math:`\mathrm{J m^3 mol^{-2}}`.
    Can be specified either as a single value **A**, or as a range of three values, where **min** is the start, 
    **max** the end, and **steps** the number of data points (including both boundaries).
    *Optional. Default: 0.0*

| **bxv = [A]**
| **bxv = [min, max, steps]**

    The exclusion volume scaling parameter :math:`b_{xv}`.
    Can be specified either as a single value **A**, or as a range of three values, where **min** is the start,
    **max** the end, and **steps** the number of data points (including both boundaries).
    *Optional. Default: 1.0*

| **amf_temp = [A]**
| **amf_temp = [min, max, steps]**

    The linear temperature dependence parameter :math:`a_{mf,temp}` of the mean field.
    The specification is similar to the one for :math:`a_{mf}`.
    This is an experimental feature and should only be used with care.
    *Optional. Default: 0.0*

| **bxv_temp = [A]**
| **bxv_temp = [min, max, steps]**

    The linear temperature dependence parameter :math:`b_{xv,temp}` of the exclusion volume.
    The specification is similar to the one for :math:`b_{xv}`.
    This is an experimental feature and should only be used with care.
    *Optional. Default: 0.0*

| **grid_iterations = N**

    The number of iterations for the parameter sampling if a sampling grid is specified.
    With each iteration, the grid center is moved to the best parameter pair and the grid size is decreased 
    with a factor of 0.2.
    *Optional. Default: 1*

| **rotor_cutoff = A**

    The cutoff frequency in :math:`cm^{-1}` at which the RRHO-correction for low frequencies will be used.
    To limit their influence on the entropy, vibrational modes with a frequency below A will be treated as 
    hindered rotations, employing a switching function to smooth the transition between harmonic oscillator 
    and rigid rotator. If set to 0, no correction will be applied.
    *Optional. Default: 0*

| **optimizer = ["amf", ~"bxv", ~"amf_temp", ~"bxv_temp"]**

    Enables the Nelder-Mead algorithm for parameter optimization.
    Possible values are:

    * "amf": Optimize only the mean field parameter :math:`a_{mf}`.
    * "bxv": Optimize only the exclusion volume scaling parameter :math:`b_{xv}`.
    * "amf_temp": Optimize only the linear temperature dependence of the mean field parameter :math:`a_{mf,temp}`.
    * "bxv_temp": Optimize only the linear temperature dependence of the exclusion volume :math:`b_{xv,temp}`.
    
    Parameters can be given in any combination and order.
    By default, no optimization is performed.

| **max_deviation = A**

    The maximum relative deviation of the Gibbs energy.
    Used to check convergence of the QCE iteration.
    A QCE cycle has converged, if 

    .. math::
        |\frac{G_{i} - G_{i-1}}{G_{i-1}}| < A .

    where :math:`G_i` is the Gibbs energy of the i-th iteration.
    *Optional. Default: 1.0e-9*

| **volume_damping_factor = A**

    The volume damping factor used to damp the initial volume guess if one of the polynomials did not converge.
    Shall be between 0 and 1.
    Damping is performed by :math:`\gamma_V = 1 \pm A`, depending on the mode of the temperature loop.
    *Optional. Default: 0.01*


.. line-block::
    **qce_iterations = N**
    The maximum number of iterations in a QCE cycle.
    *Optional. Default: 100*

    **newton_iterations = N**
    The maximum number of iterations in the Newton-Raphson cycle used to solve the n d-dimensional population
    polynomial equations.
    *Optional. Default: 100*
