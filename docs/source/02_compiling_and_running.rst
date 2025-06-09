Compiling Peacemaker 
-----------------------------
Peacemaker is a modern FORTRAN code and thus requires a modern FORTRAN compiler.
We recommend a recent version of gfortran which is used for active development. 
Before compiling Peacemaker, make sure that the following dependencies are installed:

- `Meson <https://mesonbuild.com/>`_ 
- `Ninja <https://ninja-build.org/>`_

Meson and Ninja are build systems that help to configure and compile the code efficiently.
To compile Peacemaker, follow these steps:

1. Clone the `Peacemaker <https://github.com/kirchners-manta/peacemaker4>`_ repository from GitHub:

   .. code-block:: bash

      git clone git@github.com:kirchners-manta/peacemaker4.git

2. Navigate to the cloned directory:

   .. code-block:: bash

      cd peacemaker4/pm4

3. Create a build directory using Meson:

   .. code-block:: bash

      meson setup build

4. Build the code using Ninja:

   .. code-block:: bash

      ninja -C build

5. Optionally, run tests to ensure everything is working correctly:

   .. code-block:: bash

      ninja -C build test
      build/test/tester

A run time optimized binary called **peacemaker** is created in the build directory.
In case of errors, you can adjust the ``meson.build`` to suit your system configuration.
We recommend the following compiler flags for optimal performance:

* ``-O3`` highest optimization level that guarantees standard compliance
* ``-fopenmp`` OpenMP parallelization
* ``-flto`` link-time optimization

.. note::
   Older versions of gfortran are subject to a bug which prevents OpenMP parallelization.
   If you receive the error message ``Attempting to allocate already allocated variable **ib** ``, 
   compile without OpenMP support, or upgrade to a newer compiler version.


Running Peacemaker 
-----------------------------
Two input files are required to run Peacemaker.
The first is the **QCE-input file**, which contains all necessary information about the system 
to be investigated as well as the parameters to be sampled.
The second is the **clusterset file**, which contains the paths to the structure and the frequency 
files of each cluster, as well as information about the clusters, such as their composition or energy.
More details on the structure of these files are given in the following sections.

Peacemaker can be run from the command line as follows:

   .. code-block:: bash

      peacemaker <QCE-input file> <clusterset file> 

If Peacemaker was compiled with OpenMP support, it can be run in parallel by specifying the number of 
threads:

   .. code-block:: bash

      OMP_NUM_THREADS=<number of threads> peacemaker <QCE-input file> <clusterset file> 
