Compiling and Running Peacemaker 4
#################################

Compiling Peacemaker 4
-----------------------------
Peacemaker is a modern FORTRAN code and thus requires a modern FORTRAN compiler.
We recommend a recent version of gfortran which is used for active development. 
Before compiling Peacemaker, make sure that the following dependencies are installed:

- `Meson <https://mesonbuild.com/>`_ 
- `Ninja <https://ninja-build.org/>`_

Meson and Ninja are build systems that help to configure and compile the code efficiently.
To compile Peacemaker, follow these steps:
1. Clone the Peacemaker repository from GitHub:
   ```
   git clone ``git@github.com:kirchners-manta/peacemaker4.git``
   ```
2. Navigate to the cloned directory:
   ```
   cd peacemaker4/pm4
    ``` 
3. Create a build directory:
   ```
   meson setup build
   ```
4. Build the code using Ninja:
   ```
   ninja -C build
   ```
5. Optionally, you can run tests to ensure everything is working correctly:
   ```
   ninja -C build test
   build/test/tester
   ```