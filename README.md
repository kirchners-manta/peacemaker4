# Peacemaker 4 
**The Quantum Cluster Equilibrium Approach to Liquid Phase Properties**

QCE theory applies statistical mechanics to quantum-chemically optimized clusters to obtain the partition function of the system and any quantity that can be derived therefrom. 
Peacemaker works with pure substances and multicomponent mixtures.

## Where to find the Documentation
You can either read the [online documentation](https://kirchners-mantapeacemaker4.readthedocs.io/en/latest/index.html) or download the [manual](manual_pm4/manual.pdf) in PDF format.

## Compiling Peacemaker
Before compiling Peacemaker, make sure that the following dependencies are installed:

* Meson: [Meson Build system](https://mesonbuild.com/)
* Ninja: [Ninja Build system](https://ninja-build.org/)

To compile Peacemaker, follow these steps:
1. Clone the repository:</p>
 ```git clone git@github.com:kirchners-manta/peacemaker4.git```

2. Change to the directory:</p>
 ```cd peacemaker4/pm4```

3. Create a build directory:</p>
 ```meson setup build```

4. Compile the project:</p>
 ```ninja -C build```

1. Run the automated tests and the unit tests:</p>
 ```ninja -C build test```</br>
 ```build/test/tester```

## Running Peacemaker
Peacemaker is run by

```$ peacemaker [input] [clusterset]```

where `[input]` is the location of the input file and `[clusterset]` is the location of the clusterset file. The structure of both files is explained in Section 4 of the [manual](manual_pm4/manual.pdf) or in the [online documentation](https://kirchners-mantapeacemaker4.readthedocs.io/en/latest/index.html).

## Conversion of input files from Peacemaker 3 to Peacemaker 4
Use </br>
```tools/convert2toml/clusterset2toml.py [clusterset] [clusterset.toml]```</br>
and </br>
```tools/convert2toml/QCEinput2toml.py [input] [input.toml]```</br>
to convert the clusterset and input files from Peacemaker 3 to the new TOML format.
