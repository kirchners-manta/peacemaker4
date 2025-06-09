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
- Strings are enclosed in double quotes, e.g., `"example"`.
- Numbers are written as is, e.g., `3.14`.
- Arrays are enclosed in square brackets, e.g., `[1, 2, 3]`.
