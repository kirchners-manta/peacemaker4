project('peacemaker', 'fortran',
  version : '0.1',
  default_options : ['warning_level=3'])

sources = files(
    'src/kinds.f90',
    'src/lengths.f90',
    'src/constants.f90',
    'src/iso_varying_string.f90',
    'src/auxiliary.f90',
    'src/error.f90',
    'src/atomic_data.f90',
    'src/config.f90',
    'src/input.f90',
    'src/cluster.f90',
    'src/info.f90',
    'src/polynomial.f90',
    'src/shared_data.f90',
    'src/partition_functions.f90',
    'src/thermo.f90',
    'src/qce.f90',
    'src/reader.f90',
    'main.f90'
)

# Create TOML Fortran as subproject
tomlf_dep = dependency(
  'toml-f',
  version: '>=0.2.0',
  fallback: ['toml-f', 'tomlf_dep'],
  default_options: ['default_library=static'],
)

exe = executable('peacemaker', sources, 
                 fortran_args : ['-cpp', '-fopenmp'],
                 link_args: ['-fopenmp'],
                 dependencies : [tomlf_dep],
                install : true)

test('basic', exe)
