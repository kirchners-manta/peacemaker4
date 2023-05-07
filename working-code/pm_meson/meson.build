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
    'main.f90'
)

exe = executable('peacemaker', sources, 
                 fortran_args : ['-cpp'],
                 link_args: ['-fopenmp'],
                install : true)

test('basic', exe)
