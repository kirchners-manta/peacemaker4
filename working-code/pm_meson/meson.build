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

# Create TOML Fortran as subproject
tomlf_dep = dependency(
  'toml-f',
  version: '>=0.2.0',
  fallback: ['toml-f', 'tomlf_dep'],
  default_options: ['default_library=static'],
)

fortran_args = ['-cpp', '-fopenmp', '-flto', '-O3']

link_args = ['-fopenmp', '-flto', '-O3']

dependencies = [tomlf_dep]

lib = library(
  'peacemaker',
  sources,
  version: meson.project_version(),
  dependencies: dependencies,
  fortran_args: fortran_args,
  link_args: link_args,
  link_language: 'fortran'
)

proj_dep = declare_dependency(
  version: meson.project_version(),
  link_with: lib
)

exe = executable(
  meson.project_version(),
  sources, 
  fortran_args: fortran_args,
  link_args: link_args,
  dependencies: dependencies,
  install : true
)
exe = executable('peacemaker', sources, 
                 fortran_args : ['-cpp', '-fopenmp', '-flto', '-O3'],
                 link_args: ['-fopenmp', '-flto', '-O3'],
                 dependencies : [tomlf_dep],
                install : true)

subdir('test')
