project(
  'cequalw2',
  'fortran',
  version: '5.0',
  default_options: ['fortran_std=f2008'],
)

compiler = meson.get_compiler('fortran')

sources = [
  'Diagenesis Bubbles Code 01.f90',
  'Diagenesis FFT Layer 01.f90',
  'Diagenesis Input 02.f90',
  'Diagenesis Input Files Read 01.f90',
  'Diagenesis Output 02.f90',
  'Diagenesis Sediment Flux Model 05.f90',
  'Diagenesis Sediment Model 03.f90',
  'TDGtarget.f90',
  'aerate.f90',
  'az.f90',
  'balances.F90',
  'date.f90',
  'density.f90',
  'endsimulation.F90',
  'envir_perf.f90',
  # 'fish-particle.f90',
  'fishhabitat.f90',
  'gas-transfer.f90',
  'heat-exchange.f90',
  'hydroinout.F90',
  'init-cond.F90',
  'init-geom.F90',
  'init-u-elws.f90',
  'init.F90',
  'input.F90',
  'layeraddsub.F90',
  'macrophyte-aux.f90',
  'output.f90',
  'outputa2w2tools.F90',
  'outputinitw2tools.F90',
  # 'particle.f90',
  'restart.f90',
  # 'screen_output_intel.f90', # <-- only needed if GUI is being built
  'shading.f90',
  'systdg.f90',
  'tdg.f90',
  'temperature.F90',
  'time-varying-data.f90',
  'transport.f90',
  'update.F90',
  'w2_4_win.f90',
  'w2modules.F90',
  'water-quality.f90',
  'waterbody.f90',
  'withdrawal.f90',
  'wqconstituents.F90',
  # 'Resource.fd',
]

add_project_arguments(
  [
    '-w',
    '-ffree-form',
    '-ffree-line-length-none',
    '-fallow-argument-mismatch',
    '-fdefault-real-8',
  ], language: 'fortran')

conf_data = configuration_data()
conf_data.set('MESON_VERSION', meson.version())
conf_data.set('GFORTRAN_VERSION', compiler.version())
conf_data.set('GIT_TAG', '@GIT_TAG@')

configured_file = configure_file(
  input: 'build_info.f90.in',
  output: 'build_info.f90.in2',
  configuration: conf_data,
)

vcs_tag(
  input: configured_file,
  replace_string: '@GIT_TAG@',
  output: 'build_info.f90',
  command: ['git', 'describe', '--tags', '--dirty', '--always'],
)


executable(
  'cequalw2',
  sources: sources,
)
