testdrive_dep = dependency('test-drive', fallback: ['test-drive', 'testdrive_dep'])

# Um eine neue Testsuite zu implementieren, sind die folgenden
# Schritte nötig. Dabei sei <name> der Name der Testsuite.
#   1. Schreibe '<name>', in die tests Liste in diesem File.
#   2. kopiere die Datei template.f90~ und gib ihr den Namen test_<name>.f90
#   3. ergänze Block 1 und Block 2 in test/main.f90 um eine neue Zeile
#   4. ersetze in dem neuen test_<name>.f90 file den Platzhalter <name> (mit find and replace)
#   5. schreibe deine Tests

tests = [
  'pop_polynomial',
  'vol_polynomial',
  'partition_functions',
  'process_coordinates_record',
  'qce',
  'thermo',
  'auxiliary'
]

test_srcs = files(
  'main.f90',
  'test_pop_polynomial.f90',
  'test_vol_polynomial.f90',
  'test_partition_functions.f90',
  'test_process_coordinates_record.f90',
  'test_qce.f90',
  'test_thermo.f90',
  'test_auxiliary.f90'
)

foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [proj_dep, testdrive_dep],
)

test('all tests', tester)

foreach t : tests
  test(t, tester, args: t)
endforeach
