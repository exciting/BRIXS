# Wraps the existing test/test.py regression suite (which runs the freshly
# built executables against test/data/**/*_ref.h5) as a single CTest test.
# The pybrixs_test conda environment (h5py, numpy) is expected to exist.
find_program(BASH_EXECUTABLE bash REQUIRED)

add_test(
  NAME brixs_regression
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.sh
)
set_tests_properties(brixs_regression PROPERTIES
  ENVIRONMENT "BRIXS_BIN_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
)
