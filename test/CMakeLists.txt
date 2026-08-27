# Wraps the existing test/test.py regression suite (which runs the built
# executables against test/data/**/*_ref.h5) as a single CTest test.
# The pybrixs_test conda environment (h5py, numpy) is expected to exist.
find_program(BASH_EXECUTABLE bash REQUIRED)

add_test(
  NAME brixs_regression
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.sh
)
set_tests_properties(brixs_regression PROPERTIES
  ENVIRONMENT "BRIXS_BIN_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
)

# Re-runs the suite under mpirun (2 ranks) so the #ifdef MPI code paths are
# tested. test.py's launch() reads BRIXS_MPI_LAUNCHER and prefixes it onto 
# every executable call.
if(BRIXS_ENABLE_MPI)
  add_test(
    NAME brixs_regression_mpi
    COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.sh
  )
  set_tests_properties(brixs_regression_mpi PROPERTIES
    ENVIRONMENT "BRIXS_BIN_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY};BRIXS_MPI_LAUNCHER=${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}"
  )
endif()

add_subdirectory(unit)
