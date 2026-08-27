# Fortran unit tests (vs. regression suite in test/test.py).
# Uses brixs_unit_test.F90, a small in-house assertion helper,
# following exciting's approach.
add_executable(brixs_unit_tests
  brixs_unit_test.F90
  test_mod_matmul.F90
  test_mod_blocks.F90
  test_mod_blocks_k.F90
  test_mod_io.F90
  main.F90
)
target_link_libraries(brixs_unit_tests PRIVATE brixs_core)

add_test(NAME brixs_unit_tests COMMAND brixs_unit_tests)
