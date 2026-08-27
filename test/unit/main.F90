!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! PROGRAM: brixs_unit_tests
!
! DESCRIPTION:
!> Driver for the Fortran unit tests (one call per src/ module).
!------------------------------------------------------------------------------
program brixs_unit_tests
  use test_mod_matmul, only : mod_matmul_test_driver
  use test_mod_blocks, only : mod_blocks_test_driver
  use test_mod_blocks_k, only : mod_blocks_k_test_driver
  use test_mod_io, only : mod_io_test_driver
  implicit none
  logical :: all_passed

  all_passed = .true.
  all_passed = mod_matmul_test_driver()   .and. all_passed
  all_passed = mod_blocks_test_driver()   .and. all_passed
  all_passed = mod_blocks_k_test_driver() .and. all_passed
  all_passed = mod_io_test_driver()       .and. all_passed

  if (.not. all_passed) then
    write(*,'(A)') 'One or more unit tests FAILED.'
    error stop 1
  end if
end program brixs_unit_tests
