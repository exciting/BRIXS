!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: brixs_unit_test
!
! DESCRIPTION:
!> Minimal in-house unit-test helper (assert + report), following exciting's
!> unit_test_framework.
!------------------------------------------------------------------------------
module brixs_unit_test
  implicit none
  private
  public :: unit_test_type

  type :: unit_test_type
    integer :: n_assertions = 0
    integer :: n_failures = 0
  contains
    procedure :: assert
    procedure :: report
  end type unit_test_type

contains

  !> Log one assertion. On failure, the message is printed immediately.
  subroutine assert(this, condition, message)
    class(unit_test_type), intent(inout) :: this
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    this%n_assertions = this%n_assertions + 1
    if (.not. condition) then
      this%n_failures = this%n_failures + 1
      write(*,'(A)') 'FAILED: '//trim(message)
    end if
  end subroutine assert

  !> Print a summary; returns whether every assertion passed.
  function report(this, name) result(all_passed)
    class(unit_test_type), intent(in) :: this
    character(len=*), intent(in) :: name
    logical :: all_passed

    write(*,'(A,I0,A,I0)') trim(name)//': ', &
      this%n_assertions - this%n_failures, ' / ', this%n_assertions
    all_passed = (this%n_failures == 0)
  end function report

end module brixs_unit_test
