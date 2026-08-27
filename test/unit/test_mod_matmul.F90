!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: test_mod_matmul
!
! DESCRIPTION:
!> Unit tests for matprod (mod_matmul)
!------------------------------------------------------------------------------
module test_mod_matmul
  use brixs_unit_test, only : unit_test_type
  use mod_matmul, only : matprod
  implicit none
  private
  public :: mod_matmul_test_driver

contains

  function mod_matmul_test_driver() result(all_passed)
    logical :: all_passed
    type(unit_test_type) :: test

    call test_matmat(test)
    call test_matvec(test)
    call test_vecvec(test)

    all_passed = test%report('mod_matmul')
  end function mod_matmul_test_driver

  subroutine test_matmat(test)
    type(unit_test_type), intent(inout) :: test
    complex(8) :: a(2,2), b(2,2), c(2,2)
    real(8), parameter :: tol = 1d-12

    ! a = [[1,2],[3,4]], b = [[5,6],[7,8]] => a*b = [[19,22],[43,50]]
    a = reshape([(1d0,0d0),(3d0,0d0),(2d0,0d0),(4d0,0d0)], [2,2])
    b = reshape([(5d0,0d0),(7d0,0d0),(6d0,0d0),(8d0,0d0)], [2,2])

    call matprod(a, b, c)

    call test%assert(abs(c(1,1) - (19d0,0d0)) < tol, "matprod_matmat: C(1,1)")
    call test%assert(abs(c(1,2) - (22d0,0d0)) < tol, "matprod_matmat: C(1,2)")
    call test%assert(abs(c(2,1) - (43d0,0d0)) < tol, "matprod_matmat: C(2,1)")
    call test%assert(abs(c(2,2) - (50d0,0d0)) < tol, "matprod_matmat: C(2,2)")
  end subroutine test_matmat

  subroutine test_matvec(test)
    type(unit_test_type), intent(inout) :: test
    complex(8) :: a(2,2), v(2), c(2)
    real(8), parameter :: tol = 1d-12

    ! a = [[1,2],[3,4]], v = [1,1] => a*v = [3,7]
    a = reshape([(1d0,0d0),(3d0,0d0),(2d0,0d0),(4d0,0d0)], [2,2])
    v = [(1d0,0d0),(1d0,0d0)]

    call matprod(a, v, c)

    call test%assert(abs(c(1) - (3d0,0d0)) < tol, "matprod_matvec: c(1)")
    call test%assert(abs(c(2) - (7d0,0d0)) < tol, "matprod_matvec: c(2)")
  end subroutine test_matvec

  subroutine test_vecvec(test)
    type(unit_test_type), intent(inout) :: test
    complex(8) :: a(2), b(2), c
    real(8), parameter :: tol = 1d-12

    ! matprod_vecvec computes the conjugated inner product conjg(a).b
    a = [(1d0,1d0),(2d0,-1d0)]
    b = [(2d0,0d0),(1d0,3d0)]

    call matprod(a, b, c)

    call test%assert(abs(c - (1d0,5d0)) < tol, "matprod_vecvec: conjg(a).b")
  end subroutine test_vecvec

end module test_mod_matmul
