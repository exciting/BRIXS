!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: test_mod_blocks_k
!
! DESCRIPTION:
!> Unit tests for the k-point-resolved mod_blocks transforms, used by
!> rixs_coherence for per-k evaluation: a single-k slice, and a
!> matrix -> vector -> matrix round trip.
!------------------------------------------------------------------------------
module test_mod_blocks_k
  use brixs_unit_test, only : unit_test_type
  use mod_blocks, only : block2d
  use mod_blocks_k, only : transform_matrix2matrix_k, transform_matrix2vector_k
  use mod_io, only : io
  implicit none
  private
  public :: mod_blocks_k_test_driver

contains

  function mod_blocks_k_test_driver() result(all_passed)
    logical :: all_passed
    type(unit_test_type) :: test

    call test_transform_matrix2matrix_k_single_k_slice(test)
    call test_transform_matrix2matrix_vector_k_roundtrip(test)

    all_passed = test%report('mod_blocks_k')
  end function mod_blocks_k_test_driver

  !> Block holds only k=2's transition (real callers pre-select via
  !> get_transition_range); must land in its own window and drop the k-axis.
  subroutine test_transform_matrix2matrix_k_single_k_slice(test)
    type(unit_test_type), intent(inout) :: test
    integer(4) :: koulims(4,2), smap(3,2)
    type(block2d) :: bl2d
    complex(8), allocatable :: out3d(:,:,:)
    real(8), parameter :: tol = 1d-12

    ! k=1: unocc index 1..1, occ index 1..1
    ! k=2: unocc index 3..4, occ index 2..3 -- disjoint from k=1's window
    koulims(:,1) = [1, 1, 1, 1]
    koulims(:,2) = [3, 4, 2, 3]

    ! full transition map (both k's); transition 2 is (unocc=3, occ=2, k=2)
    smap(:,1) = [1, 1, 1]
    smap(:,2) = [3, 2, 2]

    ! the block passed in holds only transition 2 (k=2's one transition)
    bl2d%nk = 1
    bl2d%blocksize = [1, 1]
    bl2d%il = 2
    bl2d%iu = 2
    bl2d%offset = [1, 0]
    allocate(bl2d%zcontent(1,1))
    bl2d%zcontent(1,1) = (2d0, 0d0)

    call transform_matrix2matrix_k(koulims, 2_4, smap, bl2d, out3d)

    call test%assert(size(out3d,1) == 4, &
      "transform_matrix2matrix_k: nu spans the union of the per-k unocc windows")
    call test%assert(size(out3d,2) == 3, &
      "transform_matrix2matrix_k: no spans the union of the per-k occ windows")
    call test%assert(size(out3d,3) == 1, &
      "transform_matrix2matrix_k: k-axis is dropped (single k-point slice)")
    call test%assert(abs(out3d(3,2,1) - (2d0,0d0)) < tol, &
      "transform_matrix2matrix_k: transition landed in k=2's own band window")
    call test%assert(abs(sum(abs(out3d)) - 2d0) < tol, &
      "transform_matrix2matrix_k: no stray values written outside the transition")
  end subroutine test_transform_matrix2matrix_k_single_k_slice

  !> transform_matrix2matrix_k and its inverse transform_matrix2vector_k
  !> should round-trip exactly.
  subroutine test_transform_matrix2matrix_vector_k_roundtrip(test)
    type(unit_test_type), intent(inout) :: test
    integer(4) :: koulims(4,2), smap(3,3)
    type(io) :: input
    type(block2d) :: bl2d_in, bl2d_out
    complex(8), allocatable :: out3d(:,:,:)
    real(8), parameter :: tol = 1d-12

    koulims(:,1) = [1, 2, 1, 1]
    koulims(:,2) = [1, 2, 1, 2]

    ! transitions 1,2 at k=1; transition 3 at k=2 (the one under test)
    smap(:,1) = [1, 1, 1]
    smap(:,2) = [2, 1, 1]
    smap(:,3) = [1, 2, 2]

    input%koulims = koulims
    input%smap = smap

    ! block holding only transition 3 (k=2), two exciton columns (lambda)
    bl2d_in%nk = 1
    bl2d_in%blocksize = [1, 2]
    bl2d_in%il = 3
    bl2d_in%iu = 3
    bl2d_in%offset = [2, 0]
    allocate(bl2d_in%zcontent(1,2))
    bl2d_in%zcontent(1,1) = (3d0, -1d0)
    bl2d_in%zcontent(1,2) = (0d0, 4d0)

    call transform_matrix2matrix_k(koulims, 2_4, smap, bl2d_in, out3d)

    bl2d_out%blocksize = [1, 2]
    bl2d_out%il = 3
    bl2d_out%iu = 3
    bl2d_out%offset = [2, 0]

    call transform_matrix2vector_k(input, 2_4, out3d, bl2d_out)

    call test%assert(abs(bl2d_out%zcontent(1,1) - bl2d_in%zcontent(1,1)) < tol, &
      "transform_matrix2{matrix,vector}_k: roundtrip lambda=1")
    call test%assert(abs(bl2d_out%zcontent(1,2) - bl2d_in%zcontent(1,2)) < tol, &
      "transform_matrix2{matrix,vector}_k: roundtrip lambda=2")
  end subroutine test_transform_matrix2matrix_vector_k_roundtrip

end module test_mod_blocks_k
