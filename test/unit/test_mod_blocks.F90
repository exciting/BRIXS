!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: test_mod_blocks
!
! DESCRIPTION:
!> Unit tests for mod_blocks: transform_matrix2matrix (koulims union-window
!> fix), block partitioning (nofblock/firstofblock/lastofblock/
!> make_block1d/make_block2d), and apply_occupation_factors.
!------------------------------------------------------------------------------
module test_mod_blocks
  use brixs_unit_test, only : unit_test_type
  use mod_blocks, only : block1d, block2d, transform_matrix2matrix, &
    &                    nofblock, firstofblock, lastofblock, &
    &                    make_block1d, make_block2d, apply_occupation_factors
  implicit none
  private
  public :: mod_blocks_test_driver

contains

  function mod_blocks_test_driver() result(all_passed)
    logical :: all_passed
    type(unit_test_type) :: test

    call test_transform_matrix2matrix_per_k_window(test)
    call test_partition_contiguous_and_covers_range(test)
    call test_partition_empty_block(test)
    call test_make_block1d_auto(test)
    call test_make_block1d_override(test)
    call test_make_block2d_override(test)
    call test_apply_occupation_factors(test)

    all_passed = test%report('mod_blocks')
  end function mod_blocks_test_driver

  subroutine test_transform_matrix2matrix_per_k_window(test)
    type(unit_test_type), intent(inout) :: test
    integer(4) :: koulims(4,2), smap(3,2)
    type(block2d) :: bl2d
    complex(8), allocatable :: outbl4d(:,:,:,:)
    real(8), parameter :: tol = 1d-12

    ! k=1: unocc index 1..1, occ index 1..1
    ! k=2: unocc index 3..4, occ index 2..3 -- disjoint from k=1's window
    koulims(:,1) = [1, 1, 1, 1]
    koulims(:,2) = [3, 4, 2, 3]

    ! transition 1: (unocc=1, occ=1, k=1); transition 2: (unocc=3, occ=2, k=2)
    smap(:,1) = [1, 1, 1]
    smap(:,2) = [3, 2, 2]

    bl2d%nk = 2
    bl2d%blocksize = [2, 1]
    bl2d%il = 1
    bl2d%iu = 2
    bl2d%offset = [0, 0]
    bl2d%k1l = 1
    allocate(bl2d%zcontent(2,1))
    bl2d%zcontent(1,1) = (1d0, 0d0)
    bl2d%zcontent(2,1) = (2d0, 0d0)

    call transform_matrix2matrix(koulims, smap, bl2d, outbl4d)

    ! union of both k windows (nu=4, no=3); k=1-only sizing (pre-fix) would
    ! give a 1x1 window and put transition 2 out of bounds.
    call test%assert(size(outbl4d,1) == 4, &
      "transform_matrix2matrix: nu must span the union of the per-k unocc windows")
    call test%assert(size(outbl4d,2) == 3, &
      "transform_matrix2matrix: no must span the union of the per-k occ windows")
    call test%assert(abs(outbl4d(1,1,1,1) - (1d0,0d0)) < tol, &
      "transform_matrix2matrix: k=1 transition landed in its own band window")
    call test%assert(abs(outbl4d(3,2,2,1) - (2d0,0d0)) < tol, &
      "transform_matrix2matrix: k=2 transition landed in its own (shifted) band window")
    call test%assert(abs(sum(abs(outbl4d)) - 3d0) < tol, &
      "transform_matrix2matrix: no stray values written outside the two transitions")
  end subroutine test_transform_matrix2matrix_per_k_window

  !> Blocks must partition [1, globalsize] exactly: no gaps/overlaps,
  !> sizes sum to globalsize.
  subroutine test_partition_contiguous_and_covers_range(test)
    type(unit_test_type), intent(inout) :: test

    call check_partition(test, 10, 3)
    call check_partition(test, 7, 7)
    call check_partition(test, 1, 4)
  end subroutine test_partition_contiguous_and_covers_range

  subroutine check_partition(test, globalsize, nblocks)
    type(unit_test_type), intent(inout) :: test
    integer(4), intent(in) :: globalsize, nblocks
    integer(4) :: i, total

    call test%assert(firstofblock(1, globalsize, nblocks) == 1, &
      "partition: first block must start at index 1")

    total = 0
    do i = 1, min(nblocks, globalsize)
      total = total + nofblock(i, globalsize, nblocks)
      if (i > 1) then
        call test%assert(firstofblock(i, globalsize, nblocks) == &
          lastofblock(i-1, globalsize, nblocks) + 1, &
          "partition: blocks must be contiguous, no gap/overlap")
      end if
    end do
    call test%assert(total == globalsize, &
      "partition: block sizes must sum to the global size")
    call test%assert(lastofblock(min(nblocks,globalsize), globalsize, nblocks) == globalsize, &
      "partition: last occupied block must end at globalsize")
  end subroutine check_partition

  !> More blocks than elements: trailing blocks are empty, using the
  !> (0, -1) sentinel.
  subroutine test_partition_empty_block(test)
    type(unit_test_type), intent(inout) :: test
    integer(4), parameter :: globalsize = 2, nblocks = 5

    call test%assert(firstofblock(1, globalsize, nblocks) == 1, "partition: block 1 first")
    call test%assert(lastofblock(1, globalsize, nblocks) == 1, "partition: block 1 last")
    call test%assert(firstofblock(2, globalsize, nblocks) == 2, "partition: block 2 first")
    call test%assert(lastofblock(2, globalsize, nblocks) == 2, "partition: block 2 last")
    call test%assert(firstofblock(3, globalsize, nblocks) == 0, &
      "partition: empty block must report first=0")
    call test%assert(lastofblock(3, globalsize, nblocks) == -1, &
      "partition: empty block must report last=-1")
  end subroutine test_partition_empty_block

  subroutine test_make_block1d_auto(test)
    type(unit_test_type), intent(inout) :: test
    type(block1d) :: b
    integer(4) :: il_expected, iu_expected

    il_expected = firstofblock(2_4, 10_4, 3_4)
    iu_expected = lastofblock(2_4, 10_4, 3_4)
    b = make_block1d(2_4, 10_4, 3_4)

    call test%assert(b%il == il_expected, "make_block1d: il from auto partition")
    call test%assert(b%iu == iu_expected, "make_block1d: iu from auto partition")
    call test%assert(b%blocksize == iu_expected - il_expected + 1, &
      "make_block1d: blocksize = iu - il + 1")
    call test%assert(b%offset == il_expected - 1, "make_block1d: offset = il - 1")
    call test%assert(b%nblocks == 3, "make_block1d: nblocks")
    call test%assert(b%global == 10, "make_block1d: global")
    call test%assert(b%id == 2, "make_block1d: id")
  end subroutine test_make_block1d_auto

  !> The il/iu override (non-equilibrium occupations) must bypass the
  !> auto-partition entirely.
  subroutine test_make_block1d_override(test)
    type(unit_test_type), intent(inout) :: test
    type(block1d) :: b

    b = make_block1d(1_4, 10_4, 3_4, il=2_4, iu=6_4)

    call test%assert(b%il == 2, &
      "make_block1d: il must come from the override, not auto partition")
    call test%assert(b%iu == 6, "make_block1d: iu must come from the override")
    call test%assert(b%blocksize == 5, "make_block1d: blocksize = iu - il + 1")
    call test%assert(b%offset == 1, "make_block1d: offset = il - 1")
  end subroutine test_make_block1d_override

  subroutine test_make_block2d_override(test)
    type(unit_test_type), intent(inout) :: test
    type(block2d) :: b
    integer(4) :: jl_expected, ju_expected

    jl_expected = firstofblock(2_4, 6_4, 2_4)
    ju_expected = lastofblock(2_4, 6_4, 2_4)
    b = make_block2d(1_4, 6_4, 2_4, 6_4, 2_4, il=1_4, iu=6_4)

    call test%assert(b%il == 1, "make_block2d: il must come from the override")
    call test%assert(b%iu == 6, "make_block2d: iu must come from the override")
    call test%assert(b%jl == jl_expected, &
      "make_block2d: jl must still come from the auto partition of the 2nd dim")
    call test%assert(b%ju == ju_expected, &
      "make_block2d: ju must still come from the auto partition of the 2nd dim")
    call test%assert(b%blocksize(1) == 6, "make_block2d: blocksize(1) = iu - il + 1")
    call test%assert(b%blocksize(2) == ju_expected - jl_expected + 1, &
      "make_block2d: blocksize(2) = ju - jl + 1")
  end subroutine test_make_block2d_override

  !> Scales rows by the *global* transition index (offset + local row);
  !> nonzero offset + decoy factors catch an off-by-offset bug.
  subroutine test_apply_occupation_factors(test)
    type(unit_test_type), intent(inout) :: test
    type(block2d) :: bl2d
    real(8) :: occupation_factors(6)
    real(8), parameter :: tol = 1d-12

    bl2d%blocksize = [3, 2]
    bl2d%offset = [2, 0]
    allocate(bl2d%zcontent(3,2))
    bl2d%zcontent(1,:) = [(1d0,0d0), (0d0,1d0)]
    bl2d%zcontent(2,:) = [(2d0,0d0), (0d0,2d0)]
    bl2d%zcontent(3,:) = [(3d0,0d0), (0d0,3d0)]

    ! indices 3,4,5 belong to this block (offset=2); 1,2,6 must never be read
    occupation_factors = [99d0, 99d0, 0.5d0, 2d0, -1d0, 99d0]

    call apply_occupation_factors(bl2d, occupation_factors)

    call test%assert(abs(bl2d%zcontent(1,1) - (0.5d0,0d0)) < tol, &
      "apply_occupation_factors: row 1 (global index 3) real component")
    call test%assert(abs(bl2d%zcontent(1,2) - (0d0,0.5d0)) < tol, &
      "apply_occupation_factors: row 1 (global index 3) imaginary component")
    call test%assert(abs(bl2d%zcontent(2,1) - (4d0,0d0)) < tol, &
      "apply_occupation_factors: row 2 (global index 4) real component")
    call test%assert(abs(bl2d%zcontent(2,2) - (0d0,4d0)) < tol, &
      "apply_occupation_factors: row 2 (global index 4) imaginary component")
    call test%assert(abs(bl2d%zcontent(3,1) - (-3d0,0d0)) < tol, &
      "apply_occupation_factors: row 3 (global index 5) real component")
    call test%assert(abs(bl2d%zcontent(3,2) - (0d0,-3d0)) < tol, &
      "apply_occupation_factors: row 3 (global index 5) imaginary component")
  end subroutine test_apply_occupation_factors

end module test_mod_blocks
