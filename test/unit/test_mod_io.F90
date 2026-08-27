!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: test_mod_io
!
! DESCRIPTION:
!> Unit tests for mod_io's: set_param, get_transition_range, get_ismap.
!------------------------------------------------------------------------------
module test_mod_io
  use brixs_unit_test, only : unit_test_type
  use mod_io, only : io, set_param, get_transition_range, get_ismap
  implicit none
  private
  public :: mod_io_test_driver

contains

  function mod_io_test_driver() result(all_passed)
    logical :: all_passed
    type(unit_test_type) :: test

    call test_set_param_per_k_window(test)
    call test_transition_range_and_ismap(test)

    all_passed = test%report('mod_io')
  end function mod_io_test_driver

  !> set_param applies the same koulims union-window logic independently
  subroutine test_set_param_per_k_window(test)
    type(unit_test_type), intent(inout) :: test
    type(io) :: object

    allocate(object%koulims(4,2))
    object%koulims(:,1) = [1, 1, 1, 1]
    object%koulims(:,2) = [3, 4, 2, 3]

    allocate(object%smap(3,4))
    object%smap(:,1) = [1, 1, 1]
    object%smap(:,2) = [1, 1, 1]
    object%smap(:,3) = [3, 2, 2]
    object%smap(:,4) = [3, 2, 2]

    call set_param(object)

    call test%assert(object%lu == 1, "set_param: lu")
    call test%assert(object%uu == 4, "set_param: uu")
    call test%assert(object%lo == 1, "set_param: lo")
    call test%assert(object%uo == 3, "set_param: uo")
    call test%assert(object%nu == 4, "set_param: nu")
    call test%assert(object%no == 3, "set_param: no")
    call test%assert(object%nk0 == 1, "set_param: nk0")
    call test%assert(object%nkmax == 2, "set_param: nkmax")
    call test%assert(object%hamsize == 4, "set_param: hamsize")
    call test%assert(object%globalk == 12, "set_param: globalk = nu*no")
    call test%assert(object%global == 4, "set_param: global = hamsize")
  end subroutine test_set_param_per_k_window

  !> Three k-points with nk0 /= 1 (global k=5..7), exercising the nk0
  !> offset in both routines.
  subroutine test_transition_range_and_ismap(test)
    type(unit_test_type), intent(inout) :: test
    type(io) :: object
    integer(4) :: il, iu

    allocate(object%koulims(4,3))
    object%koulims(:,1) = [1, 2, 1, 2]
    object%koulims(:,2) = [1, 2, 1, 2]
    object%koulims(:,3) = [1, 2, 1, 2]

    allocate(object%smap(3,6))
    object%smap(:,1) = [1, 1, 5]
    object%smap(:,2) = [2, 1, 5]
    object%smap(:,3) = [1, 2, 6]
    object%smap(:,4) = [2, 2, 6]
    object%smap(:,5) = [1, 1, 7]
    object%smap(:,6) = [2, 1, 7]

    call set_param(object)
    call test%assert(object%nk0 == 5, &
      "set_param: nk0 should be the first k-point's global index")

    ! relative k-index 2 (global k=6) -> transitions 3..4
    call get_transition_range(object, 2, 2, il, iu)
    call test%assert(il == 3, "get_transition_range: il for k=6 alone")
    call test%assert(iu == 4, "get_transition_range: iu for k=6 alone")

    ! relative k-index range 2..3 (global k=6..7) -> transitions 3..6
    call get_transition_range(object, 2, 3, il, iu)
    call test%assert(il == 3, "get_transition_range: il for k=6..7")
    call test%assert(iu == 6, "get_transition_range: iu for k=6..7")

    call get_ismap(object)
    call test%assert(size(object%ismap,1) == 2, "get_ismap: dim 1 (nu)")
    call test%assert(size(object%ismap,2) == 2, "get_ismap: dim 2 (no)")
    call test%assert(size(object%ismap,3) == 3, "get_ismap: dim 3 (nkmax)")

    ! spot-check the inverse map round-trips back to the transition index
    call test%assert(object%ismap(1,1,1) == 1, "get_ismap: (unocc=1,occ=1,k=5)")
    call test%assert(object%ismap(2,2,2) == 4, "get_ismap: (unocc=2,occ=2,k=6)")
    call test%assert(object%ismap(2,1,3) == 6, "get_ismap: (unocc=2,occ=1,k=7)")
  end subroutine test_transition_range_and_ismap

end module test_mod_io
