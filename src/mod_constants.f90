module mod_constants
  implicit none
  private

  public :: hartree_to_ev

  !> Hartree-to-eV conversion factor
  real(8), parameter :: hartree_to_ev = 27.211d0

end module mod_constants
