module mod_rixs
  use mod_io
  use mod_matmul
  implicit none
  
  public transform2vector
  public transform2matrix
  contains
  ! Methodenbereich
  !-----------------------------------------------------------------------------
  subroutine transform2matrix(inkoulims,insmap,inarray,outarray)
    implicit none
    integer(4), intent(in) :: inkoulims(:,:)
    integer(4), intent(in) :: insmap(:,:)
    complex(8), intent(in) :: inarray(:)
    complex(8), allocatable, intent(out) :: outarray(:,:,:)

    ! local variables
    integer(4) :: lu,uu,lo,uo,nu,no,nk0,nkmax,i
    integer(4) :: dim_koulims(2),dim_smap(2),hamsiz,i1,i2,i3
  
    !get shapes
    dim_koulims=shape(inkoulims)
    dim_smap=shape(insmap)
    !determine sizes
    lu=inkoulims(1,1)
    uu=inkoulims(2,1)
    lo=inkoulims(3,1)
    uo=inkoulims(4,1)
    nu=uu-lu+1
    no=uo-lo+1
    nk0=insmap(3,1)
    nkmax=dim_koulims(2)
    hamsiz=dim_smap(2)

    ! allocate the output matrix
    if (allocated(outarray)) deallocate(outarray)
    allocate(outarray(nu,no,nkmax))

    ! loop over all transitions
    do i=1, hamsiz
      i1=insmap(1,i)-lu+1
      i2=insmap(2,i)-lo+1
      i3=insmap(3,i)-nk0+1
      outarray(i1,i2,i3)=inarray(i)
    end do
  end subroutine 

  !-----------------------------------------------------------------------------
  subroutine transform2vector(inkoulims,insmap,inarray,outarray)
    implicit none
    integer(4), intent(in) :: inkoulims(:,:)
    integer(4), intent(in) :: insmap(:,:)
    complex(8), intent(in) :: inarray(:,:,:)
    complex(8), allocatable, intent(out) :: outarray(:)

    ! local variables
    integer(4) :: lu,uu,lo,uo,nu,no,nk0,i
    integer(4) :: dim_koulims(2),dim_smap(2),nkmax,hamsiz,i1,i2,i3
  
    !get shapes
    dim_koulims=shape(inkoulims)
    dim_smap=shape(insmap)
    !determine sizes
    lu=inkoulims(1,1)
    uu=inkoulims(2,1)
    lo=inkoulims(3,1)
    uo=inkoulims(4,1)
    nu=uu-lu+1
    no=uo-lo+1
    nk0=insmap(3,1)
    nkmax=dim_koulims(2)
    hamsiz=dim_smap(2)

    ! allocate the output matrix
    if (allocated(outarray)) deallocate(outarray)
    allocate(outarray(hamsiz))

    ! loop over all transitions
    do i=1, hamsiz
      i1=insmap(1,i)-lu+1
      i2=insmap(2,i)-lo+1
      i3=insmap(3,i)-nk0+1
      outarray(i)=inarray(i1,i2,i3)
    end do
  end subroutine 
end module 
