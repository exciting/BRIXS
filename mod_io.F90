module mod_io
  use mod_hdf5
  implicit none
  
  type :: io
    integer(4), allocatable :: smap(:,:), ismap(:,:,:), koulims(:,:)
    integer(4) :: hamsize, lu,uu,lo,uo, nk0, nkmax, nu, no
    real(8), allocatable :: evals(:)
    complex(8), allocatable :: eigvecs(:,:)
  end type io
  type :: input
    real(8), allocatable :: omega(:), omega2(:)
    real(8) :: broad, broad2
    logical :: oscstr
    integer :: nblocks
  end type
  

  public set_param
  public get_koulims
  public get_smap
  public get_ismap
  public get_evals
  public get_eigvecs
    
  contains
  ! Methodenbereich
  !-----------------------------------------------------------------------------
  subroutine get_koulims(object,fname)
    implicit none
    type(io), intent(inout) :: object
    character(len=1024), intent(in) :: fname
    !local variables
    integer(4) :: dims(2)
    character(len=1024) :: path, dsetname
    
    !get sizes of koulims
    path=trim(adjustl('eigvec-singlet-TDA-BAR-full/0001/parameters'))
    dsetname=trim(adjustl('koulims'))
    call hdf5_get_dims(trim(adjustl(fname)),path,dsetname,dims)
    !allocate output
    if (allocated(object%koulims)) deallocate(object%koulims)
    allocate(object%koulims(dims(1),dims(2)))
    ! get data
    call hdf5_read(trim(adjustl(fname)),path,dsetname,object%koulims(1,1),dims)
  end subroutine
  
  !-----------------------------------------------------------------------------
  subroutine get_smap(object,fname)
    implicit none
    type(io), intent(inout) :: object
    character(len=1024), intent(in) :: fname
    !local variables
    integer(4) :: dims(2)
    character(len=1024) :: path, dsetname
    
    !get sizes of koulims
    path='eigvec-singlet-TDA-BAR-full/0001/parameters'
    dsetname='smap'
    call hdf5_get_dims(fname,path,dsetname,dims)
    !allocate output
    if (allocated(object%smap)) deallocate(object%smap)
    allocate(object%smap(dims(1),dims(2)))
    ! get data
    call hdf5_read(fname,path,dsetname,object%smap(1,1),dims)
  end subroutine 

  !-----------------------------------------------------------------------------
  subroutine get_ismap(object)
    implicit none
    type(io), intent(inout) :: object
   !local variables
    integer(4) :: i, i1, i2, i3
    
    
    if (allocated(object%koulims) .and. allocated(object%smap)) then
      ! get parameters, just in case someone forgot to call it before
      call set_param(object)
      ! allocate ismap
      if (allocated(object%ismap)) deallocate(object%ismap)
      allocate(object%ismap(object%nu,object%no,object%nkmax))
      !fill in the inverse map
      do i=1,object%hamsize
        i1=object%smap(1,i)-object%lu+1
        i2=object%smap(2,i)-object%lo+1
        i3=object%smap(3,i)-object%nk0+1
        object%ismap(i1,i2,i3)=i
      end do
    end if
  end subroutine 
  !-----------------------------------------------------------------------------
  subroutine set_param(object)
    implicit none
    type(io), intent(inout) :: object
    !local variables
    integer(4), dimension(2) :: dim_koulims, dim_smap
    if ((allocated(object%koulims)) .and. (allocated(object%smap))) then
      !get shapes
      dim_koulims=shape(object%koulims)
      dim_smap=shape(object%smap)
      !determine sizes
      object%lu=object%koulims(1,1)
      object%uu=object%koulims(2,1)
      object%lo=object%koulims(3,1)
      object%uo=object%koulims(4,1)
      object%nu=object%uu-object%lu+1
      object%no=object%uo-object%lo+1
      object%nk0=object%smap(3,1)
      object%nkmax=dim_koulims(2)
      object%hamsize=dim_smap(2)
    else
      print *, 'koulims and smap have to be obtained from file before set_param can be called!'
    end if
  end subroutine 
  !-----------------------------------------------------------------------------
  subroutine get_evals(object,fname)
    implicit none
    type(io), intent(inout) :: object
    character(len=1024), intent(in) :: fname
    !local variables
    integer(4) :: dim_(1)
    character(len=1024) :: path, dsetname
    
    !get sizes of koulims
    path='eigvec-singlet-TDA-BAR-full/0001'
    dsetname='evals'
    call hdf5_get_dims(fname,path,dsetname,dim_)
    !allocate output
    if (allocated(object%evals)) deallocate(object%evals)
    allocate(object%evals(dim_(1)))
    ! get data
    call hdf5_read(fname,path,dsetname,object%evals(1),dim_)
  end subroutine
  
  !-----------------------------------------------------------------------------
  subroutine get_eigvecs(object,fname)
    implicit none
    type(io), intent(inout) :: object
    character(len=1024), intent(in) :: fname
    !local variables
    complex(8), allocatable :: eigvec_(:)
    integer(4) :: dim_(1), i, dims_(2)
    character(len=1024) :: path, dsetname
    character(256) :: ci
    
    !get size of eigvecs
    ! IMPORTANT: Here, I assume that all of the BSE eigvecs are included in the
    !            hdf5 file

    path='eigvec-singlet-TDA-BAR-full/0001'
    dsetname='evals'
    call hdf5_get_dims(fname,path,dsetname,dim_)
    !allocate output
    if (allocated(object%eigvecs)) deallocate(object%eigvecs)
    allocate(object%eigvecs(dim_(1),dim_(1)))
    ! get data
    do i=1, dim_(1)
      write(ci, '(I8.8)') i
      path='eigvec-singlet-TDA-BAR-full/0001/rvec'
      dsetname=trim(adjustl(ci))
      ! Get dimension of eigvec for given lambda
      call hdf5_get_dims(fname,path,ci,dims_)
      ! Allocate intermediate eigenvector array
      ! The first dimension is 2, since this is a complex array
      if (allocated(eigvec_)) deallocate(eigvec_)
      !allocate(eigvec_(dims_(2),dims_(3)))
      allocate(eigvec_(dims_(2)))
      ! Get data
      call hdf5_read(fname,path,dsetname,eigvec_(1),shape(eigvec_))
      ! Write data to final array
      object%eigvecs(:,i)=eigvec_(:)
    end do
    deallocate(eigvec_)
  end subroutine 
  
  !-----------------------------------------------------------------------------
  subroutine read_inputfile(object,fname)
    use modmpi
    implicit none
    type(input), intent(out) :: object
    character(*), intent(in) :: fname
    ! local variables
    integer :: line, ios, w, pos
    character(256) :: buffer,label
    integer, parameter :: fh = 15
    real(8) :: inter(3), inter2(3)
    real(8) :: broad_, broad2_
    integer :: nblocks_, ierr
    logical :: oscstr_


    ! only root reads the input file
    if (mpiglobal%rank .eq. 0) then
      ! basics taken from https://jblevins.org/log/control-file 
      line=0
      ios=0
      open(fh, file=trim(adjustl(fname)))
      do while (ios == 0)
        read(fh, '(A)', iostat=ios) buffer
        if (ios == 0) then
          line = line + 1
          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('omega')
            read(buffer, *, iostat=ios) inter
          case ('omega2')
            read(buffer, *, iostat=ios) inter2
          case ('broad')
            read(buffer, *, iostat=ios) broad_
          case ('broad2')
            read(buffer, *, iostat=ios) broad2_
          case ('do_oscstr')
            oscstr_=.true.
          case ('nblocks')
            read(buffer, *, iostat=ios) nblocks_
          case default
            print *, 'Skipping invalid label at line', line
          end select
        end if
      end do
    end if
#ifdef MPI
    ! broadcast input parameters to everybody
    call mpi_bcast(inter,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(inter2,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(broad_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(broad2_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(oscstr_,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nblocks_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#end if
    
    ! calculate input parameters from read
    object%broad=broad_
    object%broad2=broad2_
    object%nblocks=nblocks_
    if (allocated(object%omega)) deallocate(object%omega)
    if (allocated(object%omega2)) deallocate(object%omega2)
    allocate(object%omega(int(inter(3))))
    allocate(object%omega2(int(inter2(3))))
    
    do w=1,int(inter(3))
      object%omega(w)=(inter(2)-inter(1))/(inter(3)-1.0d0)*(w-1) + inter(1)
    end do
    do w=1,int(inter2(3))
      object%omega2(w)=(inter2(2)-inter2(1))/(inter2(3)-1.0d0)*(w-1) + inter2(1)
    end do
  end subroutine read_inputfile  
  !-----------------------------------------------------------------------------

  subroutine inspect_h5(object)
    implicit none
    type(io), intent(in) :: object
    ! local variables
    integer::  lo,ho,lu,hu

    ! determine lowest and highest occupied/unoccupied band included
    ! these could be different from the ones in the RIXS calculation 
    ! for metallic systems
    lu=minval(object%koulims(1,:))
    hu=maxval(object%koulims(2,:))
    lo=minval(object%koulims(3,:))
    ho=maxval(object%koulims(4,:))
    write(*,*) '************Transitions************'
    write(*,*) 'No. of k-points:', object%nkmax
    if (.not.(object%nkmax .eq. maxval(object%smap(3,:)))) then
      write(*,*) 'Not all k-points included in BSE calculation!!'
    end if
    if ((object%lo .eq. lo) .and. (object%uo .eq. ho)) then
      write(*,*) 'range of occupied bands:', object%lo, object%uo
      write(*,*) 'no. of occupied bands:', object%no
    else
      write(*,*) 'WRONG enumeration of occupied states: (', object%lo, ',', object%uo, ') != (', lo, ',', ho, ')' 
    end if
    if ((object%lu .eq. lu) .and. (object%uu .eq. hu)) then
      write(*,*) 'range of unoccupied bands:', object%lu, object%uu
      write(*,*) 'no. of unoccupied bands:', object%nu
    else
      write(*,*) 'WRONG enumeration of unoccupied states: (', object%lu, ',', object%uu, ') != (', lu, ',', hu, ')' 
    end if
    write(*,*) 'total number of transitions:', object%hamsize
    
  end subroutine
end module

