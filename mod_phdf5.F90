module mod_phdf5
  interface phdf5_write
    module procedure phdf5_write_d, &
        &            phdf5_write_z
  end interface
  interface phdf5_read
    module procedure phdf5_read_d, &
        &            phdf5_read_z, &
        &            phdf5_read_i
  end interface
  ! initialize and finalize interface
  public phdf5_initialize
  public phdf5_finalize
  ! open, create and close a file
  public phdf5_open_file
  public phdf5_create_file
  public phdf5_close_file
  ! create group
  public phdf5_create_group
  ! prepare and finalize MPI read and write of datasets
  public phdf5_setup_write
  public phdf5_cleanup
  ! determine size of dataset
  public phdf5_get_dims
contains

!-------------------------------------------------------------------------------
  subroutine phdf5_initialize
    use hdf5
    use mpi
    implicit none
    ! local variables
    integer :: ierr
    character*100 :: errmsg
    
    call h5open_f(ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_initialize): h5open_f returned ",I6)')ierr
      goto 10
    endif    
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    stop
  end subroutine

!-------------------------------------------------------------------------------
  subroutine phdf5_create_file(fname,fparallel,file_id,comm)
    use hdf5
    use mpi
    implicit none
    character(*), intent(in) :: fname
    logical, intent(in) :: fparallel
    integer(hid_t), intent(out) :: file_id
    integer, optional :: comm
    ! local variables
    integer (hid_t) :: plist_id
    integer :: ierr
    integer :: info
    character*100 :: errmsg
    ! MPI File creation & Global Access
    if (fparallel) then
      ! set mpi info object
      info=MPI_INFO_NULL
      ! create file access property list w/ parallel IO access
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id, ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5pcreate_f returned ",I6)')ierr
        goto 10
      endif    
      call h5pset_fapl_mpio_f(plist_id,comm,info,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5pset_fapl_mpio_f returned ",I6)')ierr
        goto 10
      endif    
      ! create the file collectively
      call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,file_id,ierr,access_prp=plist_id)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5fcreate_f returned ",I6)')ierr
        goto 10
      endif    
      ! close property list
      call h5pclose_f(plist_id,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5pclose_f returned ",I6)')ierr
        goto 10
      endif    
    ! Serial Access
    else
      call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,file_id,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5fcreate_f returned ",I6)')ierr
        goto 10
      endif    
    end if
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  fname : ",A)')trim(fname)
    stop
  end subroutine

!-------------------------------------------------------------------------------
  subroutine phdf5_open_file(fname,fparallel,file_id,comm)
    use hdf5
    use mpi
    implicit none
    character(*), intent(in) :: fname
    logical, intent(in) :: fparallel
    integer(hid_t), intent(out) :: file_id
    integer, optional :: comm
    ! local variables
    integer (hid_t) :: plist_id
    integer :: ierr
    integer :: info
    character*100 :: errmsg
    ! MPI File creation & Global Access
    if (fparallel) then
      ! set mpi info object
      info=MPI_INFO_NULL
      ! create file access property list w/ parallel IO access
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id, ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5pcreate_f returned ",I6)')ierr
        goto 10
      endif    
      call h5pset_fapl_mpio_f(plist_id,comm,info,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5pset_fapl_mpio_f returned ",I6)')ierr
        goto 10
      endif    
      ! create the file collectively
      call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,file_id,ierr,access_prp=plist_id)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5fcreate_f returned ",I6)')ierr
        goto 10
      endif    
      ! close property list
      call h5pclose_f(plist_id,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5pclose_f returned ",I6)')ierr
        goto 10
      endif    
    ! Serial Access
    else
      call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,file_id,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(phdf5_create_file): h5fcreate_f returned ",I6)')ierr
        goto 10
      endif    
    end if
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  fname : ",A)')trim(fname)
    stop
  end subroutine
!-------------------------------------------------------------------------------
  subroutine phdf5_close_file(file_id)
    use hdf5
    implicit none
    integer(hid_t), intent(in) :: file_id
    ! local variable
    integer :: ierr
    character*100 :: errmsg
    
    call h5fclose_f(file_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_close_file): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  file_id : ",I4)')file_id
    stop
  
  end subroutine 

!-------------------------------------------------------------------------------    
  subroutine phdf5_create_group(file_id,path,gname)
    use hdf5
    implicit none
    integer(hid_t), intent(in) :: file_id
    character(*), intent(in) :: path
    character(*), intent(in) :: gname
    integer(hid_t) :: h5_group_id,h5_new_group_id
    integer :: ierr

    call h5gopen_f(file_id,trim(path),h5_group_id,ierr)
    if (ierr.ne.0) then
      write(*,'("Error(hdf5_create_group): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gcreate_f(h5_group_id,trim(gname),h5_new_group_id,ierr)
    if (ierr.ne.0) then
      write(*,'("Error(hdf5_create_group): h5gcreate_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(h5_new_group_id,ierr)
    if (ierr.ne.0) then
      write(*,'("Error(hdf5_create_group): h5gclose_f for the new group returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(h5_group_id,ierr)
    if (ierr.ne.0) then
      write(*,'("Error(hdf5_create_group): h5gclose_f for the existing path returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'("  file_id : ",I4)')file_id
    write(*,'("  path  : ",A)')trim(path)
    write(*,'("  gname : ",A)')trim(gname)  
    stop
    end subroutine

    !-------------------------------------------------------------------------------
    logical function phdf5_exist_group(file_id,path,gname)
        ! Check if group with the given name exists
        use hdf5
        implicit none
        integer(hid_t), intent(in) :: file_id
        character(*), intent(in) :: path
        character(*), intent(in) :: gname
        integer(hid_t):: h5_group_id
        integer ierr
 
        call h5gopen_f(file_id,trim(path),h5_group_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(phdf5_exist_group): h5gopen_f returned ",I6)')ierr
          goto 10
        endif
        call h5lexists_f(h5_group_id,trim(gname),phdf5_exist_group,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(phdf5_exist_group): h5lexists_f returned ",I6)')ierr
          goto 10
        endif
        call h5gclose_f(h5_group_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(phdf5_exist_group): h5fclose_f returned ",I6)')ierr
          goto 10
        endif
        return
        10 continue
        write(*,'("  file_id : ",I4)')file_id
        write(*,'("  path  : ",A)')trim(path)
        write(*,'("  gname : ",A)')trim(gname)  
        stop
    end function phdf5_exist_group

!-------------------------------------------------------------------------------
  subroutine phdf5_finalize
    use hdf5
    implicit none
    ! local variable
    integer :: ierr
    character*100 :: errmsg

    call h5close_f(ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_finalize): h5close_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    stop
  
  end subroutine 
!-------------------------------------------------------------------------------
  subroutine phdf5_setup_write(ndims,dims,fcomplex,dname,path,file_id,dataset_id)
    use hdf5
    implicit none
    integer, intent(in) :: ndims
    integer, dimension(ndims), intent(in) :: dims
    logical, intent(in) :: fcomplex
    character(*), intent(in) :: dname, path
    integer(hid_t), intent(in) :: file_id
    integer(hid_t), intent(out) :: dataset_id
    !local variables
    integer :: ndims_
    ! HDF5 variables
    integer(hid_t) :: dataspace_id, group_id
    integer(hsize_t), allocatable :: dims_(:)
    integer :: ierr 
    character*100 :: errmsg
    ! if the dataset represents complex data, create 2D array
    if (fcomplex) then
      ndims_=ndims+1
      allocate(dims_(ndims_))
      dims_(1)=2
      dims_(2:)=dims
    else
      ndims_=ndims
      allocate(dims_(ndims_))
      dims_(:)=dims
    end if
    ! create the dataspace
    call h5screate_simple_f(ndims_,dims_,dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_setup): h5screate_simple_f returned ",I6)')ierr
      goto 10
    endif    
    ! open group
    call h5gopen_f(file_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_setup): h5gopen_f returned ",I6)')ierr
      goto 10
    endif    
    ! create the dataset
    call h5dcreate_f(group_id,trim(dname),H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(phdf5_setup): h5dcreate_f returned ",i6)')ierr
      goto 10
    endif    
    ! close dataset, group, dataspace
    call h5gclose_f(group_id,ierr)
    call h5sclose_f(dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(phdf5_setup): closing returned ",i6)')ierr
      goto 10
    endif    
    deallocate(dims_)
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims_
    write(*,'("  fname : ",A)')trim(dname)
    write(*,'("  path  : ",A)')trim(path)
    stop
  end subroutine
  
!-------------------------------------------------------------------------------
  subroutine phdf5_setup_read(ndims,dims,fcomplex,dname,path,file_id,dataset_id)
    use hdf5
    implicit none
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    logical, intent(in) :: fcomplex
    character(*), intent(in) :: dname, path
    integer(hid_t), intent(in) :: file_id
    integer(hid_t), intent(out) :: dataset_id
    !local variables
    integer :: ndims_
    ! HDF5 variables
    integer(hid_t) ::group_id
    integer(hsize_t), allocatable :: dims_(:)
    integer :: ierr 
    character*100 :: errmsg
    ! if the dataset represents complex data, create 2D array
    if (fcomplex) then
      ndims_=ndims+1
      allocate(dims_(ndims_))
      dims_(1)=2
      dims_(2:)=dims
    else
      ndims_=ndims
      allocate(dims_(ndims_))
      dims_(:)=dims(:)
    end if
    ! open group
    call h5gopen_f(file_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_setup): h5gopen_f returned ",I6)')ierr
      goto 10
    endif    
    ! open the dataset
    call h5dopen_f(group_id,trim(dname),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(phdf5_setup): h5dopen_f returned ",i6)')ierr
      goto 10
    endif    
    ! close group
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("error(phdf5_setup): closing returned ",i6)')ierr
      goto 10
    endif    
    deallocate(dims_)
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims_
    write(*,'("  fname : ",A)')trim(dname)
    write(*,'("  path  : ",A)')trim(path)
    stop
  end subroutine
!-------------------------------------------------------------------------------
  subroutine phdf5_cleanup(dataset_id)
    use hdf5
    implicit none
    integer(hid_t), intent(in) :: dataset_id
    ! local variable
    integer :: ierr

    call h5dclose_f(dataset_id,ierr)
  end subroutine

!-----------------------------------------------------------------------------
  subroutine phdf5_write_d(val,fparallel,dims,dimsg,offset,dataset_id)
    use hdf5
    implicit none
    real(8), intent(in) :: val
    logical, intent(in) :: fparallel
    integer, dimension(:), intent(in) :: dims, dimsg, offset
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer(hsize_t), allocatable, dimension(:) :: dims_, dimsg_, offset_
    integer :: ndims_
    ! get number of dimensions & allocate hdf5 size arrays
    ndims_=size(dims)
    allocate(dims_(ndims_),dimsg_(ndims_),offset_(ndims_))
    ! set local arrays
    dims_(:)=dims(:)
    dimsg_(:)=dimsg(:)
    offset_(:)=offset(:)
    ! write to hdf5
    call phdf5_write_array_d(val,fparallel,ndims_,dims_,dimsg_,offset_,dataset_id)
    !deallocate arrays
    deallocate(dims_,dimsg_,offset_)
  end subroutine

!-----------------------------------------------------------------------------
  subroutine phdf5_write_z(val,fparallel,dims,dimsg,offset,dataset_id)
    use hdf5
    implicit none
    complex(8), intent(in) :: val
    logical, intent(in) :: fparallel
    integer, dimension(:), intent(in) :: dims, dimsg, offset
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer(hsize_t), allocatable, dimension(:) :: dims_, dimsg_, offset_
    integer :: ndims_
    ! get number of dimensions & allocate hdf5 size arrays
    ndims_=size(dims)+1
    allocate(dims_(ndims_),dimsg_(ndims_),offset_(ndims_))
    ! set local arrays
    dims_(1)=2
    dims_(2:)=dims(:)
    dimsg_(1)=2
    dimsg_(2:)=dimsg(:)
    offset_(1)=0
    offset_(2:)=offset(:)
    ! write to hdf5
    call phdf5_write_array_d(val,fparallel,ndims_,dims_,dimsg_,offset_,dataset_id)
    !deallocate arrays
    deallocate(dims_,dimsg_,offset_)
  end subroutine

!-----------------------------------------------------------------------------
  subroutine phdf5_read_d(val,fparallel,dims,dimsg,offset,dataset_id)
    use hdf5
    implicit none
    real(8), intent(out) :: val
    logical, intent(in) :: fparallel
    integer, dimension(:), intent(in) :: dims, dimsg, offset
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer(hsize_t), allocatable, dimension(:) :: dims_, dimsg_, offset_
    integer :: ndims_
    ! get number of dimensions & allocate hdf5 size arrays
    ndims_=size(dims)
    allocate(dims_(ndims_),dimsg_(ndims_),offset_(ndims_))
    ! set local arrays
    dims_(:)=dims(:)
    dimsg_(:)=dimsg(:)
    offset_(:)=offset(:)
    ! write to hdf5
    call phdf5_read_array_d(val,fparallel,ndims_,dims_,dimsg_,offset_,dataset_id)
    !deallocate arrays
    deallocate(dims_,dimsg_,offset_)
  end subroutine

!-----------------------------------------------------------------------------
  subroutine phdf5_read_i(val,fparallel,dims,dimsg,offset,dataset_id)
    use hdf5
    implicit none
    integer(4), intent(out) :: val
    logical, intent(in) :: fparallel
    integer(4), dimension(:), intent(in) :: dims, dimsg, offset
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer(hsize_t), allocatable, dimension(:) :: dims_, dimsg_, offset_
    integer(4) :: ndims_
    ! get number of dimensions & allocate hdf5 size arrays
    ndims_=size(dims)
    allocate(dims_(ndims_),dimsg_(ndims_),offset_(ndims_))
    ! set local arrays
    dims_(:)=dims(:)
    dimsg_(:)=dimsg(:)
    offset_(:)=offset(:)
    ! write to hdf5
    call phdf5_read_array_i(val,fparallel,ndims_,dims_,dimsg_,offset_,dataset_id)
    !deallocate arrays
    deallocate(dims_,dimsg_,offset_)
  end subroutine
!-----------------------------------------------------------------------------
  subroutine phdf5_read_z(val,fparallel,dims,dimsg,offset,dataset_id)
    use hdf5
    implicit none
    complex(8), intent(out) :: val
    logical, intent(in) :: fparallel
    integer, dimension(:), intent(in) :: dims, dimsg, offset
    integer(hid_t), intent(in) :: dataset_id
    ! local variables
    integer(hsize_t), allocatable, dimension(:) :: dims_, dimsg_, offset_
    integer :: ndims_
    ! get number of dimensions & allocate hdf5 size arrays
    ndims_=size(dims)+1
    allocate(dims_(ndims_),dimsg_(ndims_),offset_(ndims_))
    ! set local arrays
    dims_(1)=2
    dims_(2:)=dims(:)
    dimsg_(1)=2
    dimsg_(2:)=dimsg(:)
    offset_(1)=0
    offset_(2:)=offset(:)
    ! write to hdf5
    call phdf5_read_array_d(val,fparallel,ndims_,dims_,dimsg_,offset_,dataset_id)
    !deallocate arrays
    deallocate(dims_,dimsg_,offset_)
  end subroutine

!-------------------------------------------------------------------------------
    subroutine phdf5_get_dims(file_id,path,datasetname,dims)
        use hdf5

        implicit none
        integer(hid_t) :: file_id
        character(*), intent(in) :: path, datasetname
        integer, intent(out) :: dims(:)
        
        integer(hid_t) :: dset_id, dspace_id, group_id
        integer(hsize_t), allocatable, dimension(:) :: dims_, maxdims_       
        integer :: ierr, i
        character*100 :: errmsg
        
        ! allocate dims_ and maxdims_
        allocate(dims_(size(dims)))
        allocate(maxdims_(size(dims)))
        
        ! open group
        call h5gopen_f(file_id,trim(path),group_id,ierr)
        if (ierr.ne.0) then
          write(errmsg,'("Error(phdf5_get_dims): h5gopen_f returned ",I6)')ierr
          goto 10
        endif
        ! open dataset
        call h5dopen_f(group_id,datasetname,dset_id,ierr)
        if (ierr.ne.0) then
          write(errmsg,'("Error(phdf5_get_dims): h5dopen_f returned ",I6)')ierr
          goto 10
        endif
        ! get dims from dataset
        call h5dget_space_f(dset_id,dspace_id,ierr)
        ! ierr = dataspace rank on success, so only give error message when ierr=-1
        if (ierr.ne.0) then
          write(errmsg,'("Error(phdf5_get_dims): h5dget_space_f returned ",I6)')ierr
          goto 10
        endif
        call h5sget_simple_extent_dims_f(dspace_id,dims_,maxdims_,ierr)
        if (ierr.eq.-1) then
          write(errmsg,'("Error(phdf5_get_dims): h5sget_simple_extent_dims_f returned ",I6)')ierr
          goto 10
        endif
        
        ! write dimension sizes to integer
        do i=1,size(dims)
          dims(i)=dims_(i)
        end do
        call h5sclose_f(dspace_id, ierr)
        if (ierr.ne.0) then
          write(errmsg,'("Error(hdf5_get_dims): h5sclose_f returned ",I6)')ierr
          goto 10
        endif
        call h5dclose_f(dset_id, ierr)
        if (ierr.ne.0) then
          write(errmsg,'("Error(hdf5_get_dims): h5dclose_f returned ",I6)')ierr
          goto 10
        endif
        call h5gclose_f(group_id, ierr)
        if (ierr.ne.0) then
          write(errmsg,'("Error(hdf5_get_dims): h5gclose_f returned ",I6)')ierr
          goto 10
        endif
        
        deallocate(dims_,maxdims_)
        return
        10 continue
        write(*,'(A)')trim(errmsg)
        write(*,'("  file_id : ",I4)')file_id
        write(*,'("  path  : ",A)')trim(path)
        write(*,'("  datasetname    : ",A)')trim(datasetname)
        write(*,'("  dims  : ",10I4)')dims
        stop
        
    end subroutine
end module
!-----------------------------------------------------------------------------
subroutine phdf5_write_array_d(a,fparallel,ndims,dims,dimsg,offset,dataset_id)
  use hdf5
  implicit none
  real(8), intent(in) :: a(*)
  logical, intent(in) :: fparallel
  integer, intent(in) :: ndims
  integer(hsize_t), intent(in) :: dims(ndims), dimsg(ndims), offset(ndims)
  integer(hid_t), intent(in) :: dataset_id
  ! local variables
  integer :: ierr
  integer(hid_t) :: memspace_id, dataspace_id, plist_id
  character*100 :: errmsg
  ! create memoryspace
  call h5screate_simple_f(ndims,dims,memspace_id,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_write_array_d): h5screate_f returned ",I6)')ierr
    goto 10
  endif
  ! select hyperslab in file
  call h5dget_space_f(dataset_id,dataspace_id,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_write_array_d): h5dget_space_f returned ",I6)')ierr
    goto 10
  endif
  call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,offset,dims,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_write_array_d): h5sselect_hyperslab_f returned ",I6)')ierr
    goto 10
  endif
  ! write into hyperslab
  ! MPI I/O independently
  if (fparallel) then
    ! create property list for individual dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5pcreate_f returned ",I6)')ierr
      goto 10
    endif
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5pset_dxpl_mpio_f returned ",I6)')ierr
      goto 10
    endif
    ! write the dataset
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,dimsg,ierr,memspace_id,dataspace_id,plist_id)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5dwrite_f returned ",I6)')ierr
      goto 10
    endif
    ! close the property list
    call h5pclose_f(plist_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5pclose_f returned ",I6)')ierr
      goto 10
    endif
  ! serial
  else
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,dimsg,ierr,memspace_id,dataspace_id)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5dwrite_f returned ",I6)')ierr
      goto 10
    endif
  end if
  ! close memory space and dataspace
  call h5sclose_f(dataspace_id,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_write_array_d): h5sclose_f returned ",I6)')ierr
    goto 10
  endif
  call h5sclose_f(memspace_id,ierr) 
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_write_array_d): h5sclose_f returned ",I6)')ierr
    goto 10
  endif
  
  return
  10 continue
  write(*,'(A)')trim(errmsg)
  write(*,'("  dataset_id : ",10I4)')dataset_id
  write(*,'("  offset  : ",10I4)')offset
  write(*,'("  dims  : ",10I4)')dims
  write(*,'("  dimsg  : ",10I4)')dimsg
  stop
end subroutine

!-----------------------------------------------------------------------------
subroutine phdf5_read_array_d(a,fparallel,ndims,dims,dimsg,offset,dataset_id)
  use hdf5
  implicit none
  real(8), intent(out) :: a(*)
  logical, intent(in) :: fparallel
  integer, intent(in) :: ndims
  integer(hsize_t), intent(in) :: dims(ndims), dimsg(ndims), offset(ndims)
  integer(hid_t), intent(in) :: dataset_id
  ! local variables
  integer :: ierr
  integer(hid_t) :: memspace_id, dataspace_id, plist_id
  character*100 :: errmsg
  ! create memoryspace
  call h5screate_simple_f(ndims,dims,memspace_id,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_read_array_d): h5screate_simple_f returned ",I6)')ierr
    goto 10
  endif
  ! select hyperslab in file
  call h5dget_space_f(dataset_id,dataspace_id,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_read_array_d): h5dget_space_f returned ",I6)')ierr
    goto 10
  endif
  call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,offset,dims,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_read_array_d): h5sselect_hyperslab_f returned ",I6)')ierr
    goto 10
  endif
  ! read from hyperslab
  ! MPI read
  if (fparallel) then
    ! create property list for individual dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5pcreate_f returned ",I6)')ierr
      goto 10
    endif
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5pset_dxpl_mpio_f returned ",I6)')ierr
      goto 10
    endif
    ! read dataset into memory
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,dimsg,ierr,memspace_id,dataspace_id,plist_id)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    ! close the property list
    call h5pclose_f(plist_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5pclose_f returned ",I6)')ierr
      goto 10
    endif
  ! serial read
  else
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,dimsg,ierr,memspace_id,dataspace_id)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_write_array_d): h5dread_f returned ",I6)')ierr
      goto 10
    endif
  end if
  ! close memory space and dataspace
  call h5sclose_f(dataspace_id,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_write_array_d): h5sclose_f returned ",I6)')ierr
    goto 10
  endif
  call h5sclose_f(memspace_id,ierr) 
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_write_array_d): h5sclose_f returned ",I6)')ierr
    goto 10
  endif
  
  return
  10 continue
  write(*,'(A)')trim(errmsg)
  write(*,'("  dataset_id : ",I4)')dataset_id
  write(*,'("  offset  : ",10I4)')offset
  write(*,'("  dims  : ",10I4)')dims
  write(*,'("  dimsg : ",10I4)')dimsg
  stop
end subroutine

!-----------------------------------------------------------------------------
subroutine phdf5_read_array_i(a,fparallel,ndims,dims,dimsg,offset,dataset_id)
  use hdf5
  implicit none
  integer(4), intent(out) :: a(*)
  logical, intent(in) :: fparallel
  integer(4), intent(in) :: ndims
  integer(hsize_t), intent(in) :: dims(ndims), dimsg(ndims), offset(ndims)
  integer(hid_t), intent(in) :: dataset_id
  ! local variables
  integer(4) :: ierr
  integer(hid_t) :: memspace_id, dataspace_id, plist_id
  character*100 :: errmsg

  ! create memoryspace
  call h5screate_simple_f(ndims,dims,memspace_id,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_read_array_i): h5screate_simple_f returned ",I6)')ierr
    goto 10
  endif
  ! select hyperslab in file
  call h5dget_space_f(dataset_id,dataspace_id,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_read_array_i): h5dget_space_f returned ",I6)')ierr
    goto 10
  endif
  call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,offset,dims,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_read_array_i): h5sselect_hyperslab_f returned ",I6)')ierr
    goto 10
  endif
  ! read from hyperslab
  ! MPI read
  if (fparallel) then
    ! create property list for individual dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_read_array_i): h5pcreate_f returned ",I6)')ierr
      goto 10
    endif
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_read_array_i): h5pset_dxpl_mpio_f returned ",I6)')ierr
      goto 10
    endif
    ! read dataset into memory
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,dimsg,ierr,memspace_id,dataspace_id,plist_id)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_read_array_i): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    ! close the property list
    call h5pclose_f(plist_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_read_array_i): h5pclose_f returned ",I6)')ierr
      goto 10
    endif
  ! serial read
  else
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,dimsg,ierr,memspace_id,dataspace_id)
    if (ierr.ne.0) then
      write(errmsg,'("Error(phdf5_read_array_i): h5dread_f returned ",I6)')ierr
      goto 10
    endif
  end if
  ! close memory space and dataspace
  call h5sclose_f(dataspace_id,ierr)
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_write_array_i): h5sclose_f returned ",I6)')ierr
    goto 10
  endif
  call h5sclose_f(memspace_id,ierr) 
  if (ierr.ne.0) then
    write(errmsg,'("Error(phdf5_write_array_i): h5sclose_f returned ",I6)')ierr
    goto 10
  endif
 
  return
  10 continue
  write(*,'(A)')trim(errmsg)
  write(*,'("  dataset_id : ",I4)')dataset_id
  write(*,'("  offset  : ",10I4)')offset
  write(*,'("  dims  : ",10I4)')dims
  write(*,'("  dimsg  : ",10I4)')dimsg
  stop
end subroutine
