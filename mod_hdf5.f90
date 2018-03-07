
!===============================================================================
!
! Author: Anton Kozhevnikov
! Added:  10.06.2014 by DIN
! 
! Modifications:
!   Added character type, 11.06.2014 by DIN
!   Added Reading of dimensions, 5.01.2018 by Christian Vorwerk
!===============================================================================

module mod_hdf5

    interface hdf5_write
        module procedure hdf5_write_i4, &
        &                hdf5_write_d,  &
        &                hdf5_write_z,  &
        &                hdf5_write_c,  &
        &                hdf5_write_l
    end interface

    interface hdf5_read
        module procedure hdf5_read_i4, &
        &                hdf5_read_d,  &
        &                hdf5_read_z,  &
        &                hdf5_read_c,  &
        &                hdf5_read_l
    end interface

    interface hdf5_read_block
        module procedure hdf5_read_block_i4, &
        &                hdf5_read_block_d,  &
        &                hdf5_read_block_z
    end interface

    public hdf5_initialize
    public hdf5_finalize
    public hdf5_create_file
    public hdf5_create_group

    public hdf5_get_dims

contains

!-------------------------------------------------------------------------------
    subroutine hdf5_initialize
        use hdf5
        implicit none
        integer ierr
        call h5open_f(ierr)
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_finalize
        use hdf5
        implicit none
        integer ierr
        call h5close_f(ierr)
    end subroutine

!-------------------------------------------------------------------------------
    subroutine hdf5_create_file(fname)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        integer ierr
        integer(hid_t) h5_root_id
        call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_file): h5fcreate_f returned ",I6)')ierr
          goto 10
        endif
        call h5fclose_f(h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_file): h5fclose_f returned ",I6)')ierr
          goto 10
        endif
        return
        10 continue
        write(*,'("  fname: ",A)')trim(fname)
        stop
    end subroutine
    
!-------------------------------------------------------------------------------
    logical function hdf5_exist_group(fname,path,gname)
        ! Check if group with the given name exists
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: gname
        integer(hid_t):: h5_root_id, h5_group_id
        integer ierr
 
        call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5fopen_f returned ",I6)')ierr
          goto 10
        endif
        call h5gopen_f(h5_root_id,trim(path),h5_group_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5gopen_f returned ",I6)')ierr
          goto 10
        endif
        call h5lexists_f(h5_group_id,trim(gname),hdf5_exist_group,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5lexists_f returned ",I6)')ierr
          goto 10
        endif
        call h5fclose_f(h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5fclose_f returned ",I6)')ierr
          goto 10
        endif
        return
        10 continue
        write(*,'("  fname : ",A)')trim(fname)
        write(*,'("  path  : ",A)')trim(path)
        write(*,'("  gname : ",A)')trim(gname)  
        stop
    end function hdf5_exist_group

!-------------------------------------------------------------------------------    
    subroutine hdf5_create_group(fname,path,gname)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: gname
        integer(hid_t) h5_root_id,h5_group_id,h5_new_group_id
        integer ierr

        call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5fopen_f returned ",I6)')ierr
          goto 10
        endif
        call h5gopen_f(h5_root_id,trim(path),h5_group_id,ierr)
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
        call h5fclose_f(h5_root_id,ierr)
        if (ierr.ne.0) then
          write(*,'("Error(hdf5_create_group): h5fclose_f returned ",I6)')ierr
          goto 10
        endif
        return
        10 continue
        write(*,'("  fname : ",A)')trim(fname)
        write(*,'("  path  : ",A)')trim(path)
        write(*,'("  gname : ",A)')trim(gname)  
        stop
    end subroutine

!-------------------------------------------------------------------------------
    subroutine hdf5_get_dims(fname,path,sname,dims)
        use hdf5

        implicit none
        character(*), intent(in) :: fname, path, sname
        integer, intent(out) :: dims(:)
        
        integer(hid_t) root_id, dset_id, dspace_id, group_id
        integer(hsize_t), allocatable, dimension(:) :: dims_, maxdims_       
        integer ierr, i
        character*100 errmsg
        
        ! allocate dims_ and maxdims_
        allocate(dims_(size(dims)))
        allocate(maxdims_(size(dims)))

        ! open file
        call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,root_id,ierr)
        if (ierr.ne.0) then
          write(errmsg,'("Error(hdf5_get_dims): h5fopen_f returned ",I6)')ierr
          goto 10
        endif
        
        ! open group
        call h5gopen_f(root_id,trim(path),group_id,ierr)
        if (ierr.ne.0) then
          write(errmsg,'("Error(hdf5_get_dims): h5gopen_f returned ",I6)')ierr
          goto 10
        endif

        ! open dataset
        call h5dopen_f(group_id,sname,dset_id,ierr)
        if (ierr.ne.0) then
          write(errmsg,'("Error(hdf5_get_dims): h5dopen_f returned ",I6)')ierr
          goto 10
        endif
        
        ! get dims from dataset
        call h5dget_space_f(dset_id,dspace_id,ierr)
        ! ierr = dataspace rank on success, so only give error message when ierr=-1
        if (ierr.ne.0) then
          write(errmsg,'("Error(hdf5_get_dims): h5dget_space_f returned ",I6)')ierr
          goto 10
        endif
        call h5sget_simple_extent_dims_f(dspace_id,dims_,maxdims_,ierr)
        if (ierr.eq.-1) then
          write(errmsg,'("Error(hdf5_get_dims): h5sget_simple_extent_dims_f returned ",I6)')ierr
          goto 10
        endif
        
        ! write dimension sizes to integer
        do i=1,size(dims)
          dims(i)=dims_(i)
        end do
        call h5sclose_f(dspace_id, ierr)
        if (ierr.ne.0) then
          write(errmsg,'("Error(hdf5_create_group): h5sget_space_f returned ",I6)')ierr
          goto 10
        endif
        call h5dclose_f(dset_id, ierr)
        call h5fclose_f(root_id, ierr)
        
        deallocate(dims_,maxdims_)
        return
        10 continue
        write(*,'(A)')trim(errmsg)
        write(*,'("  fname : ",A)')trim(fname)
        write(*,'("  path  : ",A)')trim(path)
        write(*,'("  sname    : ",A)')trim(sname)
        write(*,'("  dims  : ",10I4)')dims
        stop
        
    end subroutine
!-------------------------------------------------------------------------------
    subroutine hdf5_write_i4(fname,path,dname,val,dims)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        integer(4), intent(in) :: val
        integer, optional, dimension(:), intent(in) :: dims
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_write_array_i4(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_write_d(fname,path,dname,val,dims)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        real(8), intent(in) :: val
        integer, optional, dimension(:), intent(in) :: dims
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_write_array_d(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_write_z(fname,path,dname,val,dims)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        complex(8), intent(in) :: val
        integer, optional, dimension(:), intent(in) :: dims
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)+1
          allocate(dims_(ndims))
          dims_(1)=2
          dims_(2:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=2
        endif
        call hdf5_write_array_d(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
    end subroutine


!-------------------------------------------------------------------------------    
    subroutine hdf5_read_i4(fname,path,dname,val,dims)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        integer(4), intent(out) :: val
        integer(4), optional, dimension(:), intent(in) :: dims
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_read_array_i4(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_read_d(fname,path,dname,val,dims)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        real(8), intent(out) :: val
        integer, optional, dimension(:), intent(in) :: dims
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_read_array_d(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_read_z(fname,path,dname,val,dims)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        complex(8), intent(out) :: val
        integer, optional, dimension(:), intent(in) :: dims
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)+1
          allocate(dims_(ndims))
          dims_(1)=2
          dims_(2:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=2
        endif
        call hdf5_read_array_d(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
    end subroutine
    
!-------------------------------------------------------------------------------    
    subroutine hdf5_read_block_i4(fname,path,dname,val,dims, offset)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        integer(4), intent(out) :: val
        integer(4), dimension(:), intent(in) :: dims
        integer(4), dimension(:), intent(in) :: offset
        ! internal variables
        integer :: ndims
        
        ndims=size(dims)
        call hdf5_read_array_block_i4(val,ndims,dims,offset,fname,path,dname)
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_read_block_d(fname,path,dname,val,dims,offset)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        real(8), intent(out) :: val
        integer, dimension(:), intent(in) :: dims
        integer, dimension(:), intent(in) :: offset
        !internal variables
        integer :: ndims
          
        ndims=size(dims)
        call hdf5_read_array_block_d(val,ndims,dims,offset,fname,path,dname)
    end subroutine

!-------------------------------------------------------------------------------    
    subroutine hdf5_read_block_z(fname,path,dname,val,dims,offset)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        complex(8), intent(out) :: val
        integer, dimension(:), intent(in) :: dims
        integer, dimension(:), intent(in) :: offset
        ! internal variables
        integer :: ndims
        integer, allocatable :: dims_(:), offset_(:)

        ndims=size(dims)+1
        allocate(dims_(ndims))
        allocate(offset_(ndims))
        dims_(1)=2
        dims_(2:ndims)=dims(:)
        offset_(1)=0
        offset_(2:ndims)=offset(:)
        
        call hdf5_read_array_d(val,ndims,dims_,offset_,fname,path,dname)
        deallocate(dims_, offset_)
    end subroutine
!-------------------------------------------------------------------------------
    subroutine hdf5_write_c(fname,path,dname,val,dims)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        character*(*), intent(in) :: val
        integer, optional, dimension(:), intent(in) :: dims
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_write_array_c(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
    end subroutine

!-------------------------------------------------------------------------------
    subroutine hdf5_write_l(fname,path,dname,val)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        logical, intent(in) :: val
        if (val) then
          call hdf5_write_array_c("true",1,1,fname,path,dname)
        else
          call hdf5_write_array_c("false",1,1,fname,path,dname)
        end if
    end subroutine
!-------------------------------------------------------------------------------
    subroutine hdf5_read_c(fname,path,dname,val,dims)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        character*(*), intent(out) :: val
        integer, optional, dimension(:), intent(in) :: dims
        integer ndims
        integer, allocatable :: dims_(:)
        if (present(dims)) then
          ndims=size(dims)
          allocate(dims_(ndims))
          dims_(1:ndims)=dims(:)
        else
          ndims=1
          allocate(dims_(ndims))
          dims_(1)=1
        endif
        call hdf5_read_array_c(val,ndims,dims_,fname,path,dname)
        deallocate(dims_)
    end subroutine    
!-------------------------------------------------------------------------------
    subroutine hdf5_read_l(fname,path,dname,val)
        use hdf5
        implicit none
        character(*), intent(in) :: fname
        character(*), intent(in) :: path
        character(*), intent(in) :: dname
        logical, intent(out) :: val
        character(100) :: val_
        call hdf5_read_array_c(val_,1,1,fname,path,dname)
        if (trim(val_) == "true") then
          val=.TRUE.
        else
          val=.FALSE.
        end if
    end subroutine    
end module

!===============================================================================


!-------------------------------------------------------------------------------
subroutine hdf5_write_array_i4(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    integer(4), intent(in) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
    integer ierr,i
    integer(hsize_t), dimension(ndims) :: h_dims
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5screate_simple_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dcreate_f(group_id,trim(nm),H5T_NATIVE_INTEGER,dataspace_id,dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5dcreate_f returned ",I6)')ierr
      goto 10
    endif 
    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5dwrite_f returned ",I6)')ierr
      goto 10
    endif 
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5sclose_f(dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5sclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_i4): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims  : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path  : ",A)')trim(path)
    write(*,'("  nm    : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_write_array_d(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    real(8), intent(in) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
    integer ierr,i
    integer(hsize_t), dimension(ndims) :: h_dims
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5screate_simple_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dcreate_f(group_id,trim(nm),H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5dcreate_f returned ",I6)')ierr
      goto 10
    endif 
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5dwrite_f returned ",I6)')ierr
      goto 10
    endif 
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5sclose_f(dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5sclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_d): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims  : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path  : ",A)')trim(path)
    write(*,'("  nm    : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_read_array_i4(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    integer(4), intent(out) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataset_id,group_id
    integer ierr,i
    integer(HSIZE_T), dimension(ndims) :: h_dims
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
    enddo

    call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,*)
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path : ",A)')trim(path)
    write(*,'("  nm : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_read_array_d(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    real(8), intent(out) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataset_id,group_id
    integer ierr,i
    integer(HSIZE_T), dimension(ndims) :: h_dims
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5dopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,*)
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path : ",A)')trim(path)
    write(*,'("  nm : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_write_array_c(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    character(*), intent(in) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm
    integer(HID_T) :: h5_root_id,dataspace_id,dataset_id,group_id
    integer(HID_T) :: filetype
    !integer(SIZE_T) :: sdim
    integer :: ierr, i
    integer(HSIZE_T), dimension(ndims) :: h_dims
    integer(HSIZE_T), dimension(ndims+1) :: data_dims
    integer(SIZE_T), dimension(ndims) :: s_len
    character*100 :: errmsg

    do i = 1, ndims
      h_dims(i) = dims(i)
    enddo
    do i = 1, product(h_dims)
      s_len(i) = len_trim(a(i))
    end do
    data_dims(1) = len(a(1))
    do i = 1, ndims
      data_dims(i+1) = h_dims(i)
    enddo
    
    call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    !------------------------------------------------------
    ! Create variable size datatype
    !call H5Tcopy_f(H5T_FORTRAN_S1,filetype,ierr)
    !sdim = maxval(s_len)
    !call H5Tset_size_f(filetype,sdim,ierr)
    CALL H5Tcopy_f(H5T_STRING,filetype,ierr)
    CALL H5Tset_strpad_f(filetype,H5T_STR_NULLPAD_F,ierr)
    !------------------------------------------------------
    call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5screate_simple_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dcreate_f(group_id,trim(nm),filetype,dataspace_id,dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5dcreate_f returned ",I6)')ierr
      goto 10
    endif
    !------------------------------------------------------
    !call h5dwrite_f(dataset_id,filetype,a,h_dims,ierr)
    CALL h5dwrite_vl_f(dataset_id,filetype,a,data_dims,s_len,ierr,dataspace_id)
    !------------------------------------------------------
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5dwrite_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5sclose_f(dataspace_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5sclose_f returned ",I6)')ierr
      goto 10
    endif
    call H5Tclose_f(filetype,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5tclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_write_array_c): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims  : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path  : ",A)')trim(path)
    write(*,'("  nm    : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_read_array_c(a,ndims,dims,fname,path,nm)
    use hdf5
    implicit none
    character(*), intent(out) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) :: h5_root_id,dataset_id,dataspace_id,group_id
    integer(HID_T) :: filetype
    integer :: ierr, i
    integer(HSIZE_T), dimension(ndims) :: h_dims
    integer(HSIZE_T), dimension(ndims+1) :: data_dims
    integer(SIZE_T), dimension(ndims) :: s_len
    integer(HSIZE_T), dimension(1:2) :: maxdims
    character*100 :: errmsg

    do i = 1, ndims
      h_dims(i) = dims(i)
    enddo
    do i = 1, product(h_dims)
      s_len(i) = len_trim(a(i))
    end do
    data_dims(1) = len(a(1))
    do i = 1, ndims
      data_dims(i+1) = h_dims(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5dopen_f returned ",I6)')ierr
      goto 10
    endif
    !------------------------------------------------------
    ! Get the datatype
    CALL H5Dget_type_f(dataset_id,filetype,ierr)
    ! Get dataspace
    CALL H5Dget_space_f(dataset_id,dataspace_id,ierr)
    CALL H5Sget_simple_extent_dims_f(dataspace_id,h_dims,maxdims,ierr)
    !call h5dread_f(dataset_id,H5T_FORTRAN_S1,a,h_dims,ierr)
    CALL h5dread_vl_f(dataset_id,filetype,a,data_dims,s_len,ierr,dataspace_id)
    !------------------------------------------------------    
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    CALL H5Tclose_f(filetype,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5tclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_c): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,*)
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims : ",10I4)')dims
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path : ",A)')trim(path)
    write(*,'("  nm : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_read_array_block_i4(a,ndims,dims,offset,fname,path,nm)
    use hdf5
    implicit none
    integer(4), intent(out) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    integer, intent(in) :: offset(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataset_id,group_id, dataspace_id
    integer ierr,i
    integer(HSIZE_T), dimension(ndims) :: h_dims, h_offset
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
      h_offset(i)=offset(i)
    enddo

    call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dget_space_f(dataset_id,dataspace_id, ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dget_space_f returned ",I6)')ierr
      goto 10
    endif
    call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,h_offset,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5sselect_hyperslab_f returned ",I6)')ierr
      goto 10
    endif
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr,file_space_id=dataspace_id)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,*)
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims : ",10I4)')dims
    write(*,'("  offset : ",10I4)')offset
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path : ",A)')trim(path)
    write(*,'("  nm : ",A)')trim(nm)
    stop
end subroutine

!-------------------------------------------------------------------------------
subroutine hdf5_read_array_block_d(a,ndims,dims,offset,fname,path,nm)
    use hdf5
    implicit none
    real(8), intent(out) :: a(*)
    integer, intent(in) :: ndims
    integer, intent(in) :: dims(ndims)
    integer, intent(in) :: offset(ndims)
    character(*), intent(in) :: fname
    character(*), intent(in) :: path
    character(*), intent(in) :: nm

    integer(hid_t) h5_root_id,dataset_id,group_id, dataspace_id
    integer ierr,i
    integer(HSIZE_T), dimension(ndims) :: h_dims, h_offset
    character*100 errmsg

    do i=1,ndims
      h_dims(i)=dims(i)
      h_offset(i)=dims(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5dopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5dget_space_f(dataset_id,dataspace_id, ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5dget_space_f returned ",I6)')ierr
      goto 10
    endif
    call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,h_offset,h_dims,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_i4): h5sselect_hyperslab_f returned ",I6)')ierr
      goto 10
    endif
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr,file_space_id=dataspace_id)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5dread_f returned ",I6)')ierr
      goto 10
    endif
    call h5dclose_f(dataset_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5dclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    return
    10 continue
    write(*,*)
    write(*,'(A)')trim(errmsg)
    write(*,'("  ndims : ",I4)')ndims
    write(*,'("  dims : ",10I4)')dims
    write(*,'("  offset : ",10I4)')offset
    write(*,'("  fname : ",A)')trim(fname)
    write(*,'("  path : ",A)')trim(path)
    write(*,'("  nm : ",A)')trim(nm)
    stop
end subroutine
