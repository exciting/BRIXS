!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: mod_io
!
!> @author
!> Christian Vorwerk, Humboldt Universität zu Berlin.
!
! DESCRIPTION: 
!> This module provides subroutines to read the information from the HDF5 files
!> that store the information of the valence and core BSE calculations, as well 
!> as the information from the input.cfg.
!
! REVISION HISTORY:
! 09 07 2020 - Added documentation
!------------------------------------------------------------------------------
module mod_io
  implicit none
  
  type :: io
    integer(4), allocatable :: smap(:,:), ismap(:,:,:), koulims(:,:), ensortidx(:)
    integer(4) :: hamsize, lu,uu,lo,uo, nk0, nkmax, nu, no
    real(8), allocatable :: evals(:)
    complex(8), allocatable :: eigvecs(:,:)
  end type io
  type :: input
    real(8), allocatable :: omega(:)
    real(8) :: broad
    real(8) :: pol(3)
    integer :: nblocks, nstato, nstatc
    logical :: ip_c, ip_o
  end type
  

  public set_param
  public get_koulims
  public get_smap
  public get_ismap
    
  contains
  ! Methodenbereich
  !-----------------------------------------------------------------------------

  !---------------------------------------------------------------------------  
  !> @author 
  !> Christian Vorwerk, Humboldt Universität zu Berlin.
  !
  ! DESCRIPTION: 
  !> Brief description of routine. 
  !> @brief
  !> Read valence and conduction band limits for each \f$ \mathbf{k} \f$-point
  !> from HDF5 file and store in io-object
  !
  ! REVISION HISTORY:
  ! 09 07 2020 - Added documentation 
  !
  !> @param[inout] object   
  !> @param[in] file_id      
  !---------------------------------------------------------------------------  
  subroutine get_koulims(object,file_id)
    use hdf5, only: hid_t
    use mod_phdf5, only: phdf5_get_dims, phdf5_setup_read, &
     &                   phdf5_read, phdf5_cleanup
    implicit none
    type(io), intent(inout) :: object
    integer(hid_t), intent(in) :: file_id
    !local variables
    integer(4) :: dims(2), offset_(2)
    integer(hid_t) :: dataset_id
    character(len=1024) :: path, dsetname
    !set fake offset
    offset_=(/ 0, 0 /)
    !get sizes of koulims
    path=trim(adjustl('eigvec-singlet-TDA-BAR-full/0001/parameters'))
    dsetname=trim(adjustl('koulims'))
    ! get dimensions of dataset
    call phdf5_get_dims(file_id,path,dsetname,dims)
    !allocate output
    if (allocated(object%koulims)) deallocate(object%koulims)
    allocate(object%koulims(dims(1),dims(2)))
    !open dataset
    call phdf5_setup_read(2,dims,.false.,dsetname,path,file_id,dataset_id)
    !get data
    call phdf5_read(object%koulims(1,1),dims,dims,offset_,dataset_id)
    ! close dataset
    call phdf5_cleanup(dataset_id)
  end subroutine
  
  !---------------------------------------------------------------------------  
  !> @author 
  !> Christian Vorwerk, Humboldt Universität zu Berlin.
  !
  ! DESCRIPTION: 
  !> Brief description of routine. 
  !> @brief
  !> Read mapping between transition space and valence, conduction, and k-point 
  !> index from HDF5 file and store it in io-object
  !
  ! REVISION HISTORY:
  ! 09 07 2020 - Added documentation 
  !
  !> @param[inout] object   
  !> @param[in] file_id      
  !---------------------------------------------------------------------------  
  subroutine get_smap(object,file_id)
    use hdf5, only: hid_t
    use mod_phdf5, only: phdf5_get_dims, phdf5_setup_read, &
     &                   phdf5_read, phdf5_cleanup
    implicit none
    type(io), intent(inout) :: object
    integer(hid_t), intent(in) :: file_id
    !local variables
    integer(4) :: dims(2), offset_(2)
    integer(hid_t) :: dataset_id
    character(len=1024) :: path, dsetname
    ! set fake offset
    offset_=(/ 0, 0/)
    !get sizes of koulims
    path='eigvec-singlet-TDA-BAR-full/0001/parameters'
    dsetname='smap'
    ! get dimensions of dataset
    call phdf5_get_dims(file_id,path,dsetname,dims)
    !allocate output
    if (allocated(object%smap)) deallocate(object%smap)
    allocate(object%smap(dims(1),dims(2)))
    ! open dataset
    call phdf5_setup_read(2,dims,.false.,dsetname,path,file_id,dataset_id)
    ! get data
    call phdf5_read(object%smap(1,1),dims,dims,offset_,dataset_id)
    ! close dataset
    call phdf5_cleanup(dataset_id)
  end subroutine 

  !---------------------------------------------------------------------------  
  !> @author 
  !> Christian Vorwerk, Humboldt Universität zu Berlin.
  !
  ! DESCRIPTION: 
  !> @brief
  !> Read sorting vector from HDF5 file and store it in the io-object. The 
  !> vector stores the indices that sort the IP energy differences with 
  !> increasing energy
  !
  ! REVISION HISTORY:
  ! 09 07 2020 - Added documentation 
  !
  !> @param[inout] object   
  !> @param[in] file_id      
  !---------------------------------------------------------------------------  
  subroutine get_ensortidx(object,file_id)
    use hdf5, only: hid_t
    use mod_phdf5, only: phdf5_get_dims, phdf5_setup_read, &
     &                   phdf5_read, phdf5_cleanup
    implicit none
    type(io), intent(inout) :: object
    integer(hid_t), intent(in) :: file_id
    !local variables
    integer(4) :: dims(1), offset_(1)
    integer(hid_t) :: dataset_id
    character(len=1024) :: path, dsetname
    ! set fake offset
    offset_=(/ 0/)
    !get sizes of koulims
    path='eigvec-singlet-TDA-BAR-full/0001/parameters'
    dsetname='ensortidx'
    call phdf5_get_dims(file_id,path,dsetname,dims)
    !allocate output
    if (allocated(object%ensortidx)) deallocate(object%ensortidx)
    allocate(object%ensortidx(dims(1)))
    ! open dataset
    call phdf5_setup_read(1,dims,.false.,dsetname,path,file_id,dataset_id)
    ! get data
    call phdf5_read(object%ensortidx(1),dims,dims,offset_,dataset_id)
    ! close dataset
    call phdf5_cleanup(dataset_id)
  end subroutine 

  !---------------------------------------------------------------------------  
  !> @author 
  !> Christian Vorwerk, Humboldt Universität zu Berlin.
  !
  ! DESCRIPTION: 
  !> @brief
  !> Generates the inverse smap, i.e. the map that yields for each combination 
  !> of valence band index, conduction band index, and k-point index the 
  !> corresponding transition space index. The map is stored in the io-object.
  !
  ! REVISION HISTORY:
  ! 09 07 2020 - Added documentation 
  !
  !> @param[inout] object   
  !---------------------------------------------------------------------------  
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
  
  !---------------------------------------------------------------------------  
  !> @author 
  !> Christian Vorwerk, Humboldt Universität zu Berlin.
  !
  ! DESCRIPTION: 
  !> @brief
  !> Obtains several objects from object%koulims and object%smap and stores them
  !> for convenience. 
  !
  ! REVISION HISTORY:
  ! 09 07 2020 - Added documentation 
  !
  !> @param[inout] object   
  !---------------------------------------------------------------------------  
  subroutine set_param(object)
    implicit none
    type(io), intent(inout) :: object
    !local variables
    integer(4), dimension(2) :: dim_koulims, dim_smap
    if ((allocated(object%koulims)) .and. (allocated(object%smap))) then
      !get shapes
      dim_koulims=shape(object%koulims)
      dim_smap=shape(object%smap)
      !set some parameters for convenience
      object%lu=object%koulims(1,1)   ! lowest conduction band
      object%uu=object%koulims(2,1)   ! highest conduction band
      object%lo=object%koulims(3,1)   ! lowest valence band
      object%uo=object%koulims(4,1)   ! highest valence band
      object%nu=object%uu-object%lu+1 ! number of conduction bands
      object%no=object%uo-object%lo+1 ! number of valence bands
      object%nk0=object%smap(3,1)     ! index of first k-point
      object%nkmax=dim_koulims(2)     ! Number of k-points
      object%hamsize=dim_smap(2)      ! Size of BSE Hamiltonian
    else
      print *, 'koulims and smap have to be obtained from file before set_param can be called!'
    end if
  end subroutine

  !---------------------------------------------------------------------------  
  !> @author 
  !> Christian Vorwerk, Humboldt Universität zu Berlin.
  !
  ! DESCRIPTION: 
  !> @brief
  !> Reads input.cfg file and stores all information in input-object 
  !
  ! REVISION HISTORY:
  ! 09 07 2020 - Added documentation 
  !
  !> @param[out] object   
  !---------------------------------------------------------------------------  
  subroutine read_inputfile(object)
    use modmpi, only: mpiglobal, ierr
#ifdef MPI
    use mpi
#endif
    use m_config
    implicit none
    type(input), intent(out) :: object
    ! local variables
    integer, parameter :: dp=kind(0.0d0)
    type(CFG_t) :: my_cfg
    integer :: omegasize_
    real(8), allocatable :: omega_(:)
    
    integer :: line, ios, w, pos
    character(256) :: buffer,label
    integer, parameter :: fh = 15
    real(8) :: inter(3), inter2(3)
    real(8) :: broad_, broad2_
    real(8) :: pol_(3)
    integer :: nblocks_, nstato_, nstatc_
    logical :: oscstr_, vecA_
    logical :: ip_c_, ip_o_


    ! only root reads the input file
#ifdef MPI
    if (mpiglobal%rank .eq. 0) then
#endif
      !define fields and set defaults
      call CFG_add(my_cfg, 'omega', (/1.0_dp, 2.0_dp/), 'Core Frequencies', dynamic_size=.true.)
      call CFG_add(my_cfg, 'pol', (/1.0_dp, 0.0_dp, 0.0_dp/), 'Light Polarization')
      call CFG_add(my_cfg, 'broad', 0.5_dp, 'Core Broadening')
      call CFG_add(my_cfg, 'nblocks', 1, 'Number of Blocks')
      call CFG_add(my_cfg, 'eigstates_optical', 1, 'Number of eigenstates in optical BSE calculation')
      call CFG_add(my_cfg, 'eigstates_core', 1, 'Number of eigenstates in core BSE calculation')
      call CFG_add(my_cfg, 'ip_core', .false., 'IPA for core BSE calculation')
      call CFG_add(my_cfg, 'ip_optical', .false., 'IPA for optical BSE calculation')
      ! read input file
      call CFG_read_file(my_cfg, 'input.cfg')
      ! get size and values of core frequencies
      call CFG_get_size(my_cfg, 'omega', omegasize_)
      if (allocated(omega_)) deallocate(omega_)
      allocate(omega_(omegasize_))
      call CFG_get(my_cfg,'omega', omega_)
      ! get core broadening
      call CFG_get(my_cfg,'broad',broad_)
      ! get light polarization
      call CFG_get(my_cfg,'pol',pol_)
      ! get number of blocks
      call CFG_get(my_cfg,'nblocks', nblocks_)
      ! get number of optical eigenstates
      call CFG_get(my_cfg,'eigstates_optical', nstato_)
      ! get number of core eigenstates
      call CFG_get(my_cfg,'eigstates_core', nstatc_)
      ! get approximation of core calculation
      call CFG_get(my_cfg,'ip_core', ip_c_)
      ! get approximation of optical calculation
      call CFG_get(my_cfg,'ip_optical', ip_o_)
#ifdef MPI
    end if
    ! broadcast input parameters to everybody
    call mpi_bcast(omegasize_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if (.not. allocated(omega_)) allocate(omega_(omegasize_))
    call mpi_bcast(omega_,omegasize_,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(broad_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(pol_,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nblocks_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nstato_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nstatc_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ip_c_,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ip_o_,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif
    
    ! get input parameters from read
    object%broad=broad_
    object%nblocks=nblocks_
    object%nstato=nstato_
    object%nstatc=nstatc_
    object%pol=pol_(:)
    object%ip_c=ip_c_
    object%ip_o=ip_o_
    ! calculate frequency ranges
    if (allocated(object%omega)) deallocate(object%omega) 
    allocate(object%omega(omegasize_))
    object%omega(:)=omega_(:)
  end subroutine read_inputfile  
end module

