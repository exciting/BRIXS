!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: mod_io
!
!> @author
!> Christian Vorwerk, Humboldt Universität zu Berlin.
!> Elias Richter, Humboldt Universität zu Berlin.
!
! DESCRIPTION:
!> This module provides subroutines to read the information from the HDF5 files
!> that store the information of the valence and core BSE calculations, as well 
!> as the information from the input.cfg.
!
! REVISION HISTORY:
! 09 07 2020 - Added documentation
! 24 01 2025 - 2 Pol treatment
! 31 07 2026 - Use non-equilibrium occupations
!------------------------------------------------------------------------------
module mod_io
  implicit none
  
  type :: io
    integer(4), allocatable :: smap(:,:), ismap(:,:,:), koulims(:,:), ensortidx(:)
    integer(4) :: hamsize, lu,uu,lo,uo, nk0, nkmax, nu, no, global, globalk
    real(8), allocatable :: evals(:), occupation_factors(:)
    complex(8), allocatable :: eigvecs(:,:)
  end type io
  type :: input
    real(8), allocatable :: omega(:)
    real(8) :: broad
    real(8) :: pol_in(3), pol_out(3)
    integer :: nblocks, nstato, nstatc
    logical :: ip_c, ip_o, calc_incoherent, non_equilibrium
    character(1024) :: occupation_factors_core_file
    character(1024) :: occupation_factors_optical_file
  end type
  

  public set_param
  public get_koulims
  public get_smap
  public get_ismap
  public get_occupation_factors
  public get_transition_range
    
  contains
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
  !> Elias Richter, Humboldt Universität zu Berlin.
  !
  ! DESCRIPTION: 
  !> @brief
  !> Read the square root of the transition occupation difference from an
  !> exciting transitions HDF5 file. The data must use the same transition-space
  !> ordering as the BSE eigenvectors. If an smap dataset is present, that ordering 
  !> is verified against object%smap.
  !
  !> @param[inout] object   
  !> @param[in] file_id      
  !---------------------------------------------------------------------------  
  subroutine get_occupation_factors(object,file_id)
    use hdf5, only: hid_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mod_phdf5, only: phdf5_get_dims, phdf5_setup_read, &
     &                   phdf5_read, phdf5_cleanup, phdf5_exist_group
    implicit none
    type(io), intent(inout) :: object
    integer(hid_t), intent(in) :: file_id
    integer(4) :: dims(1), offset_(1), smap_dims(2), smap_offset(2), i
    integer(4), allocatable :: occupation_smap(:,:)
    integer(hid_t) :: dataset_id
    character(len=*), parameter :: path='/IQMT_000001/transitions'
    character(len=*), parameter :: dsetname='occupation_factors'
    character(len=*), parameter :: smap_dsetname='smap'

    offset_=(/ 0 /)
    call phdf5_get_dims(file_id,path,dsetname,dims)

    if (dims(1) .ne. object%hamsize) then
      write(*,'(A,I0,A,I0)') 'Error(get_occupation_factors): dataset has ', &
        & dims(1), ' entries, but the BSE transition space has ', object%hamsize
      write(*,'(A)') 'The occupation factors and BSE eigenvectors must come from matching transition selections.'
      error stop
    end if

    if (.not. phdf5_exist_group(file_id,path,smap_dsetname)) then
      write(*,'(A)') 'Warning(get_occupation_factors): optional dataset /IQMT_000001/transitions/smap is missing.'
      write(*,'(A)') 'Skipping transition-map validation; occupation factors are assumed to already ' // &
        & 'match the BSE transition ordering.'
    else
      call phdf5_get_dims(file_id,path,smap_dsetname,smap_dims)
      if (any(smap_dims .ne. shape(object%smap))) then
        write(*,'(A,2(I0,1X),A,2(I0,1X))') 'Error(get_occupation_factors): occupation smap shape ', &
          & smap_dims, 'does not match BSE smap shape ', shape(object%smap)
        error stop
      end if

      allocate(occupation_smap(smap_dims(1),smap_dims(2)))
      smap_offset=(/ 0, 0 /)
      call phdf5_setup_read(2,smap_dims,.false.,smap_dsetname,path,file_id,dataset_id)
      call phdf5_read(occupation_smap(1,1),smap_dims,smap_dims,smap_offset,dataset_id)
      call phdf5_cleanup(dataset_id)

      do i=1,object%hamsize
        if (any(occupation_smap(:,i) .ne. object%smap(:,i))) then
          write(*,'(A,I0)') 'Error(get_occupation_factors): transition-map mismatch at index ',i
          write(*,'(A,3(I0,1X))') 'Occupation file (c, v, k): ',occupation_smap(:,i)
          write(*,'(A,3(I0,1X))') 'BSE eigenvectors (c, v, k): ',object%smap(:,i)
          error stop
        end if
      end do
      deallocate(occupation_smap)
    end if

    if (allocated(object%occupation_factors)) deallocate(object%occupation_factors)
    allocate(object%occupation_factors(dims(1)))
    call phdf5_setup_read(1,dims,.false.,dsetname,path,file_id,dataset_id)
    call phdf5_read(object%occupation_factors(1),dims,dims,offset_,dataset_id)
    call phdf5_cleanup(dataset_id)

    if (any(.not. ieee_is_finite(object%occupation_factors)) .or. &
      & any(object%occupation_factors < 0.0d0)) then
      write(*,'(A)') 'Error(get_occupation_factors): occupation_factors must be finite and non-negative.'
      error stop
    end if
  end subroutine get_occupation_factors

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
      object%ismap=0
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
  !> Elias Richter, Humboldt Universität zu Berlin.
  !
  ! DESCRIPTION: 
  !> @brief
  !> Return the contiguous transition-space range belonging to a range of
  !> relative k-point indices. exciting stores smap ordered by k point.
  !    
  !---------------------------------------------------------------------------  
  subroutine get_transition_range(object,kl,ku,il,iu)
    implicit none
    type(io), intent(in) :: object
    integer(4), intent(in) :: kl, ku
    integer(4), intent(out) :: il, iu
    integer(4) :: i, ik

    il=0
    iu=-1
    do i=1,object%hamsize
      ik=object%smap(3,i)-object%nk0+1
      if ((ik >= kl) .and. (ik <= ku)) then
        if (il == 0) il=i
        iu=i
      end if
    end do

    if (il == 0) then
      write(*,'(A,I0,A,I0)') 'Error(get_transition_range): no transitions for k-point range ',kl,' to ',ku
      error stop
    end if
  end subroutine get_transition_range
  
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
  ! 31 07 2026 - Support non-equilibrium occupations
  !
  !> @param[inout] object
  !---------------------------------------------------------------------------
  subroutine set_param(object)
    implicit none
    type(io), intent(inout) :: object
    !local variables
    integer(4), dimension(2) :: dim_koulims, dim_smap
    integer(4) :: nkmax_

    if ((allocated(object%koulims)) .and. (allocated(object%smap))) then
      !get shapes
      dim_koulims=shape(object%koulims)
      dim_smap=shape(object%smap)
      !set some parameters for convenience
      object%lu=minval(object%koulims(1,:)) ! lowest conduction band
      object%uu=maxval(object%koulims(2,:)) ! highest conduction band
      object%lo=minval(object%koulims(3,:)) ! lowest valence band
      object%uo=maxval(object%koulims(4,:)) ! highest valence band
      object%nu=object%uu-object%lu+1 ! number of conduction bands
      object%no=object%uo-object%lo+1 ! number of valence bands
      object%nk0=object%smap(3,1)     ! index of first k-point
      object%nkmax=dim_koulims(2)     ! Number of k-points
      object%hamsize=dim_smap(2)      ! Size of BSE Hamiltonian
      object%globalk=object%no*object%nu
      object%global=object%hamsize
    else
      print *, 'koulims and smap have to be obtained from file before set_param can be called!'
    end if
  end subroutine

  !---------------------------------------------------------------------------  
  !> @author
  !> Christian Vorwerk, Humboldt Universität zu Berlin.
  !> Elias Richter, Humboldt Universität zu Berlin.
  !
  ! DESCRIPTION:
  !> @brief
  !> Reads input.cfg file and stores all information in input-object
  !
  ! REVISION HISTORY:
  ! 09 07 2020 - Added documentation 
  ! 24 01 2025 - 2 Pol treatment
  ! 31 07 2026 - Use non-equilibrium occupations
  !
  !> @param[out] object   
  !---------------------------------------------------------------------------  
  subroutine read_inputfile(object)
    use mod_mpi, only: mpiglobal, ierr
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
    real(8) :: pol_in_(3)
    real(8) :: pol_out_(3)
    integer :: nblocks_, nstato_, nstatc_
    logical :: oscstr_, vecA_
    logical :: ip_c_, ip_o_, calc_incoherent_
    logical :: non_equilibrium_
    character(1024) :: occupation_factors_core_file_
    character(1024) :: occupation_factors_optical_file_


    ! only root reads the input file
#ifdef MPI
    if (mpiglobal%rank .eq. 0) then
#endif
      !define fields and set defaults
      call CFG_add(my_cfg, 'omega', (/1.0_dp, 2.0_dp/), 'Core Frequencies', dynamic_size=.true.)
      call CFG_add(my_cfg, 'pol_in', (/1.0_dp, 0.0_dp, 0.0_dp/), 'Light Polarization incoming')
      call CFG_add(my_cfg, 'broad', 0.5_dp, 'Core Broadening')
      call CFG_add(my_cfg, 'nblocks', 1, 'Number of Blocks')
      call CFG_add(my_cfg, 'eigstates_optical', 1, 'Number of eigenstates in optical BSE calculation')
      call CFG_add(my_cfg, 'eigstates_core', 1, 'Number of eigenstates in core BSE calculation')
      call CFG_add(my_cfg, 'ip_core', .false., 'IPA for core BSE calculation')
      call CFG_add(my_cfg, 'ip_optical', .false., 'IPA for optical BSE calculation')
      call CFG_add(my_cfg, '_calc_incoherent_', .false., 'Calculate the incoherent contribution')
      call CFG_add(my_cfg, 'pol_out', (/1.0_dp, 0.0_dp, 0.0_dp/), 'Light Polarization outgoing')
      call CFG_add(my_cfg, 'non_equilibrium', .false., 'Use transition occupation factors')
      call CFG_add(my_cfg, 'occupation_factors_core_file', 'occupations_core.h5', &
        & 'HDF5 file containing core transition occupation factors')
      call CFG_add(my_cfg, 'occupation_factors_optical_file', 'occupations_optical.h5', &
        & 'HDF5 file containing optical transition occupation factors')
      ! read input file
      call CFG_read_file(my_cfg, 'input.cfg')
      ! get size and values of core frequencies
      call CFG_get_size(my_cfg, 'omega', omegasize_)
      if (allocated(omega_)) deallocate(omega_)
      allocate(omega_(omegasize_))
      call CFG_get(my_cfg,'omega', omega_)
      ! get core broadening
      call CFG_get(my_cfg,'broad',broad_)
      ! get light polarization incoming
      call CFG_get(my_cfg,'pol_in',pol_in_)
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
      ! determine whether the incoherent contribution is calculated (should be avoided) 
      call CFG_get(my_cfg,'_calc_incoherent_', calc_incoherent_)
      ! get light polarization outgoing
      call CFG_get(my_cfg,'pol_out',pol_out_)
      call CFG_get(my_cfg,'non_equilibrium',non_equilibrium_)
      call CFG_get(my_cfg,'occupation_factors_core_file',occupation_factors_core_file_)
      call CFG_get(my_cfg,'occupation_factors_optical_file',occupation_factors_optical_file_)
#ifdef MPI
    end if
    ! broadcast input parameters to everybody
    call mpi_bcast(omegasize_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if (.not. allocated(omega_)) allocate(omega_(omegasize_))
    call mpi_bcast(omega_,omegasize_,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(broad_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(pol_in_,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nblocks_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nstato_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nstatc_,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ip_c_,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ip_o_,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(calc_incoherent_,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(pol_out_,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(non_equilibrium_,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(occupation_factors_core_file_,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(occupation_factors_optical_file_,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif
    
    ! get input parameters from read
    object%broad=broad_
    object%nblocks=nblocks_
    object%nstato=nstato_
    object%nstatc=nstatc_
    object%pol_in=pol_in_(:)
    object%ip_c=ip_c_
    object%ip_o=ip_o_
    object%calc_incoherent=calc_incoherent_
    object%pol_out=pol_out_(:)
    object%non_equilibrium=non_equilibrium_
    object%occupation_factors_core_file=trim(occupation_factors_core_file_)
    object%occupation_factors_optical_file=trim(occupation_factors_optical_file_)
    ! calculate frequency ranges
    if (allocated(object%omega)) deallocate(object%omega) 
    allocate(object%omega(omegasize_))
    object%omega(:)=omega_(:)
  end subroutine read_inputfile  
end module
