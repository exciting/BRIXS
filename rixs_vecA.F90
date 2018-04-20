program rixs_vecA
  use mod_phdf5
  use modmpi
  use mod_io
  use mod_rixs
  use mod_matmul
  use mod_blocks
  use hdf5, only: hid_t

  implicit none
  real(8) :: broad, broad2
  real(8), allocatable :: omega(:), omega2(:)
  real(8) :: pol(3)
  integer :: nkmax, nu, no, w1
  integer :: interdim(2)
  type(io) :: optical, core
  type(input) :: inputparam
  integer(4), allocatable :: koulims_comb(:,:)
  character(1024) :: fname_core, fname_optical,fname_pmat, fname_output, &
   & fname_inter, gname, gname2, ik, datasetname 
  integer(4) :: nblocks_, blsz_,nk_, k, blsz_k
  type(block1d) :: vecA_b
  real(8) ::  test
  !MPI variables
  ! PHDF5 variables
  integer(hid_t) :: core_id, optical_id, pmat_id, inter_id
  integer(hid_t) :: dataset_id
  integer :: matsize_(1)
  !Specify file/dataset name
  fname_core='./core_output.h5'
  fname_optical='./optical_output.h5'
  fname_pmat='./pmat.h5'
  fname_inter='./data.h5'
  fname_output='./rixs.h5'
  ! initialize MPI and HDF5 interface
  call initmpi()
  call phdf5_initialize()
  ! open HDF5 files
  call phdf5_open_file(fname_core,.True.,core_id,mpiglobal%comm)
  call phdf5_open_file(fname_optical,.True.,optical_id,mpiglobal%comm)
  call phdf5_open_file(fname_pmat,.True.,pmat_id,mpiglobal%comm)
  ! create HDF5 files for output and intermediate data
  call phdf5_create_file(fname_inter,.True.,inter_id,mpiglobal%comm)

  ! initialize io objects for core and optics
  call get_koulims(optical,optical_id)
  call get_koulims(core,core_id)
  call get_smap(optical,optical_id)
  call get_smap(core,core_id)
  call set_param(optical)
  call set_param(core)
  call get_ismap(optical)
  call get_ismap(core)
  ! read input file
  call read_inputfile(inputparam,'./rixs.in')
  
  ! set parameters
  broad=inputparam%broad
  broad2=inputparam%broad2
  allocate(omega(size(inputparam%omega))) 
  allocate(omega2(size(inputparam%omega2)))
  omega(:)=inputparam%omega(:) 
  omega2(:)=inputparam%omega2(:) 
  
  ! set polarization vector
  pol(1)=1.0d0
  pol(2)=0.0d0
  pol(3)=0.0d0
  
  interdim=shape(core%koulims)
  nkmax=interdim(2)
  ! create combined map for valence-core transitions
  allocate(koulims_comb(4,nkmax))
  koulims_comb(1,:)=optical%koulims(3,:)
  koulims_comb(2,:)=optical%koulims(4,:)
  koulims_comb(3,:)=core%koulims(3,:)
  koulims_comb(4,:)=core%koulims(4,:)
  
  nu=optical%koulims(2,1)-optical%koulims(1,1)+1
  no=optical%koulims(4,1)-optical%koulims(3,1)+1
  ! test whether the blocksize is possible
  test=float(nkmax)/float(inputparam%nblocks)
  
  if ((float(floor(test)) .ne. test) .and. (mpiglobal%rank .eq. 0)) then
    print *, 'Blocksize', inputparam%nblocks, 'not compatible with ', nkmax, 'k-points'
    stop
  end if 
  ! define blocks for oscillator strength
  nk_=nkmax/inputparam%nblocks
  blsz_k=nu*no
  blsz_=nu*no*nk_
  nblocks_=inputparam%nblocks
  !-------------------------------------------------!
  !    Calculate and Store blocks of the A vector   !
  !-------------------------------------------------!
  ! open log file
  if (rank .eq. 0) then
    open(unit=7,file="log_vecA.txt",action="write",status="replace")
  end if
  ! create group in intermediate file
  call phdf5_create_group(inter_id,'/','A')
  do w1=1, size(omega)
    if (rank .eq. 0) write(7, *) 'Frequency ', w1, 'of ', size(omega) 
    write(ik, '(I4.4)') w1
    gname=trim(adjustl(ik))
    ! create group for each frequency
    if (.not. phdf5_exist_group(inter_id,'/A/',gname)) then
      call phdf5_create_group(inter_id,'/A/',gname)
    end if
    gname2=trim(adjustl('/A/'//gname//'/'))
    ! set up the dataset
    matsize_=(/ nblocks_*blsz_ /)
    datasetname='A'
    call phdf5_setup_write(1,matsize_,.true.,trim(datasetname),gname2,inter_id,dataset_id)     ! loop over blocks, now distributed over MPI ranks
    ! loop over blocks, now distributed over MPI ranks
    do k=firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
      !set-up for the blocks
      vecA_b%nblocks=nblocks_
      vecA_b%blocksize=blsz_
      vecA_b%nk=nk_
      vecA_b%il=(k-1)*blsz_+1
      vecA_b%iu=k*blsz_
      vecA_b%kl=(k-1)*nk_+1
      vecA_b%ku=k*nk_
      vecA_b%offset=(k-1)*blsz_
      vecA_b%id=k
      call generate_Avector_b(vecA_b,omega(w1),broad,core,optical,pmat_id,core_id,pol)
      call put_block1d(vecA_b,.true.,dataset_id)
    end do
    call phdf5_cleanup(dataset_id)
  end do
  call phdf5_close_file(core_id)
  call phdf5_close_file(optical_id)
  call phdf5_close_file(pmat_id)
  call phdf5_close_file(inter_id)
  ! close log file
  if (rank .eq. 0) close(7)
  !close HDF5 files
  call phdf5_finalize()
  call finitmpi()
end program

