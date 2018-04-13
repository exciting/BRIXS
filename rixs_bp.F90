program rixs_b2
  use mod_phdf5
  use modmpi
  use mod_io
  use mod_rixs
  use mod_matmul
  use mod_blocks

  implicit none
  real(8) :: broad, broad2
  real(8), allocatable :: omega(:), omega2(:)
  real(8) :: pol(3)
  integer :: nkmax, nu, no, w1
  integer :: interdim(2), lambda
  type(io) :: optical, core
  type(input) :: inputparam
  integer(4), allocatable :: koulims_comb(:,:)
  character(1024) :: fname_core, fname_optical,fname_pmat, fname_output, &
   & fname_inter, gname, gname2, gname3, ik, 
  integer(4) :: nblocks_, blsz_,nk_, k, k2, blsz_k
  type(block1d) :: oscstr_b
  type(block1d) :: vecA_b, evals_b
  type(block2d) :: evecs_b
  complex(8) :: alpha, beta
  real(8) :: start, finish, test
  !MPI variables
  ! PHDF5 variables
  integer(hid_t) :: core_id, optical_id, pmat_id, output_id, inter_id
  integer(hid_t) :: dataset_id, energy_id, oscstr_id, vecA_id
  integer :: matsize_
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
  call phdf5_create_file(fname_output,.True.,output_id,mpiglobal%comm)

  ! initialize io objects for core and optics
  call get_koulims(optical,fname_optical)
  call get_koulims(core,fname_core)
  call get_smap(optical,fname_optical)
  call get_smap(core,fname_core)
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
  ! create group in intermediate file
  call phdf5_create_group(inter_id,'/','A')
  do w1=1, size(omega)
    write(7, *) 'Frequency ', w1, 'of ', size(omega) 
    write(ik, '(I4.4)') w1
    gname=trim(adjustl(ik))
    ! create group for each frequency
    if (.not. phdf5_exist_group(inter_id,'/A/',gname)) then
      call phdf5_create_group(inter_id,'/A/',gname)
    end if
    gname2=trim(adjustl('/A/'//gname//'/'))
    ! set up the dataset
    matsize_=nblocks_*blsz_
    datasetname='A'
    call phdf5_setup_write(1,matsize_,.true.,trim(datasetname),gname2,file_id,dataset_id)
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
  
  !-------------------------------------------------!
  !    Calculation of the oscillator strength       !
  !-------------------------------------------------!
  ! create group in intermediate file
  call phdf5_create_group(output_id,'/','oscstr')
  call phdf5_create_group(output_id,'/','eval')
  gname3=trim(adjustl('/eval/'))
  do w1=1, size(omega)
    ! set up groups in output file
    write(ik, '(I4.4)') w1
    call hdf5_create_group(output_id, '/oscstr/',ik)
    gname=trim(adjustl('/oscstr/'//trim(ik)//'/'))
    gname2=trim(adjustl('/A/'//trim(ik)//'/'))
    ! set oscillator strength to zero 
    oscstr_b%zcontent(:)=0.0d0
    ! open datasets for write of energies & oscstr
    if (w1==1) then 
      call phdf5_setup_write(1,matsize_,.false.,'data',gname3,output_id,energy_id)
    endif
    call phdf5_setup_write(1,matsize_,.true.,'data',gname,output_id,oscstr_id)
    ! open dataset to read vecA
    call phdf5_setup_read(1,matsize_,.true.,'A',gname2,inter_id,vecA_id)
    ! loop over blocks
    do k=firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
      if (w1==1) then
        ! set up block for eigenvalues (needed only for file output)
        evals_b%nblocks=nblocks_
        evals_b%blocksize=blsz_
        evals_b%nk=nk_
        evals_b%il=(k-1)*blsz_+1
        evals_b%iu=k*blsz_
        evals_b%kl=(k-1)*nk_+1
        evals_b%ku=k*nk_
        evals_b%offset=(k-1)*blsz_
        evals_b%id=k
        ! generate block of eigenvalues
        call get_evals_block(evals_b,fname_optical)
        call put_block1d(evals_b,.true.,energy_id)
      end if
      ! set up block of oscillator strength
      oscstr_b%nblocks=nblocks_
      oscstr_b%blocksize=blsz_
      oscstr_b%nk=nk_
      oscstr_b%il=(k-1)*blsz_+1
      oscstr_b%iu=k*blsz_
      oscstr_b%kl=(k-1)*nk_+1
      oscstr_b%ku=k*nk_
      oscstr_b%offset=(k-1)*blsz_
      oscstr_b%id=k
      ! generate block of oscillator strength
      ! allocate content
      if (allocated(oscstr_b%zcontent)) deallocate(oscstr_b%zcontent)
      allocate(oscstr_b%zcontent(blsz_))
      oscstr_b%zcontent(:)=0.0d0 
      do k2=1, nblocks_
        ! set up block for eigenvectors
        evecs_b%nblocks=nblocks_
        evecs_b%blocksize=blsz_
        evecs_b%nk=nk_
        evecs_b%il=(k2-1)*blsz_+1
        evecs_b%iu=k2*blsz_
        evecs_b%k1l=(k2-1)*nk_+1
        evecs_b%k1u=k2*nk_
        evecs_b%jl=(k-1)*blsz_+1
        evecs_b%ju=k*blsz_
        evecs_b%k2l=(k-1)*nk_+1
        evecs_b%k2u=k*nk_
        evecs_b%offset(1)=(k2-1)*blsz_
        evecs_b%offset(2)=(k-1)*blsz_
        evecs_b%id(1)=k2
        evecs_b%id(2)=k
        !set up block for A matrix
        vecA_b%nblocks=nblocks_
        vecA_b%blocksize=blsz_
        vecA_b%nk=nk_
        vecA_b%il=(k2-1)*blsz_+1
        vecA_b%iu=k2*blsz_
        vecA_b%kl=(k2-1)*nk_+1
        vecA_b%ku=k2*nk_
        vecA_b%offset=(k2-1)*blsz_
        vecA_b%id=k2
       
        ! generate block of eigenvectors
        call get_eigvecs2D_b(evecs_b,optical_id)
        call get_block1d(vecA_b,.true.,vecA_id)
        ! generate block of oscstr
        alpha=1.0d0
        beta=1.0d0
        call zgemm('C','N',blsz_,1,blsz_,alpha,evecs_b%zcontent,blsz_,vecA_b%zcontent,blsz_,beta,oscstr_b%zcontent(:),blsz_)
      end do ! k2
      ! write oscillator strength
      call put_block1d(oscstr_b,.true.,oscstr_id)
    end do ! k
    if (w1==1) then
      call phdf5_cleanup(energy_id)
    end if
    call phdf5_cleanup(vecA_id)
    call phdf5_cleanup(oscstr_id)
  end do ! w
  !close HDF5 files
  call phdf5_close_file(fname_core)
  call phdf5_close_file(fname_optical)
  call phdf5_close_file(fname_pmat)
  call phdf5_close_file(fname_inter)
  call phdf5_close_file(fname_output)
  call phdf5_finalize()
  call finitmpi()

end program
