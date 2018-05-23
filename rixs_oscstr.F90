program rixs_oscstr
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
  type(io) :: optical
  type(input) :: inputparam
  character(1024) :: fname_core, fname_optical,fname_pmat, fname_output, &
   & fname_inter, gname, gname2, gname3, ik 
  integer(4) :: nblocks_, blsz_,nk_, k, k2, blsz_k
  type(block1d) :: oscstr_b
  type(block1d) :: vecA_b, evals_b
  type(block2d) :: evecs_b
  complex(8) :: alpha, beta
  real(8) ::  test
  !MPI variables
  ! PHDF5 variables
  integer(hid_t) :: optical_id, output_id, inter_id
  integer(hid_t) :: energy_id, oscstr_id, vecA_id
  integer :: matsize_(1), reduced_(1), matsize2_(1)
  !Specify file/dataset name
  fname_core='./core_output.h5'
  fname_optical='./optical_output.h5'
  fname_pmat='./pmat.h5'
  fname_inter='./data.h5'
  fname_output='./rixs.h5'
  ! initialize MPI and HDF5 interface
  call initmpi()
  call phdf5_initialize()

  call phdf5_open_file(fname_optical,optical_id)
  !open intermediate data hdf5 files
  call phdf5_open_file(fname_inter,inter_id)
   !create output hdf5 files
  call phdf5_create_file(fname_output,output_id)
  !initialize io objects for core and optics
  
  call get_koulims(optical,optical_id)
  call get_smap(optical,optical_id)
  call set_param(optical)
  call get_ismap(optical)
  ! read input file
  call read_inputfile(inputparam)
  
  ! set parameters
  broad=inputparam%broad
  allocate(omega(size(inputparam%omega))) 
  omega(:)=inputparam%omega(:) 
  
  ! set polarization vector
  pol(1)=1.0d0
  pol(2)=0.0d0
  pol(3)=0.0d0
  
  ! test whether the blocksize is possible
  interdim=shape(optical%koulims)
  nkmax=interdim(2)
  nu=optical%koulims(2,1)-optical%koulims(1,1)+1
  no=optical%koulims(4,1)-optical%koulims(3,1)+1
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
  matsize_=(/ nblocks_*blsz_ /)

  ! open log file
  if (rank .eq. 0) then
    open(unit=8,file="log_oscstr.txt",action="write",status="replace")
  end if
  !-------------------------------------------------!
  !    Calculation of the oscillator strength       !
  !-------------------------------------------------!
  ! create group in intermediate file
  call phdf5_create_group(output_id,'/','oscstr')
  call phdf5_create_group(output_id,'/','eval')
  gname3=trim(adjustl('/eval/'))
  do w1=1, size(omega)
    if (rank .eq. 0) write(8, *) 'Frequency ', w1, 'of ', size(omega) 
    ! set up groups in output file
    write(ik, '(I4.4)') w1
    call phdf5_create_group(output_id, '/oscstr/',ik)
    gname=trim(adjustl('/oscstr/'//trim(ik)//'/'))
    gname2=trim(adjustl('/A/'//trim(ik)//'/'))
    ! open datasets for write of energies & oscstr
    matsize2_(1)=inputparam%nstato
    if (w1==1) then 
      call phdf5_setup_write(1,matsize2_,.false.,'data',gname3,output_id,energy_id)
    endif
    call phdf5_setup_write(1,matsize2_,.true.,'data',gname,output_id,oscstr_id)
    ! open dataset to read vecA
    call phdf5_setup_read(1,matsize_,.true.,'A',gname2,inter_id,vecA_id)
    ! loop over blocks
    do k=firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
      if (w1==1) then
        ! set up block for eigenvalues (needed only for file output)
        reduced_(1)=nofblock(k, inputparam%nstato, nblocks_)
        print*, 'nstato=', inputparam%nstato
        evals_b%nblocks=nblocks_
        evals_b%blocksize=reduced_(1)
        evals_b%global=inputparam%nstato
        evals_b%il=firstofblock(k, inputparam%nstato, nblocks_)
        evals_b%iu=lastofblock(k, inputparam%nstato, nblocks_)
        evals_b%offset=firstofblock(k, inputparam%nstato, nblocks_)-1
        evals_b%id=k
        ! generate block of eigenvalues
        call get_evals_block(evals_b,optical_id)
        print *, 'il=', evals_b%il
        print *, 'iu=', evals_b%iu
        print *, 'offset(evals_b)=', evals_b%offset
        call put_block1d(evals_b,energy_id)
      end if
      reduced_(1)=nofblock(k, inputparam%nstato, nblocks_)
      ! set up block of oscillator strength
      oscstr_b%nblocks=nblocks_
      oscstr_b%blocksize=reduced_(1)
      oscstr_b%nk=nk_
      oscstr_b%il=firstofblock(k, inputparam%nstato, nblocks_)
      oscstr_b%iu=lastofblock(k, inputparam%nstato, nblocks_)
      oscstr_b%offset=firstofblock(k, inputparam%nstato, nblocks_)-1
      oscstr_b%id=k
      ! generate block of oscillator strength
      ! allocate content
      if (allocated(oscstr_b%zcontent)) deallocate(oscstr_b%zcontent)
      allocate(oscstr_b%zcontent(reduced_(1)))
      if (allocated(vecA_b%zcontent)) deallocate(vecA_b%zcontent)
      allocate(vecA_b%zcontent(blsz_))
      oscstr_b%zcontent(:)=0.0d0 
      do k2=1, nblocks_
        ! set up block for eigenvectors
        evecs_b%nblocks=nblocks_
        evecs_b%blocksize=(/blsz_, reduced_(1)/)
        evecs_b%global=(/matsize_(1), inputparam%nstato /)
        evecs_b%il=(k2-1)*blsz_+1
        evecs_b%iu=k2*blsz_
        evecs_b%jl=firstofblock(k, inputparam%nstato, nblocks_)
        evecs_b%ju=lastofblock(k, inputparam%nstato, nblocks_)
        evecs_b%offset(1)=(k2-1)*blsz_
        evecs_b%offset(2)=firstofblock(k, inputparam%nstato, nblocks_)-1
        evecs_b%id(1)=k2
        evecs_b%id(2)=k
        !set up block for A matrix
        vecA_b%nblocks=nblocks_
        vecA_b%blocksize=blsz_
        vecA_b%global=matsize_(1)
        vecA_b%nk=nk_
        vecA_b%il=(k2-1)*blsz_+1
        vecA_b%iu=k2*blsz_
        vecA_b%kl=(k2-1)*nk_+1
        vecA_b%ku=k2*nk_
        vecA_b%offset=(k2-1)*blsz_
        vecA_b%id=k2
       
        ! generate block of eigenvectors
        call get_eigvecs2D_b(evecs_b,optical_id)
        call get_block1d(vecA_b,vecA_id)
        ! generate block of oscstr
        alpha=1.0d0
        beta=1.0d0
        call zgemm('C','N',reduced_(1),1,blsz_,alpha,evecs_b%zcontent,blsz_,vecA_b%zcontent,blsz_,beta, &
          &         oscstr_b%zcontent(:),reduced_(1))
      end do ! k2
      ! write oscillator strength
      call put_block1d(oscstr_b,oscstr_id)
    end do ! k
    if (w1==1) then
      call phdf5_cleanup(energy_id)
    end if
    call phdf5_cleanup(vecA_id)
    call phdf5_cleanup(oscstr_id)
  end do ! w
  ! close HDF5 files
  call phdf5_close_file(optical_id)
  call phdf5_close_file(inter_id)
  call phdf5_close_file(output_id)
  ! close log file
  if (rank .eq. 0) close(8)
  !close HDF5 files
  call phdf5_finalize()
  call finitmpi()

end program

