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
  integer :: nkmax, w1, nw_, lambda
  integer :: interdim(2)
  type(io) :: optical
  type(input) :: inputparam
  character(1024) :: fname_core, fname_optical,fname_pmat, fname_output, &
   & fname_inter, gname_eval, gname2, gname3, cw1 
  integer(4) :: nblocks_, blsz_,nk_, k, k2, blsz_k
  integer(4) :: nu_optical, no_optical, global_optical
  type(block1d) :: oscstr_b
  type(block1d) :: t1_b, evals_b, evals2_b
  type(block2d) :: t2_b
  complex(8) :: alpha, beta
  real(8) ::  test
  !MPI variables
  ! PHDF5 variables
  integer(hid_t) :: optical_id, output_id, core_id, inter_id
  integer(hid_t) :: energy_id, t1_id, t2_id
  integer(hid_t), allocatable :: oscstr_id(:)
  integer :: matsize_(1), matsize2_(2)
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
  call phdf5_open_file(fname_core,core_id)
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
  pol=inputparam%pol
  
  ! test whether the blocksize is possible
  interdim=shape(optical%koulims)
  nkmax=interdim(2)
  nu_optical=optical%koulims(2,1)-optical%koulims(1,1)+1
  no_optical=optical%koulims(4,1)-optical%koulims(3,1)+1
  test=float(nkmax)/float(inputparam%nblocks)
  
  if ((float(floor(test)) .ne. test) .and. (mpiglobal%rank .eq. 0)) then
    print *, 'Blocksize', inputparam%nblocks, 'not compatible with ', nkmax, 'k-points'
    stop
  end if 
  ! define blocks for oscillator strength
  nk_=nkmax/inputparam%nblocks
  nblocks_=inputparam%nblocks
  global_optical=nu_optical*no_optical*nkmax

  !-------------------------------------------------!
  !    Calculation of the oscillator strength       !
  !-------------------------------------------------!
  ! create group in output file
  call phdf5_create_group(output_id,'/','oscstr')
  ! open datasets for write of energies
  matsize_=(/ inputparam%nstato /)
  call phdf5_setup_write(1,matsize_,.false.,'evals','/',output_id,energy_id)
  ! prepare datasets for oscillator strength, dataset for each frequency
  nw_=size(inputparam%omega)
  if (allocated(oscstr_id)) deallocate(oscstr_id)
  allocate(oscstr_id(nw_))
  do w1=1, nw_
    write(cw1, '(I4.4)') w1
    call phdf5_setup_write(1,matsize_,.true.,trim(adjustl(cw1)),'/oscstr/',output_id,oscstr_id(w1))
  end do
  ! open datasets for read of t(1) and t(2)
  matsize_=(/ inputparam%nstatc/)
  call phdf5_setup_read(1,matsize_,.true.,'t(1)','/',inter_id,t1_id)
  matsize2_=(/ inputparam%nstato, inputparam%nstatc/)
  call phdf5_setup_read(2,matsize2_,.true.,'t(2)','/',inter_id,t2_id)
  ! loop over blocks
  do k=firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
    ! set up block for eigenvalues (needed only for file output)
    evals_b%nblocks=nblocks_
    evals_b%blocksize=nofblock(k, inputparam%nstato, nblocks_)
    evals_b%global=inputparam%nstato
    evals_b%nk=nk_
    evals_b%il=firstofblock(k, inputparam%nstato, nblocks_)
    evals_b%iu=lastofblock(k, inputparam%nstato, nblocks_)
    evals_b%offset=firstofblock(k, inputparam%nstato, nblocks_)-1
    evals_b%id=k
    ! generate block of eigenvalues
    if (.not. inputparam%ip_o) then
      call get_evals_block(evals_b,optical_id)
    else
      call get_evalsIP_block(evals_b,optical_id)
    end if
    call put_block1d(evals_b,energy_id)
    do w1=1, nw_
      ! set up block of oscillator strength
      oscstr_b%nblocks=nblocks_
      oscstr_b%blocksize=nofblock(k, inputparam%nstato, nblocks_)
      oscstr_b%global=inputparam%nstato
      oscstr_b%nk=nk_
      oscstr_b%il=firstofblock(k, inputparam%nstato, nblocks_)
      oscstr_b%iu=lastofblock(k, inputparam%nstato, nblocks_)
      oscstr_b%offset=firstofblock(k, inputparam%nstato, nblocks_)-1
      oscstr_b%id=k
      ! allocate content of oscillator strength
      if (allocated(oscstr_b%zcontent)) deallocate(oscstr_b%zcontent)
      allocate(oscstr_b%zcontent(oscstr_b%blocksize))
      oscstr_b%zcontent(:)=cmplx(0.0d0, 0.0d0)
      do k2=1, nblocks_
        ! set up block for eigenvalues
        evals2_b%nblocks=nblocks_
        evals2_b%blocksize=nofblock(k2, inputparam%nstatc, nblocks_)
        evals2_b%global=inputparam%nstato
        evals2_b%nk=nk_
        evals2_b%il=firstofblock(k2, inputparam%nstatc, nblocks_)
        evals2_b%iu=lastofblock(k2, inputparam%nstatc, nblocks_)
        evals2_b%offset=firstofblock(k2, inputparam%nstatc, nblocks_)-1
        evals2_b%id=k2
        
        ! set up block of t(1)
        t1_b%nblocks=nblocks_
        t1_b%blocksize=nofblock(k2, inputparam%nstatc, nblocks_)
        t1_b%global=inputparam%nstatc
        t1_b%il=firstofblock(k2, inputparam%nstatc, nblocks_)
        t1_b%iu=lastofblock(k2, inputparam%nstatc, nblocks_)
        t1_b%nk=nk_
        t1_b%offset=firstofblock(k2, inputparam%nstatc, nblocks_)-1
        t1_b%id=k2
        
        !set up block for t(2) matrix
        t2_b%nblocks=nblocks_
        t2_b%blocksize=(/ nofblock(k, inputparam%nstato, nblocks_), nofblock(k2, inputparam%nstatc, nblocks_) /)
        t2_b%global=(/ inputparam%nstato, inputparam%nstatc /)
        t2_b%nk=nk_
        t2_b%il=firstofblock(k, inputparam%nstato, nblocks_)
        t2_b%iu=lastofblock(k, inputparam%nstato, nblocks_)
        t2_b%jl=firstofblock(k2, inputparam%nstatc, nblocks_)
        t2_b%ju=lastofblock(k2, inputparam%nstatc, nblocks_)
        t2_b%offset(1)=t2_b%il-1
        t2_b%offset(2)=t2_b%jl-1
        t2_b%id=(/ k, k2 /)
        
        ! generate block of core eigenvalues
        if (.not. inputparam%ip_c) then
          call get_evals_block(evals2_b,core_id)
        else
          call get_evalsIP_block(evals2_b,core_id)
        end if
        ! prepare content for t(1) and t(2)
        if (allocated(t1_b%zcontent)) deallocate(t1_b%zcontent)
        allocate(t1_b%zcontent(t1_b%blocksize))
        if (allocated(t2_b%zcontent)) deallocate(t2_b%zcontent)
        allocate(t2_b%zcontent(t2_b%blocksize(1), t2_b%blocksize(2)))
        ! generate block of t(1) and t(2)
        call get_block1d(t1_b,t1_id)
        call get_block2d(t2_b,t2_id)
        ! adjust t(1) by multiplication with frequency-dependent prefactor
        do lambda=1, t1_b%blocksize
          t1_b%zcontent(lambda)=(-1.0d0/(evals2_b%dcontent(lambda)*27.211d0-inputparam%omega(w1) &
            &+cmplx(0.0d0,inputparam%broad)))*t1_b%zcontent(lambda)
        end do
        ! generate block of oscstr
        alpha=1.0d0
        beta=1.0d0
        call zgemm('N', 'N', t2_b%blocksize(1), 1, t2_b%blocksize(2), alpha, t2_b%zcontent, t2_b%blocksize(1), t1_b%zcontent, &
          & t1_b%blocksize, beta, oscstr_b%zcontent, oscstr_b%blocksize)
      end do ! k2
      ! write oscillator strength
      call put_block1d(oscstr_b,oscstr_id(w1))
    end do !w1
  end do ! k
  call phdf5_cleanup(energy_id)
  call phdf5_cleanup(t1_id)
  call phdf5_cleanup(t2_id)
  do w1=1, nw_
    call phdf5_cleanup(oscstr_id(w1))
  end do
  ! close HDF5 files
  call phdf5_close_file(optical_id)
  call phdf5_close_file(inter_id)
  call phdf5_close_file(output_id)
  !close HDF5 files
  call phdf5_finalize()
  call finitmpi()

end program

