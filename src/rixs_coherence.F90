program rixs_coherence
  use mod_phdf5
  use mod_mpi
  use mod_io
  use mod_matmul
  use mod_blocks
  use mod_blocks_k
  use hdf5, only: hid_t

  implicit none
  real(8) :: broad
  real(8) :: pol(3)
  integer :: nkmax, w1, nw_, lambda
  integer :: interdim(2)
  type(io) :: optical, core
  type(input) :: inputparam
  character(1024) :: fname_core, fname_optical,fname_pmat, fname_output, &
   & gname_c, gname_ic, gname_w, cw1
  integer(4) :: nblocks_, ik, ik1, ik2, koulims_comb(4)
  integer(4) :: blocks_, blocks2_
  type(block1d) :: oscstr_b, evals_b, evals2_b, t1_b
  type(block2d) :: t2_b, tprime_b
  complex(8) :: alpha, beta
  !MPI variables
  ! PHDF5 variables
  integer(hid_t) :: optical_id, output_id, core_id, energy_id, pmat_id
  integer(hid_t), allocatable :: coherent_id(:), incoherent_id(:)
  integer :: matsize_(1)
  !Specify file/dataset name
  fname_core='./core_output.h5'
  fname_optical='./optical_output.h5'
  fname_pmat='./pmat.h5'
  fname_output='./rixs.h5'
  ! initialize MPI and HDF5 interface
  call initmpi()
  call phdf5_initialize()
    
  ! open BSE and pmat files
  call phdf5_open_file(fname_optical,optical_id)
  call phdf5_open_file(fname_core,core_id)
  call phdf5_open_file(fname_pmat,pmat_id)
  !create output hdf5 files
  call phdf5_create_file(fname_output,output_id)
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
  call read_inputfile(inputparam)
  
  ! set parameters
  broad=inputparam%broad
  pol=inputparam%pol
  
  ! test whether the blocksize is possible
  interdim=shape(optical%koulims)
  nkmax=interdim(2)
  
  ! define blocks for oscillator strength
  nblocks_=inputparam%nblocks

  !-------------------------------------------------!
  !    Calculation of the oscillator strength       !
  !-------------------------------------------------!
  ! create group in output file
  
  call phdf5_create_group(output_id,'/','oscstr')
  matsize_=(/ inputparam%nstato /)
  call phdf5_setup_write(1,matsize_,.false.,'evals','/',output_id,energy_id)
  
  gname_c='coherent'
  gname_ic='incoherent'
  ! prepare datasets for coherent oscillator-strength dataset for each frequency
  nw_=size(inputparam%omega)
  if (allocated(coherent_id)) deallocate(coherent_id)
  allocate(coherent_id(nw_))
  if (inputparam%calc_incoherent) then
    if (allocated(incoherent_id)) deallocate(incoherent_id)
    allocate(incoherent_id(nw_))
  end if

  do w1=1, nw_
    write(cw1, '(I4.4)') w1
    gname_w='/oscstr/'//cw1//'/'
    ! create group for each frequency
    if (.not. phdf5_exist_group(output_id,'/oscstr/',trim(cw1))) then
      call phdf5_create_group(output_id,'/oscstr/',trim(cw1))
    end if
    matsize_=(/ inputparam%nstato/)
    call phdf5_setup_write(1,matsize_,.true.,'coherent',trim(gname_w),output_id,coherent_id(w1))
    if (inputparam%calc_incoherent) then
      call phdf5_setup_write(1,matsize_,.true.,'incoherent',trim(gname_w),output_id,incoherent_id(w1))
    end if
  end do

  !----------------------------------------------------!
  !    Write optical BSE eigenvalues E_{\lambda}       !
  !----------------------------------------------------!
  if (mpiglobal%rank .eq. 0) then
    do blocks_=1, nblocks_
      ! set up block for eigenvalues (needed only for file output)
      evals_b%nblocks=nblocks_
      evals_b%blocksize=nofblock(blocks_, inputparam%nstato, nblocks_)
      evals_b%global=inputparam%nstato
      evals_b%il=firstofblock(blocks_, inputparam%nstato, nblocks_)
      evals_b%iu=lastofblock(blocks_, inputparam%nstato, nblocks_)
      evals_b%offset=firstofblock(blocks_, inputparam%nstato, nblocks_)-1
      evals_b%id=blocks_
      
      call get_evals(evals_b,optical_id)
      call put_block1d(evals_b,energy_id)
    end do
  end if
  call phdf5_cleanup(energy_id)
  
  !----------------------------------------------------!
  !    Write coherent oscillator strength |t^{3}_c|    !
  !----------------------------------------------------!
  do w1=1, nw_
    ! loop over kpoints
    do blocks_=firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
      ! set up block of oscillator strength
      oscstr_b%nblocks=nblocks_
      oscstr_b%blocksize=nofblock(blocks_, inputparam%nstato, nblocks_)
      oscstr_b%global=inputparam%nstato
      oscstr_b%il=firstofblock(blocks_, inputparam%nstato, nblocks_)
      oscstr_b%iu=lastofblock(blocks_, inputparam%nstato, nblocks_)
      oscstr_b%offset=firstofblock(blocks_, inputparam%nstato, nblocks_)-1
      oscstr_b%id=blocks_

      ! allocate content of oscillator strength
      if (allocated(oscstr_b%zcontent)) deallocate(oscstr_b%zcontent)
      allocate(oscstr_b%zcontent(oscstr_b%blocksize))
      oscstr_b%zcontent(:)=cmplx(0.0d0, 0.0d0)
      
      ! sum over k-points
      do ik=1, nkmax
        ! generate combined koulims map
        koulims_comb(1)=optical%koulims(3,ik)
        koulims_comb(2)=optical%koulims(4,ik)
        koulims_comb(3)=core%koulims(3,ik)
        koulims_comb(4)=core%koulims(4,ik)
        
        ! generate block of tprime
        tprime_b%nblocks=nblocks_
        tprime_b%id=ik
        call generate_tprime_k(tprime_b, ik, inputparam%pol,koulims_comb, pmat_id)
        
        do blocks2_=1, nblocks_
          ! set up block for eigenvalues
          evals2_b%nblocks=nblocks_
          evals2_b%blocksize=nofblock(blocks2_, inputparam%nstatc, nblocks_)
          evals2_b%global=inputparam%nstatc
          evals2_b%il=firstofblock(blocks2_, inputparam%nstatc, nblocks_)
          evals2_b%iu=lastofblock(blocks2_, inputparam%nstatc, nblocks_)
          evals2_b%offset=firstofblock(blocks2_, inputparam%nstatc, nblocks_)-1
          evals2_b%id=blocks2_
          
          ! set up block of t(1)
          t1_b%nblocks=nblocks_
          t1_b%blocksize=nofblock(blocks2_, inputparam%nstatc, nblocks_)
          t1_b%global=inputparam%nstatc
          t1_b%il=firstofblock(blocks2_, inputparam%nstatc, nblocks_)
          t1_b%iu=lastofblock(blocks2_, inputparam%nstatc, nblocks_)
          t1_b%offset=firstofblock(blocks2_, inputparam%nstatc, nblocks_)-1
          t1_b%id=blocks2_
          
          !set up block for t(2) matrix
          t2_b%nblocks=nblocks_
          t2_b%blocksize=(/ nofblock(blocks_ , inputparam%nstato, nblocks_), &
           & nofblock(blocks2_, inputparam%nstatc, nblocks_) /)
          t2_b%global=(/ inputparam%nstato, inputparam%nstatc /)
          t2_b%il=firstofblock(blocks_, inputparam%nstato, nblocks_)
          t2_b%iu=lastofblock(blocks_, inputparam%nstato, nblocks_)
          t2_b%jl=firstofblock(blocks2_, inputparam%nstatc, nblocks_)
          t2_b%ju=lastofblock(blocks2_, inputparam%nstatc, nblocks_)
          t2_b%offset(1)=t2_b%il-1
          t2_b%offset(2)=t2_b%jl-1
          t2_b%id=(/ blocks_, blocks2_ /)
          
          ! generate block of core eigenvalues
          call get_evals(evals2_b,core_id)
          ! prepare content for t(1) and t(2)
          if (allocated(t1_b%zcontent)) deallocate(t1_b%zcontent)
          allocate(t1_b%zcontent(t1_b%blocksize))
          if (allocated(t2_b%zcontent)) deallocate(t2_b%zcontent)
          allocate(t2_b%zcontent(t2_b%blocksize(1), t2_b%blocksize(2)))
          ! generate block of t(1) and t(2)
          call gen_t1_k(t1_b, ik, core, core_id, pmat_id, inputparam)
          call gen_t2_k(t2_b, ik, tprime_b, core, optical, core_id, &
            optical_id, inputparam)
          
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

        end do ! blocks2_
      end do ! ik
      ! write oscillator strength
      call put_block1d(oscstr_b,coherent_id(w1))
    end do !blocks_
  end do ! w1
  
  !----------------------------------------------------!
  !    Write coherent oscillator strength |t^{3}_ic|   !
  !----------------------------------------------------!
 if (inputparam%calc_incoherent) then
   ! loop over frequencies
   do w1=1, nw_
      do blocks_=firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
        ! set up block of oscillator strength
        oscstr_b%nblocks=nblocks_
        oscstr_b%blocksize=nofblock(blocks_, inputparam%nstato, nblocks_)
        oscstr_b%global=inputparam%nstato
        oscstr_b%il=firstofblock(blocks_, inputparam%nstato, nblocks_)
        oscstr_b%iu=lastofblock(blocks_, inputparam%nstato, nblocks_)
        oscstr_b%offset=firstofblock(blocks_, inputparam%nstato, nblocks_)-1
        oscstr_b%id=blocks_

        ! allocate content of oscillator strength
        if (allocated(oscstr_b%zcontent)) deallocate(oscstr_b%zcontent)
        allocate(oscstr_b%zcontent(oscstr_b%blocksize))
        oscstr_b%zcontent(:)=cmplx(0.0d0, 0.0d0)
        
        ! 1st loop over k-points
        do ik1=1, nkmax
          ! generate combined koulims map
          koulims_comb(1)=optical%koulims(3,ik1)
          koulims_comb(2)=optical%koulims(4,ik1)
          koulims_comb(3)=core%koulims(3,ik1)
          koulims_comb(4)=core%koulims(4,ik1)
          ! generate block of tprime
          tprime_b%nblocks=nblocks_
          tprime_b%id=ik
          call generate_tprime_k(tprime_b, ik1, inputparam%pol,koulims_comb, pmat_id)
          ! 2nd loop over k-poins
          do ik2=1, nkmax
            ! ik1=ik2 is the coherent contribution
            if (ik1 .ne. ik2) then
              do blocks2_=1, nblocks_
                ! set up block for eigenvalues
                evals2_b%nblocks=nblocks_
                evals2_b%blocksize=nofblock(blocks2_, inputparam%nstatc, nblocks_)
                evals2_b%global=inputparam%nstatc
                evals2_b%il=firstofblock(blocks2_, inputparam%nstatc, nblocks_)
                evals2_b%iu=lastofblock(blocks2_, inputparam%nstatc, nblocks_)
                evals2_b%offset=firstofblock(blocks2_, inputparam%nstatc, nblocks_)-1
                evals2_b%id=blocks2_
                
                !set up block for t(2) matrix
                t2_b%nblocks=nblocks_
                t2_b%blocksize=(/ nofblock(blocks_ , inputparam%nstato, nblocks_), &
                 & nofblock(blocks2_, inputparam%nstatc, nblocks_) /)
                t2_b%global=(/ inputparam%nstato, inputparam%nstatc /)
                t2_b%il=firstofblock(blocks_, inputparam%nstato, nblocks_)
                t2_b%iu=lastofblock(blocks_, inputparam%nstato, nblocks_)
                t2_b%jl=firstofblock(blocks2_, inputparam%nstatc, nblocks_)
                t2_b%ju=lastofblock(blocks2_, inputparam%nstatc, nblocks_)
                t2_b%offset(1)=t2_b%il-1
                t2_b%offset(2)=t2_b%jl-1
                t2_b%id=(/ blocks_, blocks2_ /)
                
                ! set up block of t(1)
                t1_b%nblocks=nblocks_
                t1_b%blocksize=nofblock(blocks2_, inputparam%nstatc, nblocks_)
                t1_b%global=inputparam%nstatc
                t1_b%il=firstofblock(blocks2_, inputparam%nstatc, nblocks_)
                t1_b%iu=lastofblock(blocks2_, inputparam%nstatc, nblocks_)
                t1_b%offset=firstofblock(blocks2_, inputparam%nstatc, nblocks_)-1
                t1_b%id=blocks2_
                
                
                ! generate block of core eigenvalues
                call get_evals(evals2_b,core_id)
                ! prepare content for t(1) and t(2)
                if (allocated(t1_b%zcontent)) deallocate(t1_b%zcontent)
                allocate(t1_b%zcontent(t1_b%blocksize))
                if (allocated(t2_b%zcontent)) deallocate(t2_b%zcontent)
                allocate(t2_b%zcontent(t2_b%blocksize(1), t2_b%blocksize(2)))
                ! generate block of t(1) and t(2)
                call gen_t2_k(t2_b, ik1, tprime_b, core, optical, core_id, &
                  optical_id, inputparam)
                call gen_t1_k(t1_b, ik2, core, core_id, pmat_id, inputparam)
                ! adjust t(1) by multiplication with frequency-dependent prefactor
                do lambda=1, t1_b%blocksize
                  t1_b%zcontent(lambda)=(-1.0d0/(evals2_b%dcontent(lambda)*27.211d0-inputparam%omega(w1) &
                    &+cmplx(0.0d0,inputparam%broad)))*t1_b%zcontent(lambda)
                end do
                ! generate block of oscstr
                alpha=1.0d0
                beta=1.0d0
                call zgemm('N', 'N', t2_b%blocksize(1), 1, t2_b%blocksize(2), &
                 & alpha, t2_b%zcontent, t2_b%blocksize(1), t1_b%zcontent, &
                 & t1_b%blocksize, beta, oscstr_b%zcontent, oscstr_b%blocksize)
              end do ! blocks2_
            end if 
          end do ! ik2
        end do ! ik1
        call put_block1d(oscstr_b,incoherent_id(w1))
      end do ! nw_
    end do ! blocks_
  end if 
  do w1=1, nw_
    call phdf5_cleanup(coherent_id(w1))
    if (inputparam%calc_incoherent) call phdf5_cleanup(incoherent_id(w1))
  end do
  ! close HDF5 files
  call phdf5_close_file(optical_id)
  call phdf5_close_file(core_id)
  call phdf5_close_file(output_id)
  call phdf5_close_file(pmat_id)
  !close HDF5 files
  call phdf5_finalize()
  call finitmpi()

end program rixs_coherence

