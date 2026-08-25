!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! PROGRAM: rixs_coherence
!
!> @author
!> Christian Vorwerk, Humboldt Universität zu Berlin.
!> Elias Richter, Humboldt Universität zu Berlin.
!
! DESCRIPTION:
!> This program computes the coherent and incoherent contributions to the
!> RIXS oscillator strength from the matrices \f$ t^{(1)} \f$ and
!> \f$ t^{(2)} \f$ stored in data.h5, and writes the resulting \f$ t^{(3)} \f$ 
!> to rixs.h5.
!
! REVISION HISTORY:
! 15 09 2025 - 2 Pol treatment
! 31 07 2026 - Use non-non_equilibrium occupations
!
!------------------------------------------------------------------------------
program rixs_coherence
  use mod_phdf5
  use mod_mpi
  use mod_io
  use mod_matmul
  use mod_blocks
  use mod_blocks_k
  use mod_constants, only: hartree_to_ev
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
  type(block1d) :: oscstr_b, evalsv_b, evalsc_b, evals2_b, t1_b
  type(block2d) :: t2_b, tprime_out_b
  complex(8) :: alpha, beta
  !MPI variables
  ! PHDF5 variables
  integer(hid_t) :: optical_id, output_id, core_id, energyv_id, energyc_id, pmat_id
  integer(hid_t) :: occupations_core_id, occupations_optical_id
  integer(hid_t) :: omega_id
  integer(hid_t), allocatable :: coherent_id(:), incoherent_id(:)
  integer :: matsize_(1), offset_(1)
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
  if (inputparam%non_equilibrium) then
    call phdf5_open_file(trim(inputparam%occupation_factors_core_file),occupations_core_id)
    call phdf5_open_file(trim(inputparam%occupation_factors_optical_file),occupations_optical_id)
    call get_occupation_factors(core,occupations_core_id)
    call get_occupation_factors(optical,occupations_optical_id)
    call phdf5_close_file(occupations_core_id)
    call phdf5_close_file(occupations_optical_id)
  end if
  
  ! set parameters
  broad=inputparam%broad
  pol=inputparam%pol_in
  
  ! test whether the blocksize is possible
  interdim=shape(optical%koulims)
  nkmax=interdim(2)
  
  ! define blocks for oscillator strength
  nblocks_=inputparam%nblocks

  !-------------------------------------------------!
  !    Calculation of the oscillator strength       !
  !-------------------------------------------------!
  ! create groups in output file
  
  call phdf5_create_group(output_id,'/','oscstr')
  call phdf5_create_group(output_id,'/','omega')
  ! open datasets for write of valence excitation energies
  matsize_=(/ inputparam%nstato /)
  call phdf5_setup_write(1,matsize_,.false.,'vevals','/',output_id,energyv_id)
  ! open datasets for write of core excitation energies
  matsize_=(/ inputparam%nstatc /)
  call phdf5_setup_write(1,matsize_,.false.,'cevals','/',output_id,energyc_id)
  gname_c='coherent'
  gname_ic='incoherent'
  ! store the incident frequencies in the same order as the oscstr groups
  nw_=size(inputparam%omega)
  matsize_=(/ nw_ /)
  offset_=(/ 0 /)
  call phdf5_setup_write(1,matsize_,.false.,'values','/omega/',output_id,omega_id)
  if (mpiglobal%rank .eq. 0) then
    call phdf5_write(inputparam%omega(1),matsize_,matsize_,offset_,omega_id)
  end if
  call barrier(mpiglobal)
  call phdf5_cleanup(omega_id)
  ! prepare datasets for coherent oscillator-strength dataset for each frequency
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
  !    Write BSE eigenvalues E_{\lambda}       !
  !----------------------------------------------------!
  !--------------------------------------------!
  ! write valence eigenvalues (vevals)         !
  !--------------------------------------------!
  if (mpiglobal%rank .eq. 0) then
    do blocks_=1, nblocks_
      evalsv_b = make_block1d(blocks_, inputparam%nstato, nblocks_)

      call get_evals(evalsv_b, optical_id)
      call put_block1d(evalsv_b, energyv_id)
    end do
  end if
  call phdf5_cleanup(energyv_id)

  !--------------------------------------------
  ! write core eigenvalues (cevals)
  !--------------------------------------------
  if (mpiglobal%rank .eq. 0) then
    do blocks_=1, nblocks_
      evalsc_b = make_block1d(blocks_, inputparam%nstatc, nblocks_)

      call get_evals(evalsc_b, core_id)
      call put_block1d(evalsc_b, energyc_id)
    end do
  end if
  call phdf5_cleanup(energyc_id)  
  
  !----------------------------------------------------!
  !    Write coherent oscillator strength |t^{3}_c|    !
  !----------------------------------------------------!
  do w1=1, nw_
    ! loop over kpoints
    do blocks_=firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
      ! set up block of oscillator strength
      oscstr_b = make_block1d(blocks_, inputparam%nstato, nblocks_)

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
        
        ! generate block of tprime with pol_out
        tprime_out_b%nblocks=nblocks_
        tprime_out_b%id=ik
        call generate_tprime_k(tprime_out_b, ik, inputparam%pol_out, koulims_comb, pmat_id)
        
        do blocks2_=1, nblocks_
          ! set up block for eigenvalues
          evals2_b = make_block1d(blocks2_, inputparam%nstatc, nblocks_)

          ! set up block of t(1)
          t1_b = make_block1d(blocks2_, inputparam%nstatc, nblocks_)

          !set up block for t(2) matrix
          t2_b = make_block2d(blocks_, inputparam%nstato, blocks2_, inputparam%nstatc, nblocks_)
          
          ! generate block of core eigenvalues
          call get_evals(evals2_b,core_id)
          ! prepare content for t(1) and t(2)
          if (allocated(t1_b%zcontent)) deallocate(t1_b%zcontent)
          allocate(t1_b%zcontent(t1_b%blocksize))
          if (allocated(t2_b%zcontent)) deallocate(t2_b%zcontent)
          allocate(t2_b%zcontent(t2_b%blocksize(1), t2_b%blocksize(2)))
          ! generate block of t(1) and t(2)
          call gen_t1_k(t1_b, ik, core, core_id, pmat_id, inputparam)
          call gen_t2_k(t2_b, ik, tprime_out_b, core, optical, core_id, &
            optical_id, inputparam)
          
          ! adjust t(1) by multiplication with frequency-dependent prefactor
          do lambda=1, t1_b%blocksize
            t1_b%zcontent(lambda)=(-1.0d0/(evals2_b%dcontent(lambda)*hartree_to_ev-inputparam%omega(w1) &
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
  !    Write incoherent oscillator strength |t^{3}_ic|   !
  !----------------------------------------------------!
 if (inputparam%calc_incoherent) then
   ! loop over frequencies
   do w1=1, nw_
      do blocks_=firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
        ! set up block of oscillator strength
        oscstr_b = make_block1d(blocks_, inputparam%nstato, nblocks_)

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
          ! generate block of tprime with pol_out
          tprime_out_b%nblocks=nblocks_
          tprime_out_b%id=ik
          call generate_tprime_k(tprime_out_b, ik1, inputparam%pol_out, koulims_comb, pmat_id)
          ! 2nd loop over k-poins
          do ik2=1, nkmax
            ! ik1=ik2 is the coherent contribution
            if (ik1 .ne. ik2) then
              do blocks2_=1, nblocks_
                ! set up block for eigenvalues
                evals2_b = make_block1d(blocks2_, inputparam%nstatc, nblocks_)

                !set up block for t(2) matrix
                t2_b = make_block2d(blocks_, inputparam%nstato, blocks2_, inputparam%nstatc, nblocks_)

                ! set up block of t(1)
                t1_b = make_block1d(blocks2_, inputparam%nstatc, nblocks_)
                
                
                ! generate block of core eigenvalues
                call get_evals(evals2_b,core_id)
                ! prepare content for t(1) and t(2)
                if (allocated(t1_b%zcontent)) deallocate(t1_b%zcontent)
                allocate(t1_b%zcontent(t1_b%blocksize))
                if (allocated(t2_b%zcontent)) deallocate(t2_b%zcontent)
                allocate(t2_b%zcontent(t2_b%blocksize(1), t2_b%blocksize(2)))
                ! generate block of t(1) and t(2)
                call gen_t2_k(t2_b, ik1, tprime_out_b, core, optical, core_id, &
                  optical_id, inputparam)
                call gen_t1_k(t1_b, ik2, core, core_id, pmat_id, inputparam)
                ! adjust t(1) by multiplication with frequency-dependent prefactor
                do lambda=1, t1_b%blocksize
                  t1_b%zcontent(lambda)=(-1.0d0/(evals2_b%dcontent(lambda)*hartree_to_ev-inputparam%omega(w1) &
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
