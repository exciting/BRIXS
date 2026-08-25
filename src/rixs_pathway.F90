!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! PROGRAM: rixs_pathway
!
!> @author
!> Christian Vorwerk, Humboldt Universität zu Berlin.
!> Elias Richter, Humboldt Universität zu Berlin.
!
! DESCRIPTION:
!> This program determines the matrices \f$ t^{(1)} \f$ and \f$ t^{(2)} \f$ from
!> the results of two BSE calculations stored in core_output.h5 and 
!> optical_output.h5. The results are stored in data.h5.
! 
! REVISION HISTORY:
! 09 07 2020 - Added documentation
! 24 01 2025 - 2 Pol treatment
!------------------------------------------------------------------------------
program rixs_pathway
  use mod_phdf5
  use mod_mpi
  use mod_io
  use mod_matmul
  use mod_blocks
  use hdf5, only: hid_t

  implicit none
  integer :: nkmax, nu_optical, no_optical, nu_core, no_core, global_core, global_optical
  integer :: interdim(2)
  type(io) :: optical, core
  type(input) :: inputparam
  integer(4), allocatable :: koulims_comb(:,:)
  integer(4) :: blocks_, blocks2_, blocks3_
  character(1024) :: fname_core, fname_optical,fname_pmat, fname_output, &
   & fname_inter, gname, gname2, ik, datasetname 
  integer(4) :: nblocks_, blsz_,nk_, k, blsz_k
  type(block1d) :: t1_b, t_b, evalsc_b
  type(block2d) :: eigvec_b, t2_b, prod_b
  real(8) ::  test
  complex(8) :: alpha, beta
  !MPI variables
  ! PHDF5 variables
  integer(hid_t) :: core_id, optical_id, pmat_id, inter_id
  integer(hid_t) :: t1_id, t2_id, evalsc_id
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
  ! open HDF5 files
  call phdf5_open_file(fname_core,core_id)
  call phdf5_open_file(fname_optical,optical_id)
  call phdf5_open_file(fname_pmat,pmat_id)
  ! create HDF5 files for output and intermediate data
  call phdf5_create_file(fname_inter,inter_id)

  ! initialize io objects for core and optics
  call get_koulims(optical,optical_id)
  call get_koulims(core,core_id)
  call get_smap(optical,optical_id)
  call get_smap(core,core_id)
  call get_ensortidx(core,core_id)
  call get_ensortidx(optical,optical_id)
  call set_param(optical)
  call set_param(core)
  call get_ismap(optical)
  call get_ismap(core)
  ! read input file
  call read_inputfile(inputparam)
  ! get number of k-grid (has to be the same for optical and core calculation)
  interdim=shape(core%koulims)
  nkmax=interdim(2)
  ! get global sizes
  nu_optical=optical%koulims(2,1)-optical%koulims(1,1)+1
  no_optical=optical%koulims(4,1)-optical%koulims(3,1)+1
  global_optical=nu_optical*no_optical*nkmax

  nu_core=core%koulims(2,1)-core%koulims(1,1)+1
  no_core=core%koulims(4,1)-core%koulims(3,1)+1
  global_core=nu_core*no_core*nkmax
  ! create combined map for valence-core transitions
  allocate(koulims_comb(4,nkmax))
  koulims_comb(1,:)=optical%koulims(3,:)
  koulims_comb(2,:)=optical%koulims(4,:)
  koulims_comb(3,:)=core%koulims(3,:)
  koulims_comb(4,:)=core%koulims(4,:)
  
  ! test whether the blocksize is possible
  test=float(nkmax)/float(inputparam%nblocks)
  
  if (float(floor(test)) .ne. test) then
    ! every rank evaluates this (purely local, deterministic) check and
    ! terminates together; only rank 0 prints, avoiding a deadlock where
    ! rank 0 would stop here while other ranks proceed into collective calls
    if (mpiglobal%rank .eq. 0) then
      print *, 'Blocksize', inputparam%nblocks, 'not compatible with ', nkmax, 'k-points'
    end if
    call terminate
  end if
  ! get number of k-points per block 
  nk_=nkmax/inputparam%nblocks
  !blsz_k=nu*no
  !blsz_=nu*no*nk_
  nblocks_=inputparam%nblocks
  ! set up the dataset for t(1)
  matsize_=(/ inputparam%nstatc /)
  datasetname='t(1)'
  call phdf5_setup_write(1,matsize_,.true.,trim(datasetname),'/',inter_id,t1_id)
  datasetname='evals'  
  call phdf5_setup_write(1,matsize_,.false.,trim(datasetname),'/',inter_id,evalsc_id)     
  !----------------------------------------------------!
  !    Calculate and Store blocks of the t(1) vector   !
  !----------------------------------------------------!
  do blocks_= firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
    !set-up for the blocks of t(1)
    t1_b = make_block1d(blocks_, inputparam%nstatc, nblocks_)
    !set-up for the block of core eigenvalues
    evalsc_b = make_block1d(blocks_, inputparam%nstatc, nblocks_)

    ! get block of core eigenvalues and write it to intermediate file
    if (.not. inputparam%ip_c) then
      call get_evals(evalsc_b,core_id)
    else
      call get_evalsIP(evalsc_b,core_id)
    end if
    call put_block1d(evalsc_b,evalsc_id)
    ! prepare output
    if (allocated(t1_b%zcontent)) deallocate(t1_b%zcontent)
    allocate(t1_b%zcontent(t1_b%blocksize))
    t1_b%zcontent(:)=cmplx(0.0d0, 0.0d0)

    ! 2nd loop over the blocks
    do blocks2_= 1, nblocks_
      ! set-up for the blocks of core eigenvectors
      eigvec_b = make_block2d(blocks2_, global_core, blocks_, inputparam%nstatc, nblocks_)
      ! set-up for the blocks of t
      t_b = make_block1d(blocks2_, global_core, nblocks_)
      t_b%nk=nk_
      t_b%kl=(blocks2_-1)*nk_+1
      t_b%ku=blocks2_*nk_
      
      ! generate block of X
      if (inputparam%ip_c) then
        call get_eigvecsIP(eigvec_b, core)
      else
        call get_eigvecs(eigvec_b, core_id)
      end if
      ! generate block of t
      call generate_t(t_b, core%koulims, core%smap, core%ismap, inputparam%pol_in, pmat_id)
      ! matrix-vector multiplication
      alpha=cmplx(1.0d0, 0.0d0)
      beta=cmplx(1.0d0, 0.0d0)
      call zgemm('T','N', eigvec_b%blocksize(2), 1, eigvec_b%blocksize(1), alpha, eigvec_b%zcontent, &
        & eigvec_b%blocksize(1), t_b%zcontent, t_b%blocksize, beta, t1_b%zcontent, t1_b%blocksize)
    end do
    call put_block1d(t1_b,t1_id)
  end do
  call phdf5_cleanup(t1_id)
  call phdf5_cleanup(evalsc_id)
  !----------------------------------------------------!
  !    Calculate and Store blocks of the t(2) vector   !
  !----------------------------------------------------!
  matsize2_=(/ inputparam%nstato, inputparam%nstatc /)
  datasetname='t(2)'
  call phdf5_setup_write(2,matsize2_,.true.,trim(datasetname),'/',inter_id,t2_id)     
  do blocks_= firstofset(mpiglobal%rank, nblocks_), lastofset(mpiglobal%rank, nblocks_)
    do blocks2_=1, nblocks_
      ! set-up block for t(2) matrix
      t2_b = make_block2d(blocks_, inputparam%nstato, blocks2_, inputparam%nstatc, nblocks_)
      t2_b%nk=nk_
      t2_b%k1l=(blocks_-1)*nk_+1
      t2_b%k1u=blocks_*nk_
      t2_b%k2l=(blocks2_-1)*nk_+1
      t2_b%k2u=blocks2_*nk_
     
      ! prepare output array
      if (allocated(t2_b%zcontent)) deallocate(t2_b%zcontent)
      allocate(t2_b%zcontent(t2_b%blocksize(1), t2_b%blocksize(2)))
      t2_b%zcontent(:,:)=cmplx(0.0d0, 0.0d0)
      do blocks3_=1, nblocks_
        ! set-up block for prod vector
        prod_b = make_block2d(blocks3_, global_optical, blocks2_, inputparam%nstatc, nblocks_)
        prod_b%nk=nk_
        prod_b%k1l=(blocks3_-1)*nk_+1
        prod_b%k1u=blocks3_*nk_

        ! set up block of optical eigenvectors
        eigvec_b = make_block2d(blocks3_, global_optical, blocks_, inputparam%nstato, nblocks_)
        ! generate block of eigenvectors
        if (inputparam%ip_o) then
          call get_eigvecsIP(eigvec_b, optical)
        else
          call get_eigvecs(eigvec_b, optical_id)
        end if
        ! generate block of intermediate product
        call generate_product(prod_b, inputparam, core, optical, core_id, pmat_id)

        !matrix-matrix multiplication
        alpha=cmplx(1.0d0, 0.0d0)
        beta=cmplx(1.0d0, 0.0d0)
        call zgemm('T', 'N', eigvec_b%blocksize(2), prod_b%blocksize(2), eigvec_b%blocksize(1), alpha, eigvec_b%zcontent, &
        & eigvec_b%blocksize(1), prod_b%zcontent, prod_b%blocksize(1), beta, t2_b%zcontent, t2_b%blocksize(1))
      end do ! blocks3_
      ! write block of t(2) into file
      call put_block2d(t2_b, t2_id)
    end do ! blocks2_
  end do ! blocks_ 
  call phdf5_cleanup(t2_id)
  
  call phdf5_close_file(core_id)
  call phdf5_close_file(optical_id)
  call phdf5_close_file(pmat_id)
  call phdf5_close_file(inter_id)
  ! close log file
  if (rank .eq. 0) close(7)
  !close HDF5 files
  call phdf5_finalize()
  call finitmpi()
end program rixs_pathway

