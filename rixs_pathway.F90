program rixs_pathway
  use mod_phdf5
  use modmpi
  use mod_io
  use mod_rixs
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
  
  if ((float(floor(test)) .ne. test) .and. (mpiglobal%rank .eq. 0)) then
    print *, 'Blocksize', inputparam%nblocks, 'not compatible with ', nkmax, 'k-points'
    stop
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
    t1_b%nblocks=nblocks_
    t1_b%blocksize=nofblock(blocks_, inputparam%nstatc, nblocks_)
    t1_b%global=inputparam%nstatc
    t1_b%il=firstofblock(blocks_, inputparam%nstatc, nblocks_)
    t1_b%iu=lastofblock(blocks_, inputparam%nstatc, nblocks_)
    t1_b%offset=firstofblock(blocks_, inputparam%nstatc, nblocks_)-1
    t1_b%id=blocks_
    !set-up for the block of core eigenvalues
    evalsc_b%nblocks=nblocks_
    evalsc_b%blocksize=nofblock(blocks_, inputparam%nstatc, nblocks_)
    evalsc_b%global=inputparam%nstatc
    evalsc_b%il=firstofblock(blocks_, inputparam%nstatc, nblocks_)
    evalsc_b%iu=lastofblock(blocks_, inputparam%nstatc, nblocks_)
    evalsc_b%offset=firstofblock(blocks_, inputparam%nstatc, nblocks_)-1
    evalsc_b%id=blocks_

    ! get block of core eigenvalues and write it to intermediate file
    call get_evals_block(evalsc_b,core_id)
    call put_block1d(evalsc_b,evalsc_id)
    ! prepare output
    if (allocated(t1_b%zcontent)) deallocate(t1_b%zcontent)
    allocate(t1_b%zcontent(t1_b%blocksize))
    t1_b%zcontent(:)=cmplx(0.0d0, 0.0d0)

    ! 2nd loop over the blocks
    do blocks2_= 1, nblocks_
      ! set-up for the blocks of core eigenvectors
      eigvec_b%nblocks=nblocks_
      eigvec_b%blocksize=(/ nofblock(blocks2_, global_core, nblocks_), t1_b%blocksize /)
      eigvec_b%global=(/ global_core, inputparam%nstatc /)
      eigvec_b%il=firstofblock(blocks2_, global_core, nblocks_)
      eigvec_b%iu=lastofblock(blocks2_, global_core, nblocks_)
      eigvec_b%jl=t1_b%il
      eigvec_b%ju=t1_b%iu
      eigvec_b%offset(1)=firstofblock(blocks2_, global_core, nblocks_)-1
      eigvec_b%offset(2)=t1_b%offset
      eigvec_b%id=(/ blocks2_, blocks_ /)
      ! set-up for the blocks of t
      t_b%nblocks=nblocks_
      t_b%blocksize=nofblock(blocks2_, global_core, nblocks_)
      t_b%global=global_core
      t_b%il=firstofblock(blocks2_, global_core, nblocks_)
      t_b%iu=lastofblock(blocks2_, global_core, nblocks_)
      t_b%nk=nk_
      t_b%kl=(blocks2_-1)*nk_+1
      t_b%ku=blocks2_*nk_
      t_b%offset=firstofblock(blocks2_, global_core, nblocks_)-1
      t_b%id=blocks2_
      
      ! generate block of X
      call get_eigvecs2D_b(eigvec_b, core_id)
      ! generate block of t
      call generate_tblock(t_b, core%koulims, core%smap, core%ismap, inputparam%pol, pmat_id)
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
      t2_b%nblocks=nblocks_
      t2_b%blocksize=(/ nofblock(blocks_, inputparam%nstato, nblocks_), nofblock(blocks2_, inputparam%nstatc, nblocks_) /)
      t2_b%global=(/ inputparam%nstato, inputparam%nstatc /)
      t2_b%nk=nk_
      t2_b%il=firstofblock(blocks_, inputparam%nstato, nblocks_)
      t2_b%iu=lastofblock(blocks_, inputparam%nstato, nblocks_)
      t2_b%jl=firstofblock(blocks2_, inputparam%nstatc, nblocks_)
      t2_b%ju=lastofblock(blocks2_, inputparam%nstatc, nblocks_)
      t2_b%k1l=(blocks_-1)*nk_+1
      t2_b%k1u=blocks_*nk_
      t2_b%k2l=(blocks2_-1)*nk_+1
      t2_b%k2u=blocks2_*nk_
      t2_b%offset=(/ t2_b%il-1, t2_b%jl-1 /)
      t2_b%id=(/ blocks_, blocks2_ /)
     
      ! prepare output array
      if (allocated(t2_b%zcontent)) deallocate(t2_b%zcontent)
      allocate(t2_b%zcontent(t2_b%blocksize(1), t2_b%blocksize(2)))
      t2_b%zcontent(:,:)=cmplx(0.0d0, 0.0d0)
      do blocks3_=1, nblocks_
        ! set-up block for prod vector
        prod_b%nblocks=nblocks_
        prod_b%blocksize=(/ nofblock(blocks3_, global_optical, nblocks_), nofblock(blocks2_, inputparam%nstatc, nblocks_)  /)
        prod_b%global=(/ global_optical, inputparam%nstatc /)
        prod_b%nk=nk_
        prod_b%il=firstofblock(blocks3_, global_optical, nblocks_)
        prod_b%iu=lastofblock(blocks3_, global_optical, nblocks_)
        prod_b%jl=firstofblock(blocks2_, inputparam%nstatc, nblocks_)
        prod_b%ju=lastofblock(blocks2_, inputparam%nstatc, nblocks_)
        prod_b%k1l=(blocks3_-1)*nk_+1
        prod_b%k1u=blocks3_*nk_
        prod_b%offset=(/ prod_b%il-1, prod_b%jl-1 /)
        prod_b%id=(/ blocks3_, blocks2_ /)

        ! set up block of optical eigenvectors
        eigvec_b%nblocks=nblocks_
        eigvec_b%blocksize=(/ nofblock(blocks3_, global_optical, nblocks_), t2_b%blocksize(1) /)
        eigvec_b%global=(/ global_optical, inputparam%nstato /)
        eigvec_b%il=firstofblock(blocks3_, global_optical, nblocks_)
        eigvec_b%iu=lastofblock(blocks3_, global_optical, nblocks_)
        eigvec_b%jl=t2_b%il
        eigvec_b%ju=t2_b%iu
        eigvec_b%offset(1)=firstofblock(blocks3_, global_optical, nblocks_)-1
        eigvec_b%offset(2)=t2_b%offset(1)
        eigvec_b%id=(/ blocks3_, blocks_ /)
        ! generate block of eigenvectors
        call get_eigvecs2D_b(eigvec_b, optical_id)
        ! generate block of intermediate product
        call gen_prod_b(prod_b, inputparam, core, optical, core_id, pmat_id)

        !matrix-matrix multiplication
        alpha=cmplx(1.0d0, 0.0d0)
        beta=cmplx(1.0d0, 0.0d0)
        call zgemm('C', 'N', eigvec_b%blocksize(2), prod_b%blocksize(2), eigvec_b%blocksize(1), alpha, eigvec_b%zcontent, &
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

