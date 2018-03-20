program rixs
  use mod_hdf5
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
  character(1024) :: fname_core, fname_optical,fname_pmat
  integer(4) :: nblocks_, blsz_,nk_, k, k2, blsz_k
  complex(8), allocatable :: oscstr_b(:,:)
  type(block1d) :: vecA_b, evals_b
  type(block2d) :: evecs_b
  complex(8) :: alpha, beta
  real(8) :: start, finish, test
  call hdf5_initialize()
  
  !Specify file/dataset name
  fname_core='./core_output.h5'
  fname_optical='./optical_output.h5'
  fname_pmat='./pmat.h5'
  ! initialize io objects for core and optics
  call get_koulims(optical,fname_optical)
  call get_koulims(core,fname_core)
  call get_smap(optical,fname_optical)
  call get_smap(core,fname_core)
  call set_param(optical)
  call set_param(core)
  call get_ismap(optical)
  call get_ismap(core)
  ! inspect input objects
  write(*,*) '**********************************'
  write(*,*) '*      X-ray BSE Calculation     *'
  write(*,*) '**********************************'
  call inspect_h5(core)
  write(*,*) '**********************************'
  write(*,*) '*     Optics BSE Calculation     *'
  write(*,*) '**********************************'
  call inspect_h5(optical)
  ! read input file
  call read_inputfile(inputparam,'./rixs.in')
  ! set parameters
  ! for now, only ONE broadening parameter is use
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
  allocate(koulims_comb(4,nkmax))
  koulims_comb(1,:)=optical%koulims(3,:)
  koulims_comb(2,:)=optical%koulims(4,:)
  koulims_comb(3,:)=core%koulims(3,:)
  koulims_comb(4,:)=core%koulims(4,:)
  nu=optical%koulims(2,1)-optical%koulims(1,1)+1
  no=optical%koulims(4,1)-optical%koulims(3,1)+1
  ! test whether the blocksize is possible
  test=float(nkmax)/float(inputparam%nblocks)
  if (float(floor(test)) .ne. test) then
    print *, 'Blocksize', inputparam%nblocks, 'not compatible with ', nkmax, 'k-points'
    stop
  end if 
  ! define blocks for oscillator strenght
  nk_=nkmax/inputparam%nblocks
  blsz_k=nu*no
  blsz_=nu*no*nk_
  nblocks_=inputparam%nblocks
  print *, 'calculating ', nblocks_, 'blocks of size ', blsz_
  allocate(oscstr_b(blsz_,size(omega)))
  oscstr_b(:,:)=0.0d0
  ! open file
  open(unit=6,file="exciton.txt",action="write",status="replace")
  print *, 'start loop over blocks'
  !do k=1, nblocks_
  do k=1, nblocks_
    print *, 'calculating block ', k, '/', nblocks_
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
    ! generate block of oscillator strength
    oscstr_b(:,:)=0.0d0 
    print *, 'start loop over subblocks'
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
      call get_eigvecs2D_b(evecs_b,fname_optical)
      do w1=1, size(omega)
        ! generate block of A vector
        call cpu_time(start)
        call generate_Avector_b(vecA_b,omega(w1),broad,core,optical,fname_pmat,fname_core,pol)
        call cpu_time(finish)
        ! generate block of oscstr
        alpha=1.0d0
        beta=1.0d0
        call zgemm('C','N',blsz_,1,blsz_,alpha,evecs_b%zcontent,blsz_,vecA_b%zcontent,blsz_,beta,oscstr_b(:,w1),blsz_)
      end do ! w1
    end do ! k2
    do lambda=1, blsz_
      write(6,'( *(2X, f14.6)\ )') evals_b%dcontent(lambda)*27.211d0,(abs(oscstr_b(lambda,w1))**2, w1=1,size(omega))
    end do ! lambda
  end do ! k  
  !close file and hdf5 interface
  close(6)
  call hdf5_finalize()

end program
