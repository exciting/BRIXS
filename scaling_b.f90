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
  integer :: nkmax, nu, no
  integer :: interdim(2)
  type(io) :: optical, core
  type(input) :: inputparam
  integer(4), allocatable :: koulims_comb(:,:)
  character(1024) :: fname_core, fname_optical,fname_pmat
  integer(4) :: nblocks_, blsz_, blsz2_, k, k2, k3, i, j
  type(block2d) ::chi_core_b, eigvec_b, matB_b, tprime_b, matA_b
  type(block1d) :: eval_b, t_block, vecB_b, vecA_b
  complex(8), allocatable :: chi_core(:,:,:), pmat(:), B(:,:)
  complex(8), allocatable :: B_matrix(:,:,:), tprime(:,:,:)
  integer(8) :: dims(2)
  real(8) :: start, finish, cummulant, cummulant2, cummulant3

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
  ! define blocks for oscillator strenght
  nblocks_=nkmax
  blsz_=nu*no
  nu=core%koulims(2,1)-core%koulims(1,1)+1
  no=core%koulims(4,1)-core%koulims(3,1)+1
  blsz2_=nu*no
  print *, 'calculating ', nblocks_, 'blocks of size ', blsz_
  ! calculate block (k,k2) block of chi_core
  k=1
  k2=5
  do k3=1,10
  blsz2_=k3*10
  cummulant=0.0d0
  cummulant2=0.0d0
  cummulant3=0.0d0
  do k2=1, nblocks_
  do k=1, nblocks_
    ! set eigvec_b block
    eigvec_b%nblocks=nblocks_
    eigvec_b%blocksize=blsz2_
    eigvec_b%il=1
    eigvec_b%iu=blsz2_
    eigvec_b%jl=1
    eigvec_b%ju=blsz2_
    eigvec_b%offset=(/0, 0/)
    eigvec_b%id=(/1,1/)
    call cpu_time(start)
    call get_eigvecs2D_b_quick(eigvec_b,fname_core)
    call cpu_time(finish)
    cummulant=cummulant+finish-start
    call cpu_time(start)
    call get_eigvecs2D_b(eigvec_b,fname_core)
    call cpu_time(finish)
    cummulant2=cummulant2+finish-start
    call cpu_time(start)
    call get_eigvecs2D_b_quick2(eigvec_b,fname_core)
    call cpu_time(finish)
    cummulant3=cummulant3+finish-start
  end do
  end do
  print *, 'average run time:', cummulant/(nblocks_*nblocks_)
  print *, 'average run time2:', cummulant2/(nblocks_*nblocks_)
  print *, 'average run time3:', cummulant3/(nblocks_*nblocks_)
  end do
 
  !end do
  ! calculate block of chi_core for omega(1)
  !call generate_chi_block(chi_core_b,omega(1),broad,fname_core)
  ! generate block of B vector
  !call generate_Bvector_b(vecB_b,omega(1),broad,core,fname_pmat,fname_core,pol)
  ! generate block of B matrix
  !call transform2matrix_b(core%koulims,core%smap,vecB_b,matB_b)
  ! generate block of tprime
  !call generate_tprime_block(tprime_b,pol,koulims_comb,fname_pmat)
  
  ! generate full t vector
  !call generate_t(core%koulims,core%smap,core%ismap,pol,fname_pmat,pmat)
  ! get the full list of eigenvalues
  !call get_evals(core, fname_core)
  ! get full matrix of eigenvectors
  !do k=1,nblocks_
  call cpu_time(start)
  call get_eigvecs(core,fname_core)
  call cpu_time(finish)
  print *,'global eigvecs:', (finish-start)*1000.0d0, 'msec'
  !end do
  ! calculate the full chi_core
  !call generate_chi(omega,broad,core,fname_core,chi_core)
  ! calculate the full B vector
  !call generate_Bvector(omega,broad,core,fname_pmat,fname_core,pol,B)
  ! generate full B matrix
  !call transform2matrix(core%koulims,core%smap,B(:,1),B_matrix)
  ! generate tprime matrix
  !call generate_tprime(pol, koulims_comb,fname_pmat,tprime)
  ! print difference
  !do i=1, blsz2_
  !  print *, 'eval(',i,')=', eval_b%dcontent(i)-core%evals(i+eval_b%offset)
  !  do j=1,blsz2_
  !    print *, 'eigvec(', i, ',', j,')=', abs(eigvec_b%zcontent(i,j)-core%eigvecs(i+chi_core_b%offset(1),j+chi_core_b%offset(2)))
  !    print *, 'chi_core(', i, ',', j,')=', abs(chi_core_b%zcontent(i,j)-chi_core(i+chi_core_b%offset(1),j+chi_core_b%offset(2),1))
  !    print *, '(', i,')=', eval_b%dcontent(i),core%evals(i+eval_b%offset), eval_b%dcontent(i)-core%evals(i+eval_b%offset)
  !    print *, 'B_matrix(', i, ',', j,')=', abs(matB_b%zcontent(i,j)-B_matrix(i+matB_b%offset(1),j+matB_b%offset(2),k))
  !  end do
  !end do
  !dims=shape(tprime_b%zcontent)
  !do i=1, dims(1)
    !do j=1, dims(2)
    !  print *, 'tprime(', i, ',', j,')=', abs(tprime_b%zcontent(i,j)-tprime(i,j,k))
    !end do
  !end do
  !print *, 'shape tests'
  !print *, shape(tprime_b%zcontent),'=', shape(tprime(:,:,1))
  !print *, shape(matB_b%zcontent),'=', shape(B_matrix(:,:,1))
  !do i=1, blsz2_
  !    print *, 'pmat(', i, ')=', abs(t_block%zcontent(i)-pmat(i+t_block%offset))
  !    print *, 'vecB(', i, ')=', abs(vecB_b%zcontent(i)-B(i+vecB_b%offset,1))
  !end do
  !print *, 'complex(8) epsilon=', epsilon(eval_b%dcontent(1)) 
  call hdf5_finalize()

end program
