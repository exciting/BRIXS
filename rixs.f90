program rixs
  use mod_hdf5
  use mod_io
  use mod_rixs
  use mod_matmul

  implicit none
  real(8) :: broad, broad2
  real(8), allocatable :: omega(:), omega2(:)
  real(8) :: pol(3)
  integer :: w, nkmax, nu, no, k, w1, w2
  integer :: interdim(2), hamsiz, lambda
  type(io) :: optical, core
  type(input) :: inputparam
  complex(8), allocatable :: chi_optical(:,:,:), B(:,:), tprime(:,:,:)
  complex(8), allocatable :: B_matrix(:,:,:),A(:,:), A_matrix(:,:,:),A_inter(:)
  complex(8), allocatable :: inter(:), oscstr(:,:)
  complex(8) :: inter2
  real(8), allocatable :: results(:,:)
  integer(4), allocatable :: koulims_comb(:,:)
  character(1024) :: fname_core, fname_optical,fname_pmat
  
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
  ! generate optical chi
  call generate_chi(omega2,broad2,optical,fname_optical,chi_optical)
  ! generate B_vector
  call generate_Bvector(omega,broad,core,fname_pmat,fname_core,pol,B)
  ! generate a combined maps for t'
  interdim=shape(core%koulims)
  nkmax=interdim(2)
  allocate(koulims_comb(4,nkmax))
  koulims_comb(1,:)=optical%koulims(3,:)
  koulims_comb(2,:)=optical%koulims(4,:)
  koulims_comb(3,:)=core%koulims(3,:)
  koulims_comb(4,:)=core%koulims(4,:)
  ! generate t' matrix
  call generate_tprime(pol,koulims_comb,fname_pmat,tprime)
  !allocate A matrix
  ! get sizes from optical io type
  nu=optical%koulims(2,1)-optical%koulims(1,1)+1
  no=optical%koulims(4,1)-optical%koulims(3,1)+1
  interdim=shape(optical%koulims)
  nkmax=interdim(2)

  interdim=shape(optical%smap)
  hamsiz=interdim(2)
  allocate(A_matrix(nu,no,nkmax))
  allocate(A(hamsiz,size(omega)))
  do w=1,size(omega)
    ! generate B matrix from vector
    call transform2matrix(core%koulims,core%smap,B(:,w),B_matrix)
    do k=1, nkmax
      call matprod(B_matrix(:,:,k),tprime(:,:,k),A_matrix(:,:,k))
    end do
    ! generate A vector from matrix
    call transform2vector(optical%koulims,optical%smap,A_matrix,A_inter)
    A(:,w)=A_inter(:)
  end do
  deallocate(A_matrix,A_inter)
  ! allocate final spectrum
  allocate(inter(hamsiz))
  allocate(results(size(omega),size(omega2)))
  allocate(oscstr(hamsiz,size(omega)))
 
  !generate final spectrum
  do w1=1,size(omega)
    do w2=1,size(omega2)
      call matprod(chi_optical(:,:,w2),A(:,w1),inter)
      call matprod(A(:,w1),inter(:),inter2)
      results(w1,w2)=aimag(inter2)
    end do
  end do
  
  ! write results to file
  open(unit=4,file="results.txt",action="write",status="replace")
  write(4,'( *(2X, f14.6)\ )') (omega(w1), w1=1,size(omega))
  do w2=1,size(omega2)
    write(4,'( *(2X, f14.6)\ )') omega2(w2),(results(w1,w2), w1=1,size(omega))
  end do
  close(4)
  
  if (inputparam%oscstr) then 
    ! generate oscillator strenght
    call generate_oscstr(optical,A,oscstr)
    ! write oscillator strength to file
    open(unit=6,file="exciton.txt",action="write",status="replace")
    do lambda=1,hamsiz
      write(6,'( *(2X, f14.6)\ )') optical%evals(lambda)*27.211d0,(abs(oscstr(lambda,w1)), w1=1,size(omega))
    end do
    close(6)
  end if
 
  call hdf5_finalize()

end program
