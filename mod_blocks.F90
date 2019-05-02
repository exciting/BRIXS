module mod_blocks
  implicit none

  type :: block1d
    integer(4) :: nblocks     ! number of blocks
    integer(4) :: blocksize   ! length of each block
    integer(4) :: global      ! global vector size
    integer(4) :: nk          ! number of k-pts per block
    integer(4) :: il, iu      ! absolute first & last transition index within block
    integer(4) :: kl, ku      ! absolute first & last k-index within block
    integer(4) :: offset      ! offset
    integer(4) :: id  ! ID (in this case number) of the block
    complex(8), allocatable :: zcontent(:)
    real(8), allocatable :: dcontent(:)
  end type block1d
  
  type :: block2d
    integer(4) :: nblocks     ! number of block
    integer(4), dimension(2) :: blocksize   ! non-square 2D blocksize
    integer(4), dimension(2) :: global      ! global 2D size
    integer(4) :: nk          ! number of k-pts per block
    integer(4) :: il, iu, jl, ju ! absolute transition index ranges within the block  
    integer(4) :: k1l, k1u, k2l, k2u ! absolute k-index within the block
    integer(4), dimension(2) :: offset    ! offset
    integer(4), dimension(2) :: id    ! 2-D ID of block
    complex(8), allocatable :: zcontent(:,:)
    real(8), allocatable :: dcontent(:,:)
  end type block2d

  type :: block3d
    integer(4) :: nblocks                 ! number of blocks
    integer(4) :: nk                      ! number of k-pts / block
    integer(4), dimension(3) :: blocksize ! non-square 3D blocksize
    integer(4) :: kl, ku                  ! absolute k-index within block
    integer(4) :: id                      ! store the id of the block it was generated from
    complex(8), allocatable :: zcontent(:,:,:)
  end type block3d
  contains
  ! Methodenbereich
    !-----------------------------------------------------------------------------
    subroutine transform2matrix_b(inkoulims,insmap,inbl1d,outbl3d)
      implicit none
      integer(4), intent(in) :: inkoulims(:,:)
      integer(4), intent(in) :: insmap(:,:)
      type(block1d), intent(in) :: inbl1d
      type(block3d), intent(out) :: outbl3d

      ! local variables
      integer(4) :: lu,uu,lo,uo,nu,no,nk0,nkmax,i
      integer(4) :: dim_koulims(2),dim_smap(2),hamsiz,i1,i2, i3
  
      !get shapes
      dim_koulims=shape(inkoulims)
      dim_smap=shape(insmap)
      !determine sizes
      lu=inkoulims(1,1)
      uu=inkoulims(2,1)
      lo=inkoulims(3,1)
      uo=inkoulims(4,1)
      nu=uu-lu+1
      no=uo-lo+1
      nk0=insmap(3,1)
      nkmax=dim_koulims(2)
      hamsiz=dim_smap(2)
      ! generate the block output
      outbl3d%nblocks=inbl1d%nblocks
      outbl3d%nk=inbl1d%nk
      outbl3d%blocksize=(/nu,no,inbl1d%nk/)
      outbl3d%kl=inbl1d%kl
      outbl3d%ku=inbl1d%ku
      outbl3d%id=inbl1d%id
      ! allocate the output matrix
      if (allocated(outbl3d%zcontent)) deallocate(outbl3d%zcontent)
      allocate(outbl3d%zcontent(nu, no, inbl1d%nk))
      ! loop over all transitions
      do i=inbl1d%il, inbl1d%iu
        i1=insmap(1,i)-lu+1
        i2=insmap(2,i)-lo+1
        i3=insmap(3,i)-nk0+1
        outbl3d%zcontent(i1,i2,i3-inbl1d%kl+1)=inbl1d%zcontent(i-inbl1d%offset)
      end do
    end subroutine 
    !-----------------------------------------------------------------------------
    subroutine transform_matrix2matrix(inkoulims,insmap,inbl2d,out4d)
      implicit none
      integer(4), intent(in) :: inkoulims(:,:)
      integer(4), intent(in) :: insmap(:,:)
      type(block2d), intent(in) :: inbl2d
      complex(8), allocatable, intent(out) :: out4d(:,:,:,:)

      ! local variables
      integer(4) :: lu,uu,lo,uo,nu,no,nk0,nkmax,i, lambda
      integer(4) :: dim_koulims(2),dim_smap(2),hamsiz,i1,i2, i3
  
      !get shapes
      dim_koulims=shape(inkoulims)
      dim_smap=shape(insmap)
      !determine sizes
      lu=inkoulims(1,1)
      uu=inkoulims(2,1)
      lo=inkoulims(3,1)
      uo=inkoulims(4,1)
      nu=uu-lu+1
      no=uo-lo+1
      nk0=insmap(3,1)
      nkmax=dim_koulims(2)
      hamsiz=dim_smap(2)
      ! allocate the output matrix
      if (allocated(out4d)) deallocate(out4d)
      allocate(out4d(nu, no, inbl2d%nk,inbl2d%blocksize(2)))
      ! loop over all excitons
      do lambda=1, inbl2d%blocksize(2)
        ! loop over all transitions
        do i=inbl2d%il, inbl2d%iu
          i1=insmap(1,i)-lu+1
          i2=insmap(2,i)-lo+1
          i3=insmap(3,i)-nk0+1
          out4d(i1,i2,i3-inbl2d%k1l+1,lambda)=inbl2d%zcontent(i-inbl2d%offset(1),lambda)
        end do !i
      end do !lambda
    end subroutine 
    !-----------------------------------------------------------------------------
    subroutine transform2vector_b(inkoulims,insmap,inbl3d,outbl1d)
      implicit none
      integer(4), intent(in) :: inkoulims(:,:)
      integer(4), intent(in) :: insmap(:,:)
      type(block3d), intent(in) :: inbl3d
      type(block1d), intent(inout) :: outbl1d

      ! local variables
      integer(4) :: lu,uu,lo,uo,nu,no,nk0,i, blsz_
      integer(4) :: dim_koulims(2),dim_smap(2),nkmax,hamsiz,i1,i2, i3
  
      !get shapes
      dim_koulims=shape(inkoulims)
      dim_smap=shape(insmap)
      !determine sizes
      lu=inkoulims(1,1)
      uu=inkoulims(2,1)
      lo=inkoulims(3,1)
      uo=inkoulims(4,1)
      nu=uu-lu+1
      no=uo-lo+1
      nk0=insmap(3,1)
      nkmax=dim_koulims(2)
      hamsiz=dim_smap(2)
      ! store shorthand for blocksize
      blsz_=inbl3d%blocksize(1)*inbl3d%blocksize(2)*inbl3d%blocksize(3)
      !generate the block output
      outbl1d%nblocks=inbl3d%nblocks
      outbl1d%blocksize=blsz_
      outbl1d%global=blsz_*inbl3d%nblocks
      outbl1d%nk=inbl3d%nk
      outbl1d%il=(inbl3d%id-1)*blsz_+1
      outbl1d%iu=(inbl3d%id)*blsz_
      outbl1d%kl=inbl3d%kl
      outbl1d%ku=inbl3d%ku
      outbl1d%offset=(inbl3d%id-1)*blsz_
      outbl1d%id=inbl3d%id

      ! allocate the output matrix
      if (allocated(outbl1d%zcontent)) deallocate(outbl1d%zcontent)
      allocate(outbl1d%zcontent(outbl1d%blocksize))

      ! loop over all transitions
      do i=outbl1d%il, outbl1d%iu
        i1=insmap(1,i)-lu+1
        i2=insmap(2,i)-lo+1
        i3=insmap(3,i)-nk0+1-inbl3d%kl+1
        outbl1d%zcontent(i-outbl1d%offset)=inbl3d%zcontent(i1,i2,i3)
      end do
    end subroutine 

    !-----------------------------------------------------------------------------
    subroutine transform_matrix2vector(input, in4d, outbl2d)
      use mod_io, only: io
      implicit none
      type(io), intent(in) :: input
      complex(8), allocatable, intent(in) :: in4d(:,:,:,:)
      type(block2d), intent(inout) :: outbl2d

      ! local variables
      integer(4) :: matsize_(4)
      integer(4) :: lu,uu,lo,uo,nu,no,nk0,i, blsz_, lambda
      integer(4) :: dim_koulims(2),dim_smap(2),nkmax,hamsiz,i1,i2, i3
  
      !get shapes
      dim_koulims=shape(input%koulims)
      dim_smap=shape(input%smap)
      !determine sizes
      lu=input%koulims(1,1)
      uu=input%koulims(2,1)
      lo=input%koulims(3,1)
      uo=input%koulims(4,1)
      nu=uu-lu+1
      no=uo-lo+1
      nk0=input%smap(3,1)
      nkmax=dim_koulims(2)
      hamsiz=dim_smap(2)
      
      !generate the block output
      if (allocated(outbl2d%zcontent)) deallocate(outbl2d%zcontent)
      allocate(outbl2d%zcontent(outbl2d%blocksize(1), outbl2d%blocksize(2)))

      ! loop over all transitions
      do lambda=1, outbl2d%blocksize(2)
        do i=outbl2d%il, outbl2d%iu
          i1=input%smap(1,i)-lu+1
          i2=input%smap(2,i)-lo+1
          i3=input%smap(3,i)-nk0+1-outbl2d%k1l+1
          outbl2d%zcontent(i-outbl2d%offset(1), lambda)=in4d(i1,i2,i3,lambda)
        end do !i
      end do !lambda
    end subroutine 

    !-----------------------------------------------------------------------------
    subroutine get_evals_block(inblock1d,file_id)
      use mod_phdf5, only: phdf5_setup_read, phdf5_cleanup, &
        &                  phdf5_read
      use hdf5, only: hid_t
      implicit none
      type(block1d), intent(inout) :: inblock1d
      integer(hid_t), intent(in) :: file_id
      !local variables
      character(len=1024) :: path, dsetname
      integer, dimension(1) :: dims_, offset_, dimsg_
      integer(hid_t) :: dataset_id
      !allocate output
      if (allocated(inblock1d%dcontent)) deallocate(inblock1d%dcontent)
      allocate(inblock1d%dcontent(inblock1d%blocksize))
    
      path='/eigvec-singlet-TDA-BAR-full/0001'
      dsetname='evals'
      ! get data
      dims_(1)=inblock1d%blocksize
      dimsg_(1)=inblock1d%global
      offset_(1)=inblock1d%offset
      ! open the dataset
      call phdf5_setup_read(1,dims_,.false.,dsetname,path,file_id,dataset_id)
      ! read data
      call phdf5_read(inblock1d%dcontent(1),dims_,dimsg_,offset_,dataset_id)
      ! close dataset
      call phdf5_cleanup(dataset_id)
    end subroutine
    
    !-----------------------------------------------------------------------------
    subroutine get_evalsIP_block(inblock1d,file_id)
      use mod_phdf5, only: phdf5_setup_read, phdf5_cleanup, &
        &                  phdf5_read
      use hdf5, only: hid_t
      implicit none
      type(block1d), intent(inout) :: inblock1d
      integer(hid_t), intent(in) :: file_id
      !local variables
      character(len=1024) :: path, dsetname
      integer, dimension(1) :: dims_, offset_, dimsg_
      integer(hid_t) :: dataset_id
      !allocate output
      if (allocated(inblock1d%dcontent)) deallocate(inblock1d%dcontent)
      allocate(inblock1d%dcontent(inblock1d%blocksize))
    
      path='/eigvec-singlet-TDA-BAR-full/0001'
      dsetname='evalsIP'
      ! get data
      dims_(1)=inblock1d%blocksize
      dimsg_(1)=inblock1d%global
      offset_(1)=inblock1d%offset
      ! open the dataset
      call phdf5_setup_read(1,dims_,.false.,dsetname,path,file_id,dataset_id)
      ! read data
      call phdf5_read(inblock1d%dcontent(1),dims_,dimsg_,offset_,dataset_id)
      ! close dataset
      call phdf5_cleanup(dataset_id)
  end subroutine

  !-----------------------------------------------------------------------------
  subroutine get_eigvecs_b(inblock2d,file_id)
    use mod_phdf5, only:  phdf5_setup_read, phdf5_read, phdf5_cleanup 
    use hdf5, only: hid_t
    implicit none
    type(block2d), intent(inout) :: inblock2d
    integer(hid_t), intent(in) :: file_id
    !local variables
    complex(8), allocatable :: eigvec_(:)
    integer(4) ::  i, offset_(1), dim_(1), dimsg_(1)
    integer(hid_t) :: dataset_id
    character(len=1024) :: path, dsetname
    character(256) :: ci
    ! allocate output
    if (allocated(inblock2d%zcontent)) deallocate(inblock2d%zcontent)
    allocate(inblock2d%zcontent(inblock2d%blocksize(1),inblock2d%blocksize(2))) 
    ! get data
    do i=1, inblock2d%blocksize(2)
      write(ci, '(I8.8)') i+inblock2d%offset(2)
      path='eigvec-singlet-TDA-BAR-full/0001/rvec'
      dsetname=trim(adjustl(ci))
      ! Get dimension of eigvec for given lambda
      !call hdf5_get_dims(fname,path,ci,dims_)
      ! Allocate intermediate eigenvector array
      if (allocated(eigvec_)) deallocate(eigvec_)
      allocate(eigvec_(inblock2d%blocksize(1)))
      !allocate(eigvec_(dims_(2)))
      ! Get data
      offset_(1)=inblock2d%offset(1)
      dim_(1)=inblock2d%blocksize(1)
      dimsg_(1)=inblock2d%global(1)
      !open dataset
      call phdf5_setup_read(1,dim_,.true.,dsetname,path,file_id,dataset_id)
      !read data
      call phdf5_read(eigvec_(1),dim_,dimsg_,offset_,dataset_id)
      !close dataset
      call phdf5_cleanup(dataset_id)
      ! Write data to final array
      inblock2d%zcontent(:,i)=eigvec_(:)
    end do
    deallocate(eigvec_)
  end subroutine 

  !-----------------------------------------------------------------------------
  subroutine get_eigvecsIP_b(inblock2d, io_in)
    use mod_io, only: io
    implicit none
    type(block2d), intent(inout) :: inblock2d
    type(io), intent(in) :: io_in
    !local variables
    integer(4) :: lambda, pos, i
    ! allocate output
    if (allocated(inblock2d%zcontent)) deallocate(inblock2d%zcontent)
    allocate(inblock2d%zcontent(inblock2d%blocksize(1),inblock2d%blocksize(2))) 
    inblock2d%zcontent(:,:)=cmplx(0.0d0,0.0d0)
    ! generate data
    do i=1, inblock2d%blocksize(2)
      !absolute excitonic index
      lambda= i+inblock2d%offset(2)
      !position of transition in transition space
      pos=io_in%ensortidx(lambda)
      if ((pos-inblock2d%offset(1)>0) .and. (pos-inblock2d%offset(1)<inblock2d%blocksize(1))) then
        inblock2d%zcontent(pos-inblock2d%offset(1),i)=cmplx(1.0d0,0.0d0)

      end if
    end do
  end subroutine

  !-----------------------------------------------------------------------------
  subroutine generate_tblock(inblock1d,koulims,smap,ismap,pol,file_id)
    use hdf5, only: hid_t
    use mod_phdf5
    implicit none
    type(block1d), intent(inout) :: inblock1d
    integer(4), intent(in) :: koulims(:,:)
    integer(4), intent(in) :: smap(:,:)
    integer(4), intent(in) :: ismap(:,:,:)
    real(8), intent(in) :: pol(3)
    integer(hid_t), intent(in) :: file_id
    !internal variables
    integer(4), dimension(2) :: dim_koulims, dim_smap
    integer(4) :: lu, uu, lo, uo, nu, no, nk0
    integer(4) :: k,i,j,dimensions(4), dimsg_(3), offset_(3)
    character(len=1024) :: path, dsetname, cik
    complex(8), allocatable :: pmat_(:,:,:)
    !complex(8) :: pmat_(3,2,35)
    integer :: stat_var
    integer(hid_t) :: dataset_id

    ! allocate output
    if (allocated(inblock1d%zcontent)) deallocate(inblock1d%zcontent)
    allocate(inblock1d%zcontent(inblock1d%blocksize))
    dsetname=trim(adjustl('pmat'))
    !get shapes
    dim_koulims=shape(koulims)
    dim_smap=shape(smap)
    !determine sizes
    lu=koulims(1,1)
    uu=koulims(2,1)
    lo=koulims(3,1)
    uo=koulims(4,1)
    nu=uu-lu+1
    no=uo-lo+1
    nk0=smap(3,1)
    do k=inblock1d%kl, inblock1d%ku
      write(cik, '(I4.4)') k
      !determine size of matrix in hdf5 file
      path=trim(adjustl('/pmat/'//trim(adjustl(cik))))
      call phdf5_get_dims(file_id,path,trim(adjustl(dsetname)),dimensions)
      dimsg_=(/dimensions(2), dimensions(3), dimensions(4)/)
      offset_=(/0, 0, 0/)
      !allocate intermediate transition matrix for each k-point
      if (allocated(pmat_)) deallocate(pmat_)
      allocate(pmat_(dimensions(2),dimensions(3),dimensions(4)), stat=stat_var)
      call phdf5_setup_read(3,dimsg_,.true.,dsetname,path,file_id,dataset_id)
      call phdf5_read(pmat_(1,1,1),dimsg_,dimsg_,offset_,dataset_id)
      call phdf5_cleanup(dataset_id)
      ! write  transition matrix into file for the states included 
      ! in the BSE calculation
      do i=1, no
        do j=1, nu
          inblock1d%zcontent(ismap(j,i,k)-inblock1d%offset)=conjg(pmat_(1,i+lo-1,j+lu-1))*pol(1)+&
            & conjg(pmat_(2,i+lo-1,j+lu-1))*pol(2)+&
            & conjg(pmat_(3,i+lo-1,j+lu-1))*pol(3)
        end do
      end do
    end do
    deallocate(pmat_)
  end subroutine 
!-----------------------------------------------------------------------------
  subroutine generate_tprime_block(in3d,pol,koulims,file_id)
    use mod_phdf5
    use hdf5, only: hid_t
    implicit none
    real(8), intent(in) :: pol(3)
    integer(4), intent(in) :: koulims(:,:)
    integer(hid_t), intent(in) :: file_id
    type(block3d), intent(inout) :: in3d
    ! local variables
    integer(4) :: lu, uu, lo, uo, nu, no, nkmax, i,j
    integer(4) :: dimensions(4), inter(2), nk_, k
    integer :: dimsg_(3), offset_(3)
    integer(hid_t) :: dataset_id
    character(1024) :: cik, path, dsetname
    complex(8), allocatable :: pmat_(:,:,:)
    !complex(8) :: pmat_(3,2,35)
    !determine sizes
    lu=koulims(1,1)
    uu=koulims(2,1)
    lo=koulims(3,1)
    uo=koulims(4,1)
    nu=uu-lu+1
    no=uo-lo+1
    inter=shape(koulims)
    nkmax=inter(2)
    nk_=in3d%nk
    ! allocate output array
    if (allocated(in3d%zcontent)) deallocate(in3d%zcontent)
    allocate(in3d%zcontent(no,nu,nk_))
    in3d%blocksize=(/no, nu, nk_/)
    dsetname='pmat'
    do k=in3d%kl, in3d%ku
      write(cik, '(I4.4)') k
      !determine size of matrix in hdf5 file
      path='/pmat/'//trim(adjustl(cik))
      call phdf5_get_dims(file_id,path,trim(adjustl(dsetname)),dimensions)
      dimsg_=(/dimensions(2), dimensions(3), dimensions(4)/)
      offset_=(/0, 0, 0/)
      !allocate intermediate transition matrix for each k-point
      if (allocated(pmat_)) deallocate(pmat_)
      allocate(pmat_(dimensions(2),dimensions(3),dimensions(4)))
      !open dataset
      call phdf5_setup_read(3,dimsg_,.true.,trim(adjustl(dsetname)),path,file_id,dataset_id)
      !read data
      call phdf5_read(pmat_(1,1,1),dimsg_,dimsg_,offset_,dataset_id)
      ! close dataset
      call phdf5_cleanup(dataset_id)
      do i=1,no
        do j=1,nu
          in3d%zcontent(i,j,k-in3d%kl+1)=pol(1)*pmat_(1,i+lo-1,j+lu-1)+&
            & pol(2)*pmat_(2,i+lo-1,j+lu-1)+&
            & pol(3)*pmat_(3,i+lo-1,j+lu-1)
        end do
      end do
     end do
    deallocate(pmat_)
  end subroutine 

  !-----------------------------------------------------------------------------
  subroutine generate_chi_block(inblock2d,omega,inputparam,file_id)
    use mod_io, only: io, input
    use hdf5, only: hid_t
    implicit none
    real(8), intent(in) :: omega
    type(input), intent(in) :: inputparam
    integer(hid_t), intent(in) :: file_id
    type(block2d), intent(inout) :: inblock2d
    ! internal variables
    type(block2d) :: block1, block2
    type(block1d) :: block3
    complex(8), allocatable :: inter(:,:), inter2(:,:)
    !complex(8) :: inter
    integer(4) :: blsz,i, k, reduced_
    complex(8) :: alpha, beta
    
    ! shorthand notation for blocksize
    blsz=inblock2d%blocksize(1)
    
    ! allocate chi_b
    if (allocated(inblock2d%zcontent)) deallocate(inblock2d%zcontent)
    allocate(inblock2d%zcontent(inblock2d%blocksize(1),inblock2d%blocksize(2)))
    inblock2d%zcontent(:,:)=0.0d0
    ! chi_block is calculated as follows:
    ! $\chi^{ij}=\sum_k \hat{X}^{ik}\times\hat{Y}^{ii}\times[\hat{X}^{ki}]$
    ! I use non-square blocks here, such that I can calculate the density-density
    ! response function from subsets of excitonic eigenstates. In this case, the 
    ! blocksize and number of blocks is not uniquely defined any more. 
    ! Here, I chose to use the same number of blocks as in the overall calculation
    do k=1, inblock2d%nblocks
      ! local matrix sizes
      reduced_=nofblock(k, inputparam%nstatc, inblock2d%nblocks)
      ! define the blocks for the summation
      ! block1 and block2 are blocks of the eigenvectors
      block1%nblocks=inblock2d%nblocks
      block1%blocksize=(/blsz, reduced_ /)
      block1%global=(/blsz*inblock2d%nblocks, inputparam%nstatc/)
      block1%il=inblock2d%il
      block1%iu=inblock2d%iu
      block1%jl=firstofblock(k, inputparam%nstatc, inblock2d%nblocks)
      block1%ju=lastofblock(k, inputparam%nstatc, inblock2d%nblocks)
      block1%offset(1)=inblock2d%offset(1)
      block1%offset(2)=firstofblock(k, inputparam%nstatc, inblock2d%nblocks)-1
      block1%id(1)=inblock2d%id(1)
      block1%id(2)=k

      block2%nblocks=inblock2d%nblocks
      block2%blocksize=(/blsz, reduced_ /)
      block2%global=(/blsz*inblock2d%nblocks, inputparam%nstatc/)
      block2%il=inblock2d%jl
      block2%iu=inblock2d%ju
      block2%jl=firstofblock(k, inputparam%nstatc, inblock2d%nblocks)
      block2%ju=lastofblock(k, inputparam%nstatc, inblock2d%nblocks)
      block2%k2l=(k-1)*block2%nk+1
      block2%k2u=k*block2%nk
      block2%offset(1)=inblock2d%offset(2)
      block2%offset(2)=firstofblock(k, inputparam%nstatc, inblock2d%nblocks)-1
      block2%id(1)=inblock2d%id(2)
      block2%id(2)=k
      ! block3 is a 1D block of the eigenvectors
      block3%nblocks=inblock2d%nblocks
      block3%blocksize=reduced_
      block3%global=inputparam%nstatc
      block3%il=firstofblock(k, inputparam%nstatc, inblock2d%nblocks)
      block3%iu=lastofblock(k, inputparam%nstatc, inblock2d%nblocks)
      block3%offset=firstofblock(k, inputparam%nstatc, inblock2d%nblocks) -1
      block3%id=k
       
      !get eigenvalues and eigenvectors
      call get_evals_block(block3,file_id)
      call get_eigvecs_b(block1,file_id)
      call get_eigvecs_b(block2,file_id)
    
      ! allocate array for hermetian conjugate of eigvecs
      if (allocated(inter)) deallocate(inter)
      if (allocated(inter2)) deallocate(inter2)
      allocate(inter(reduced_, reduced_))
      allocate(inter2(reduced_,blsz))
      inter(:,:)=0.0d0
      do i=1,reduced_
        inter(i,i)=-1.0d0/(block3%dcontent(i)*27.211d0-omega+cmplx(0.0d0,inputparam%broad))
      end do
      ! generate intermediate matrix
      !do i=1,blsz
      !  inter=-1.0d0/(block3%dcontent(i)*27.211d0-omega+cmplx(0.0d0,broad))
      !   do j=1,blsz
      !    block2%zcontent(i,j)=inter*block2%zcontent(i,j)
      !  end do
      !end do
      alpha=1.0d0
      beta=0.0d0
      call zgemm('n','c',reduced_,blsz,reduced_,alpha,inter,reduced_,block2%zcontent,&
          &        blsz,beta,inter2,reduced_)
      alpha=1.0d0
      beta=1.0d0
      call zgemm('N','N',blsz,blsz,reduced_,alpha,block1%zcontent,blsz,inter2,&
        &        reduced_,beta,inblock2d%zcontent(1:blsz,1:blsz),blsz)
      !alpha=1.0d0
      !beta=1.0d0
      !call zgemm('N','N',blsz,blsz,blsz,alpha,block1%zcontent(1:blsz,1:blsz),blsz,block2%zcontent,&
      !  &        blsz,beta,inblock2d%zcontent(1:blsz,1:blsz),blsz)
    end do ! loop over blocks
    deallocate(block1%zcontent,block2%zcontent,block3%dcontent)
    deallocate(inter, inter2)
    
  end subroutine generate_chi_block
  !-----------------------------------------------------------------------------
  subroutine generate_Bvector_b(inbl,omega,inputparam,object,p_file,c_file, pol)
    use mod_io, only: io, input
    use mod_matmul, only: matprod
    use hdf5, only: hid_t
    implicit none
    type(block1d), intent(inout) :: inbl
    real(8), intent(in) :: omega, pol(3)
    type(input), intent(in) :: inputparam
    type(io), intent(in) :: object
    integer(hid_t), intent(in) :: p_file, c_file
    !internal variables
    type(block2d) :: bl2d_
    type(block1d) :: bl1d_
    complex(8) :: alpha, beta
    integer(4) :: k
    alpha=1.0d0
    beta=1.0d0
    ! allocate output
    if (allocated(inbl%zcontent)) deallocate(inbl%zcontent)
    allocate(inbl%zcontent(inbl%blocksize))
    inbl%zcontent(:)=0.0d0
    ! loop over all blocks
    do k=1, inbl%nblocks
      ! create blocks
      bl2d_%nblocks=inbl%nblocks
      bl2d_%blocksize=(/inbl%blocksize, inbl%blocksize/)
      bl2d_%global=(/inbl%nblocks*inbl%blocksize, inbl%nblocks*inbl%blocksize/)
      bl2d_%nk=inbl%nk
      bl2d_%il=inbl%il
      bl2d_%iu=inbl%iu
      bl2d_%k1l=inbl%kl
      bl2d_%k1u=inbl%ku
      bl2d_%jl=(k-1)*inbl%blocksize+1
      bl2d_%ju=k*inbl%blocksize
      bl2d_%k2l=(k-1)*bl2d_%nk
      bl2d_%k2u=k*bl2d_%nk
      bl2d_%offset(1)=inbl%offset
      bl2d_%offset(2)=(k-1)*inbl%blocksize
      bl2d_%id(1)=inbl%id
      bl2d_%id(2)=k

      bl1d_%nblocks=inbl%nblocks
      bl1d_%blocksize=inbl%blocksize
      bl1d_%global=inbl%nblocks*inbl%blocksize
      bl1d_%nk=inbl%nk
      bl1d_%il=(k-1)*inbl%blocksize+1
      bl1d_%iu=k*inbl%blocksize
      bl1d_%kl=(k-1)*bl1d_%nk+1
      bl1d_%ku=k*bl1d_%nk
      bl1d_%offset=(k-1)*inbl%blocksize
      bl1d_%id=k
      ! generate core chi block
      call generate_chi_block(bl2d_,omega,inputparam,c_file)
      ! generate t vector block
      call generate_tblock(bl1d_,object%koulims,object%smap,object%ismap,pol,p_file)
      ! generate B vector block
      call zgemm('n','n',inbl%blocksize,1,inbl%blocksize,alpha,bl2d_%zcontent,inbl%blocksize,bl1d_%zcontent &
        &  ,inbl%blocksize,beta,inbl%zcontent,inbl%blocksize)
      deallocate(bl2d_%zcontent,bl1d_%zcontent)
    end do ! loop over blocks
    
  end subroutine 
  !-----------------------------------------------------------------------------
  subroutine generate_Avector_b(inbl,omega,inputparam,core,optical,p_file,c_file,dataset,pol)
    use hdf5, only: hid_t
    use mod_io, only: io, input
    use mod_matmul, only: matprod
    use modmpi, only: mpiglobal
    implicit none
    type(block1d), intent(inout) :: inbl
    real(8), intent(in) :: omega, pol(3)
    type(io), intent(in) :: core, optical
    type(input), intent(in) :: inputparam
    integer(hid_t) :: p_file, c_file, dataset
    ! internal variables
    type(block1d) :: vecB_b
    type(block3d)  :: tprime_b, matB_b,  matA_b
    integer(4) :: interdim(2), nkmax, nu, no, nk, k
    integer(4) :: id_, blsz_
    integer(4), allocatable :: koulims_comb(:,:)
    complex(8), allocatable :: matB_(:,:,:)
    ! generate a combined map for t'
    interdim=shape(core%koulims)
    nkmax=interdim(2)
    allocate(koulims_comb(4,nkmax))
    koulims_comb(1,:)=optical%koulims(3,:)
    koulims_comb(2,:)=optical%koulims(4,:)
    koulims_comb(3,:)=core%koulims(3,:)
    koulims_comb(4,:)=core%koulims(4,:)
    
    ! note that B vector block and A vector block can have different size
    ! set temporary id and size
    id_=inbl%id
    blsz_=core%no*core%nu*inbl%nk
    ! set up block for B vector
    vecB_b%nblocks=inbl%nblocks
    vecB_b%blocksize=blsz_
    vecB_b%global=blsz_*inbl%nblocks
    vecB_b%nk=inbl%nk
    vecB_b%il=(id_-1)*blsz_+1
    vecB_b%iu=id_*blsz_
    vecB_b%kl=inbl%kl
    vecB_b%ku=inbl%ku
    vecB_b%offset=(id_-1)*blsz_
    vecB_b%id=id_
    ! set up block for t' matrix
    tprime_b%nblocks=inbl%nblocks
    tprime_b%nk=inbl%nk
    tprime_b%kl=inbl%kl
    tprime_b%ku=inbl%ku
    tprime_b%id=inbl%id
    ! set up block for A matrix
    matA_b%nblocks=inbl%nblocks
    matA_b%nk=inbl%nk
    matA_b%kl=inbl%kl
    matA_b%ku=inbl%ku
    matA_b%id=inbl%id
    ! generate block of B vector
    call generate_Bvector_b(vecB_b,omega,inputparam,core,p_file,c_file, pol)
    call put_block1d(vecB_b,dataset)
    ! generate block of B matrix
    call transform2matrix_b(core%koulims,core%smap,vecB_b,matB_b)
    ! generate block of tprimegenerate_Bvector_b
    call generate_tprime_block(tprime_b,pol,koulims_comb,p_file)
    ! allocate block of A matrix
    nu=optical%koulims(2,1)-optical%koulims(1,1)+1
    no=optical%koulims(4,1)-optical%koulims(3,1)+1

    nk=inbl%nk
    if (allocated(matA_b%zcontent)) deallocate(matA_b%zcontent)
    allocate(matA_b%zcontent(nu,no,nk))
    matA_b%nblocks=inbl%nblocks
    matA_b%nk=inbl%nk
    matA_b%blocksize=(/nu, no, nk/)
    matA_b%kl=inbl%kl
    matA_b%ku=inbl%ku
    matA_b%id=inbl%id
    ! in case optical & core calculations have different numbers of empty states, the matrices have to be adjusted
    if (allocated(matB_)) deallocate(matB_)
    allocate(matB_(nu, core%no, inbl%nk))
    if (nu .gt. core%nu) then
      ! more empty states in optical calculation than in core one
      matB_(:,:,:)=0.0d0
      matB_(1:core%nu,:,:)=matB_b%zcontent
    elseif (nu .lt. core%nu) then
      ! less empty states in optical calculation than in core one
      matB_(:,:,:)=matB_b%zcontent(1:nu,:,:)
    else
      matB_=matB_b%zcontent
    end if
    do k=1, nk
      ! generate block of A matrix
      call matprod(matB_(:,:,k),tprime_b%zcontent(:,:,k),matA_b%zcontent(:,:,k))
    end do
    ! generate block of A vector
    call transform2vector_b(optical%koulims,optical%smap,matA_b,inbl)
    deallocate(koulims_comb) 
    deallocate(vecB_b%zcontent,matB_b%zcontent,matA_b%zcontent)
    deallocate(matB_)
  end subroutine
  
  !-----------------------------------------------------------------------------
  subroutine gen_prod_b(inbl, inputparam, core, optical, core_id, pmat_id)
    use mod_io, only: io, input
    use hdf5, only: hid_t
    implicit none
    type(block2d), intent(inout) :: inbl
    type(input), intent(in) :: inputparam
    type(io), intent(in) :: core, optical
    integer(hid_t), intent(in) :: core_id, pmat_id
    !local variables
    type(block2d) :: eigvec
    type(block3d) :: tprime_b
    integer(4), allocatable :: koulims_comb(:,:)
    complex(8), allocatable :: eigvec_matrix(:,:,:,:)
    complex(8) :: alpha, beta
    complex(8), allocatable :: inter(:,:), inter2(:,:), inter3(:,:)
    complex(8), allocatable :: prod_prime(:,:,:,:), prod_matrix(:,:,:,:)
    integer(4) :: interdim(2), nkmax, id_(2), blsz_, nk_, global_
    integer(4) :: lambda, ik

    ! generate combined koulims index range
    interdim=shape(core%koulims)
    nkmax=interdim(2)
    allocate(koulims_comb(4,nkmax))
    koulims_comb(1,:)=optical%koulims(3,:)
    koulims_comb(2,:)=optical%koulims(4,:)
    koulims_comb(3,:)=core%koulims(3,:)
    koulims_comb(4,:)=core%koulims(4,:)
    
    ! set temporary id and size
    ! core eigenvector block does not have
    ! the same size as the inblock
    id_=inbl%id
    blsz_=core%no*core%nu*inbl%nk
    global_=core%no*core%nu*nkmax
    ! set up block for core eigenstates
    eigvec%nblocks=inbl%nblocks
    eigvec%blocksize=(/ blsz_, inbl%blocksize(2) /)
    eigvec%global=global_
    eigvec%nk=inbl%nk
    eigvec%il=(inbl%id(1)-1)*blsz_+1
    eigvec%iu=inbl%id(1)*blsz_
    eigvec%jl=inbl%jl
    eigvec%ju=inbl%ju
    eigvec%k1l=inbl%k1l
    eigvec%k1u=inbl%k1u
    eigvec%offset=(/ eigvec%il-1, inbl%offset(2) /)
    eigvec%id=inbl%id
    ! set up block for t' matrix
    tprime_b%nblocks=inbl%nblocks
    tprime_b%nk=inbl%nk
    tprime_b%kl=inbl%k1l
    tprime_b%ku=inbl%k1u
    tprime_b%id=inbl%id(1)
    ! get block of core eigenstates
    if (inputparam%ip_c) then
      call get_eigvecsIP_b(eigvec, core)
    else
      call get_eigvecs_b(eigvec, core_id)
    end if
    
    ! generate block of B matrix
    call transform_matrix2matrix(core%koulims,core%smap,eigvec,eigvec_matrix)
    
    ! generate block of tprime
    call generate_tprime_block(tprime_b,inputparam%pol,koulims_comb,pmat_id)
    
    nk_=inbl%nk
    ! in case optical & core calculations have different numbers of empty states, the matrices have to be adjusted
    if (allocated(prod_prime)) deallocate(prod_prime)
    allocate(prod_prime(core%nu, optical%no, inbl%nk, inbl%blocksize(2)))
    prod_prime(:,:,:,:)=cmplx(0.0d0, 0.0d0)
    allocate(prod_matrix(optical%nu, optical%no, inbl%nk, inbl%blocksize(2)))
    prod_matrix(:,:,:,:)=cmplx(0.0d0, 0.0d0)

    ! loop over block ofcore excitons
    do lambda=1, inbl%blocksize(2)
      ! loop over k-points in input block
      do ik=1, inbl%nk
        alpha=1.0d0
        beta=1.0d0
        call zgemm('N', 'N', core%nu, optical%no, core%no, alpha, eigvec_matrix(:,:,ik,lambda), & 
          &core%nu, tprime_b%zcontent(:,:,ik), core%no, beta, prod_prime(:,:,ik,lambda), core%nu)
      end do
    end do
    if (optical%nu .gt. core%nu) then
      ! more empty states in optical calculation than in core one
      prod_matrix(:,:,:,:)=0.0d0
      prod_matrix(1:core%nu,:,:,:)=prod_prime(:,:,:,:)
    elseif (optical%nu .lt. core%nu) then
      ! less empty states in optical calculation than in core one
      prod_matrix(:,:,:,:)=prod_prime(1:optical%nu,:,:,:)
    else
      prod_matrix(:,:,:,:)=prod_prime(:,:,:,:)
    end if
    ! generate block of product vector
    call transform_matrix2vector(optical,prod_prime,inbl)
    
    deallocate(koulims_comb) 
    deallocate(prod_prime, prod_matrix, eigvec_matrix)

  end subroutine
  !-----------------------------------------------------------------------------
# ifdef DEBUG
  subroutine generate_oscstr_b(nblocks, blsz, nk, k, omega, inputparam, core, optical, pmat_id, core_id, &
     & optical_id, pol, oscstr_b)
    use mod_io, only: io, input
    use hdf5, only: hid_t
    integer, intent(in) :: nblocks, blsz, nk, k
    real(8), intent(in) :: omega(:)
    type(input), intent(in) :: inputparam
    real(8), intent(in) :: pol(3)
    type(io), intent(in) :: core, optical
    integer(hid_t), intent(in) :: pmat_id, core_id, optical_id
    complex(8), intent(out) :: oscstr_b(:,:)
    ! internal variables
    type(block2d) :: evecs_b
    type(block1d) :: vecA_b
    integer :: k2, w
    complex(8) :: alpha, beta
    
    ! loop over blocks
    do k2=1, nblocks
      ! set up block for eigenvectors
      evecs_b%nblocks=nblocks
      evecs_b%blocksize=(/blsz, blsz/)
      evecs_b%nk=nk
      evecs_b%il=(k2-1)*blsz+1
      evecs_b%iu=k2*blsz
      evecs_b%k1l=(k2-1)*nk+1
      evecs_b%k1u=k2*nk
      evecs_b%jl=(k-1)*blsz+1
      evecs_b%ju=k*blsz
      evecs_b%k2l=(k-1)*nk+1
      evecs_b%k2u=k*nk
      evecs_b%offset(1)=(k2-1)*blsz
      evecs_b%offset(2)=(k-1)*blsz
      evecs_b%id(1)=k2
      evecs_b%id(2)=k
      !set up block for A matrix
      vecA_b%nblocks=nblocks
      vecA_b%blocksize=blsz
      vecA_b%nk=nk
      vecA_b%il=(k2-1)*blsz+1
      vecA_b%iu=k2*blsz
      vecA_b%kl=(k2-1)*nk+1
      vecA_b%ku=k2*nk
      vecA_b%offset=(k2-1)*blsz
      vecA_b%id=k2
      
      ! generate block of eigenvectors
      call get_eigvecs_b(evecs_b,optical_id)
      do w=1, size(omega)
        ! generate block of A vector
        call generate_Avector_b(vecA_b,omega(w),inputparam,core,optical,pmat_id,core_id,pol)
        ! generate block of oscstr
        alpha=1.0d0
        beta=1.0d0
        call zgemm('C','N',blsz,1,blsz,alpha,evecs_b%zcontent,blsz,vecA_b%zcontent,blsz,beta,oscstr_b(:,w),blsz)
      end do ! w
    end do ! k2
  end subroutine
#endif 
!-----------------------------------------------------------------------------
  subroutine put_block1d(in1d,dataset_id)
    use hdf5, only: hid_t
    use mod_phdf5
    implicit none
    type(block1d), intent(inout) :: in1d
    integer(hid_t) :: dataset_id
    ! local variables
    integer, dimension(1) :: dims_, dimsg_, offset_
    ! set dimension & offset
    dims_(1)=in1d%blocksize
    dimsg_(1)=in1d%global
    offset_(1)=in1d%offset
    ! if allocated, write the dcontent
    if (allocated(in1d%dcontent)) then
      call phdf5_write(in1d%dcontent(1),dims_,dimsg_,offset_,dataset_id)
    elseif (allocated(in1d%zcontent)) then
      call phdf5_write(in1d%zcontent(1),dims_,dimsg_,offset_,dataset_id)
    end if
  end subroutine

!-----------------------------------------------------------------------------
  subroutine get_block1d(in1d,dataset_id)
    use hdf5, only: hid_t
    use mod_phdf5
    implicit none
    type(block1d), intent(inout) :: in1d
    integer(hid_t) :: dataset_id
    ! local variables
    integer, dimension(1) :: dims_, dimsg_, offset_
    ! set dimension & offset
    dims_(1)=in1d%blocksize
    dimsg_(1)=in1d%global
    offset_(1)=in1d%offset
    ! if allocated, write the dcontent
    if (allocated(in1d%dcontent)) then
      call phdf5_read(in1d%dcontent(1),dims_,dimsg_,offset_,dataset_id)
    elseif (allocated(in1d%zcontent)) then
      call phdf5_read(in1d%zcontent(1),dims_,dimsg_,offset_,dataset_id)
    end if
  end subroutine

!-----------------------------------------------------------------------------
  subroutine put_block2d(in2d,dataset_id)
    use hdf5, only: hid_t
    use mod_phdf5
    implicit none
    type(block2d), intent(inout) :: in2d
    integer(hid_t) :: dataset_id
    ! local variables
    integer, dimension(2) :: dims_, dimsg_, offset_
    ! set dimension & offset
    dims_=in2d%blocksize
    dimsg_=in2d%global
    offset_=in2d%offset
    ! if allocated, write the dcontent
    if (allocated(in2d%dcontent)) then
      call phdf5_write(in2d%dcontent(1,1),dims_,dimsg_,offset_,dataset_id)
    elseif (allocated(in2d%zcontent)) then
      call phdf5_write(in2d%zcontent(1,1),dims_,dimsg_,offset_,dataset_id)
    end if
  end subroutine

!-----------------------------------------------------------------------------
  subroutine get_block2d(in2d,dataset_id)
    use hdf5, only: hid_t
    use mod_phdf5
    implicit none
    type(block2d), intent(inout) :: in2d
    integer(hid_t) :: dataset_id
    ! local variables
    integer, dimension(2) :: dims_, dimsg_, offset_
    ! set dimension & offset
    dims_=(/ in2d%blocksize(1), in2d%blocksize(2)/)
    dimsg_=(/ in2d%blocksize(1)*in2d%nblocks, in2d%blocksize(2)*in2d%nblocks /) 
    offset_=in2d%offset
    ! if allocated, write the dcontent
    if (allocated(in2d%dcontent)) then
      call phdf5_read(in2d%dcontent(1,1),dims_,dimsg_,offset_,dataset_id)
    elseif (allocated(in2d%zcontent)) then
      call phdf5_read(in2d%zcontent(1,1),dims_,dimsg_,offset_,dataset_id)
    end if
  end subroutine
 
  !-----------------------------------------------------------------------------
  subroutine print_properties1d(block)
    implicit none
    type(block1d), intent(in) :: block

    print *, '  blocksize=', block%blocksize
    print *, '  nblocks=', block%nblocks
    print *, '  il=', block%il
    print *, '  iu=', block%iu
    print *, '  offset=', block%offset
    print *, '  id=', block%id
    if (allocated(block%dcontent)) print *, ' dcontent=', size(block%dcontent)
    if (allocated(block%zcontent)) print *, ' zcontent=', size(block%zcontent)
  end subroutine

  !-----------------------------------------------------------------------------
  subroutine print_properties2d(block)
    implicit none
    type(block2d), intent(in) :: block

    print *, '  blocksize=', block%blocksize
    print *, '  nblocks=', block%nblocks
    print *, '  il=', block%il
    print *, '  iu=', block%iu
    print *, '  jl=', block%jl
    print *, '  ju=', block%ju
    print *, '  offset=', block%offset
    print *, '  id=', block%id
    if (allocated(block%dcontent)) print *, ' dcontent=', shape(block%dcontent)
    if (allocated(block%zcontent)) print *, ' zcontent=', shape(block%zcontent)
  end subroutine

  !-----------------------------------------------------------------------------
  function nofblock(block, globalsize, nblocks)
    implicit none
    integer :: nofblock
    integer, intent(in) :: block, globalsize, nblocks

    nofblock= globalsize / nblocks
    if (mod(globalsize,nblocks)> block-1) nofblock= nofblock+1
  
  end function nofblock

  !-----------------------------------------------------------------------------
  function firstofblock(block, globalsize, nblocks)
    implicit none
    integer :: firstofblock
    integer, intent(in) :: block, globalsize, nblocks
    ! local variable
    integer :: i
    
    firstofblock=1
    do i=1, min(block-1, globalsize-1)
      firstofblock=firstofblock+nofblock(i, globalsize, nblocks)
    end do
    if (globalsize < block) firstofblock=0

    end function firstofblock

  !-----------------------------------------------------------------------------
  function lastofblock(block, globalsize, nblocks)
    integer :: lastofblock
    integer, intent(in) :: block, globalsize, nblocks
    ! local variables
    integer :: i
    lastofblock=0
    do i=1, min(block,globalsize)
      lastofblock=lastofblock+nofblock(i, globalsize, nblocks)
    end do
    if (globalsize < block) lastofblock=-1

    end function lastofblock

end module mod_blocks
