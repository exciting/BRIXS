module mod_blocks
  implicit none

  type :: block1d
    integer(4) :: nblocks     ! number of blocks
    integer(4) :: blocksize   ! length of each block
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
    integer(4) :: blocksize   ! size of each block is blocksize x blocksize
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
    subroutine get_evals_block(inblock1d,fname)
      use mod_hdf5, only: hdf5_read_block
      implicit none
      type(block1d), intent(inout) :: inblock1d
      character(len=1024), intent(in) :: fname
      !local variables
      character(len=1024) :: path, dsetname
      integer, dimension(1) :: dims_, offset_ 
      !allocate output
      if (allocated(inblock1d%dcontent)) deallocate(inblock1d%dcontent)
      allocate(inblock1d%dcontent(inblock1d%blocksize))
    
      path='eigvec-singlet-TDA-BAR-full/0001'
      dsetname='evals'
      ! get data
      dims_(1)=inblock1d%blocksize
      offset_(1)=inblock1d%offset
      call hdf5_read_block(fname,path,dsetname,inblock1d%dcontent(1),dims_,offset_)
  end subroutine
  !-----------------------------------------------------------------------------
  subroutine get_eigvecs2D_b(inblock2d,fname)
    use mod_hdf5, only: hdf5_get_dims, hdf5_read_block, hdf5_read
    implicit none
    type(block2d), intent(inout) :: inblock2d
    character(len=1024), intent(in) :: fname
    !local variables
    complex(8), allocatable :: eigvec_(:)
    integer(4) ::  i, offset_(1), dim_(1)
    character(len=1024) :: path, dsetname
    character(256) :: ci
    ! allocate output
    if (allocated(inblock2d%zcontent)) deallocate(inblock2d%zcontent)
    allocate(inblock2d%zcontent(inblock2d%blocksize,inblock2d%blocksize)) 
    ! get data
    do i=1, inblock2d%blocksize
      write(ci, '(I8.8)') i+inblock2d%offset(2)
      path='eigvec-singlet-TDA-BAR-full/0001/rvec'
      dsetname=trim(adjustl(ci))
      ! Get dimension of eigvec for given lambda
      !call hdf5_get_dims(fname,path,ci,dims_)
      ! Allocate intermediate eigenvector array
      if (allocated(eigvec_)) deallocate(eigvec_)
      allocate(eigvec_(inblock2d%blocksize))
      !allocate(eigvec_(dims_(2)))
      ! Get data
      offset_(1)=inblock2d%offset(1)
      dim_(1)=inblock2d%blocksize
      !call hdf5_read_block(fname,path,dsetname,eigvec_(1),dim_,offset_)
      call hdf5_read_block(fname,path,dsetname,eigvec_(1),dim_,offset_)
      ! Write data to final array
      inblock2d%zcontent(:,i)=eigvec_(:)
    end do
    deallocate(eigvec_)
  end subroutine 
  !-----------------------------------------------------------------------------
  subroutine generate_tblock(inblock1d,koulims,smap,ismap,pol,hdf5_file)
    use mod_hdf5
    implicit none
    type(block1d), intent(inout) :: inblock1d
    integer(4), intent(in) :: koulims(:,:)
    integer(4), intent(in) :: smap(:,:)
    integer(4), intent(in) :: ismap(:,:,:)
    real(8), intent(in) :: pol(3)
    character(*), intent(in) :: hdf5_file
    !internal variables
    integer(4), dimension(2) :: dim_koulims, dim_smap
    integer(4) :: lu, uu, lo, uo, nu, no, nk0
    integer(4) :: k,i,j,dimensions(4)
    character(len=1024) :: path, dsetname, cik
    complex(8), allocatable :: pmat_(:,:,:)
    !complex(8) :: pmat_(3,2,35)
    integer :: stat_var

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
      call hdf5_get_dims(hdf5_file,path,trim(adjustl(dsetname)),dimensions)
      !allocate intermediate transition matrix for each k-point
      if (allocated(pmat_)) deallocate(pmat_)
      allocate(pmat_(dimensions(2),dimensions(3),dimensions(4)), stat=stat_var)
      call hdf5_read(hdf5_file,path,dsetname,pmat_(1,1,1),shape(pmat_))
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
  subroutine generate_tprime_block(in3d,pol,koulims,fname)
    use mod_hdf5
    implicit none
    real(8), intent(in) :: pol(3)
    integer(4), intent(in) :: koulims(:,:)
    character(1024), intent(in) :: fname
    type(block3d), intent(inout) :: in3d
    ! local variables
    integer(4) :: lu, uu, lo, uo, nu, no, nkmax, i,j
    integer(4) :: dimensions(4), inter(2), nk_, k
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
      call hdf5_get_dims(fname,path,trim(adjustl(dsetname)),dimensions)
      !allocate intermediate transition matrix for each k-point
      if (allocated(pmat_)) deallocate(pmat_)
      allocate(pmat_(dimensions(2),dimensions(3),dimensions(4)))
      call hdf5_read(fname,path,dsetname,pmat_(1,1,1),shape(pmat_))
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
  subroutine generate_chi_block(inblock2d,omega,broad,fname)
    use mod_io, only: io
    implicit none
    real(8), intent(in) :: omega
    real(8), intent(in) :: broad
    character(1024), intent(in) :: fname
    type(block2d), intent(inout) :: inblock2d
    ! internal variables
    type(block2d) :: block1, block2
    type(block1d) :: block3
    complex(8), dimension(inblock2d%blocksize,inblock2d%blocksize) :: inter, inter2
    !complex(8) :: inter
    integer(4) :: blsz,i, k
    complex(8) :: alpha, beta
    
    ! shorthand notation for blocksize
    blsz=inblock2d%blocksize
    
    ! allocate chi_b
    if (allocated(inblock2d%zcontent)) deallocate(inblock2d%zcontent)
    allocate(inblock2d%zcontent(inblock2d%blocksize,inblock2d%blocksize))
    inblock2d%zcontent(:,:)=0.0d0
    ! chi_block is calculated as follows:
    ! $\chi^{ij}=\sum_k \hat{X}^{ik}\times\hat{Y}^{ii}\times[\hat{X}^{ki}]$
    do k=1, inblock2d%nblocks
      ! define the blocks for the summation
      ! block1 and block2 are blocks of the eigenvectors
      block1%nblocks=inblock2d%nblocks
      block1%blocksize=blsz
      block1%nk=inblock2d%nk
      block1%il=inblock2d%il
      block1%iu=inblock2d%iu
      block1%k1l=inblock2d%k1l
      block1%k1u=inblock2d%k1u
      block1%jl=(k-1)*blsz+1
      block1%ju=k*blsz
      block1%k2l=(k-1)*block1%nk+1
      block1%k2u=k*block1%nk
      block1%offset(1)=inblock2d%offset(1)
      block1%offset(2)=(k-1)*blsz
      block1%id(1)=inblock2d%id(1)
      block1%id(2)=k

      block2%nblocks=inblock2d%nblocks
      block2%blocksize=blsz
      block2%nk=inblock2d%nk
      block2%il=inblock2d%jl
      block2%iu=inblock2d%ju
      block2%k1l=inblock2d%k2l
      block2%k1u=inblock2d%k2u
      block2%jl=(k-1)*blsz+1
      block2%ju=k*blsz
      block2%k2l=(k-1)*block2%nk+1
      block2%k2u=k*block2%nk
      block2%offset(1)=inblock2d%offset(2)
      block2%offset(2)=(k-1)*blsz
      block2%id(1)=inblock2d%id(2)
      block2%id(2)=k
      ! block3 is a 1D block of the eigenvectors
      block3%nblocks=inblock2d%nblocks
      block3%blocksize=blsz
      block3%nk=inblock2d%nk
      block3%il=(k-1)*blsz+1
      block3%iu=k*blsz
      block3%kl=(k-1)*block3%nk+1
      block3%ku=k*block3%nk
      block3%offset=(k-1)*blsz
      block3%id=k
       
      !get eigenvalues and eigenvectors
      call get_evals_block(block3,fname)
      call get_eigvecs2D_b(block1,fname)
      call get_eigvecs2D_b(block2,fname)
    
      ! allocate array for hermetian conjugate of eigvecs
      inter(:,:)=0.0d0
      do i=1,blsz
        inter(i,i)=-1.0d0/(block3%dcontent(i)*27.211d0-omega+cmplx(0.0d0,broad))
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

      call zgemm('n','c',blsz,blsz,blsz,alpha,inter,blsz,block2%zcontent(1:blsz,1:blsz),&
          &        blsz,beta,inter2,blsz)
      alpha=1.0d0
      beta=1.0d0
      call zgemm('N','N',blsz,blsz,blsz,alpha,block1%zcontent(1:blsz,1:blsz),blsz,inter2,&
        &        blsz,beta,inblock2d%zcontent(1:blsz,1:blsz),blsz)
      !alpha=1.0d0
      !beta=1.0d0
      !call zgemm('N','N',blsz,blsz,blsz,alpha,block1%zcontent(1:blsz,1:blsz),blsz,block2%zcontent,&
      !  &        blsz,beta,inblock2d%zcontent(1:blsz,1:blsz),blsz)
    end do ! loop over blocks
    deallocate(block1%zcontent,block2%zcontent,block3%dcontent)
    
  end subroutine generate_chi_block
  !-----------------------------------------------------------------------------
  subroutine generate_Bvector_b(inbl,omega,broad,object,p_file,c_file, pol)
    use mod_io, only: io
    use mod_matmul, only: matprod
    implicit none
    type(block1d), intent(inout) :: inbl
    real(8), intent(in) :: omega, broad, pol(3)
    type(io) :: object
    character(1024), intent(in) :: p_file, c_file
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
      bl2d_%blocksize=inbl%blocksize
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
      bl1d_%nk=inbl%nk
      bl1d_%il=(k-1)*inbl%blocksize+1
      bl1d_%iu=k*inbl%blocksize
      bl1d_%kl=(k-1)*bl1d_%nk+1
      bl1d_%ku=k*bl1d_%nk
      bl1d_%offset=(k-1)*inbl%blocksize
      bl1d_%id=k
      ! generate core chi block
      !call cpu_time(start)
      call generate_chi_block(bl2d_,omega,broad,c_file)
      !call cpu_time(finish)
      !print *, '        generate_chi_block:', finish-start, 'seconds'
      ! generate t vector block
      !call cpu_time(start)
      call generate_tblock(bl1d_,object%koulims,object%smap,object%ismap,pol,p_file)
      !call cpu_time(finish)
      !print *, '        generate_tblock:', finish-start, 'seconds'
      ! generate B vector block
      call zgemm('n','n',inbl%blocksize,1,inbl%blocksize,alpha,bl2d_%zcontent,inbl%blocksize,bl1d_%zcontent &
        &  ,inbl%blocksize,beta,inbl%zcontent,inbl%blocksize)
      deallocate(bl2d_%zcontent,bl1d_%zcontent)
    end do ! loop over blocks
    
  end subroutine 
  !-----------------------------------------------------------------------------
  subroutine generate_Avector_b(inbl,omega,broad,core,optical,p_file,c_file,pol)
    use mod_io, only: io
    use mod_matmul, only: matprod
    implicit none
    type(block1d), intent(inout) :: inbl
    real(8), intent(in) :: omega, broad, pol(3)
    type(io) :: core, optical
    character(1024), intent(in) :: p_file, c_file
    ! internal variables
    type(block1d) :: vecB_b
    type(block3d)  :: tprime_b, matB_b,  matA_b
    integer(4) :: interdim(2), nkmax, nu, no, nk, k
    integer(4), allocatable :: koulims_comb(:,:)
    ! generate a combined map for t'
    interdim=shape(core%koulims)
    nkmax=interdim(2)
    allocate(koulims_comb(4,nkmax))
    koulims_comb(1,:)=optical%koulims(3,:)
    koulims_comb(2,:)=optical%koulims(4,:)
    koulims_comb(3,:)=core%koulims(3,:)
    koulims_comb(4,:)=core%koulims(4,:)
    
    ! set up block for B vector
    vecB_b%nblocks=inbl%nblocks
    vecB_b%blocksize=core%no*core%nu*inbl%nk
    vecB_b%nk=inbl%nk
    vecB_b%il=inbl%il
    vecB_b%iu=inbl%iu
    vecB_b%kl=inbl%kl
    vecB_b%ku=inbl%ku
    vecB_b%offset=inbl%offset
    vecB_b%id=inbl%id
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
    call generate_Bvector_b(vecB_b,omega,broad,core,p_file,c_file, pol)
    
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

    do k=1, nk
      ! generate block of A matrix
      call matprod(matB_b%zcontent(:,:,k),tprime_b%zcontent(:,:,k),matA_b%zcontent(:,:,k))
    end do
    
      ! generate block of A vector
    call transform2vector_b(optical%koulims,optical%smap,matA_b,inbl)
    deallocate(koulims_comb) 
    deallocate(vecB_b%zcontent,matB_b%zcontent,matA_b%zcontent)
  end subroutine
  
  !-----------------------------------------------------------------------------
  subroutine generate_oscstr_b(nblocks, blsz, nk, k, omega, broad, core, optical, fname_pmat, fname_core, &
     & fname_optical, pol, oscstr_b)
    use mod_io, only: io
    integer, intent(in) :: nblocks, blsz, nk, k
    real(8), intent(in) :: omega(:)
    real(8), intent(in) :: broad, pol(3)
    type(io), intent(in) :: core, optical
    character(1024) :: fname_pmat, fname_core, fname_optical
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
      evecs_b%blocksize=blsz
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
      call get_eigvecs2D_b(evecs_b,fname_optical)
      do w=1, size(omega)
        ! generate block of A vector
        call generate_Avector_b(vecA_b,omega(w),broad,core,optical,fname_pmat,fname_core,pol)
        ! generate block of oscstr
        alpha=1.0d0
        beta=1.0d0
        call zgemm('C','N',blsz,1,blsz,alpha,evecs_b%zcontent,blsz,vecA_b%zcontent,blsz,beta,oscstr_b(:,w),blsz)
      end do ! w
    end do ! k2
  end subroutine
  
  !-----------------------------------------------------------------------------
  subroutine put_block1d(in1d,fname,groupname)
    use mod_hdf5
    implicit none
    type(block1d), intent(in) :: in1d
    character(len=1024), intent(in) :: fname, groupname
    ! local variables
    character(len=1024) :: gname_, id_
    ! create group if needed
    write(id_, '(I4.4)') in1d%id
    ! create group for block
    call hdf5_create_group(fname,groupname,id_)
    gname_=trim(adjustl(groupname))//'/'//trim(adjustl(id_))//"/"
    ! write metadata
    call hdf5_write(fname,gname_,"nblocks",in1d%nblocks)
    call hdf5_write(fname,gname_,"blocksize",in1d%blocksize)
    call hdf5_write(fname,gname_,"nk",in1d%nk)
    call hdf5_write(fname,gname_,"il",in1d%il)
    call hdf5_write(fname,gname_,"iu",in1d%iu)
    call hdf5_write(fname,gname_,"kl",in1d%kl)
    call hdf5_write(fname,gname_,"ku",in1d%ku)
    call hdf5_write(fname,gname_,"offset",in1d%offset)
    ! write content
    if (allocated(in1d%zcontent)) then
      call hdf5_write(fname,gname_,"zcontent",in1d%zcontent(1),shape(in1d%zcontent))
    else
      call hdf5_write(fname,gname_,"dcontent",in1d%dcontent(1),shape(in1d%dcontent))
    end if

  end subroutine
  
  !-----------------------------------------------------------------------------
  subroutine get_block1d(in1d,fname,groupname)
    use mod_hdf5
    implicit none
    type(block1d), intent(inout) :: in1d
    character(len=1024), intent(in) :: fname, groupname
    ! local variables
    character(len=1024) :: gname_, id_, group_
    ! get groupname 
    write(id_, '(I4.4)') in1d%id
    gname_=trim(adjustl(groupname))//'/'//trim(adjustl(id_))//"/"
    
    ! read metadata
    call hdf5_read(fname,gname_,"nblocks",in1d%nblocks)
    call hdf5_read(fname,gname_,"blocksize",in1d%blocksize)
    call hdf5_read(fname,gname_,"nk",in1d%nk)
    call hdf5_read(fname,gname_,"il",in1d%il)
    call hdf5_read(fname,gname_,"iu",in1d%iu)
    call hdf5_read(fname,gname_,"kl",in1d%kl)
    call hdf5_read(fname,gname_,"ku",in1d%ku)
    call hdf5_read(fname,gname_,"offset",in1d%offset)
    ! allocate array if necessary
    if (allocated(in1d%zcontent)) deallocate(in1d%zcontent)
    allocate(in1d%zcontent(in1d%blocksize))
    ! read content
    call hdf5_read(fname,gname_,"zcontent",in1d%zcontent(1),shape(in1d%zcontent))

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
  subroutine get_eigvecs2D_b_quick(inblock2d,fname)
    use hdf5
    implicit none
    type(block2d), intent(inout) :: inblock2d
    character(len=1024), intent(in) :: fname
    !local variables
    !complex(8), allocatable :: eigvec_(:)
    integer(hid_t) :: h5_root_id,dataset_id,group_id, dataspace_id, memspace_id
    integer :: ierr,i,j    
    integer(HSIZE_T), dimension(2) :: h_dims, h_offset
    character*100 errmsg
    real(8), allocatable :: eigvec_(:,:)
    integer(4) ::  dims_(2), offset_(2), ndims
    character(len=1024) :: path, dsetname
    character(256) :: ci
    ! allocate output
    if (allocated(inblock2d%zcontent)) deallocate(inblock2d%zcontent)
    allocate(inblock2d%zcontent(inblock2d%blocksize,inblock2d%blocksize)) 
    !allocate intermediate real matrix
    allocate(eigvec_(2,inblock2d%blocksize))
    eigvec_(:,:)=0.0d0
    ! determine dimensions
    ndims=2
    dims_=(/2,inblock2d%blocksize/)
    offset_=(/0,inblock2d%offset(1)/)
    path='eigvec-singlet-TDA-BAR-full/0001/rvec'
    do i=1,2
      h_dims(i)=dims_(i)
      h_offset(i)=offset_(i)
    enddo
    call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5fopen_f returned ",I6)')ierr
      goto 10
    endif
    call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5gopen_f returned ",I6)')ierr
      goto 10
    endif
    do i=1, inblock2d%blocksize
      write(ci, '(I8.8)') i+inblock2d%offset(2)
      path='eigvec-singlet-TDA-BAR-full/0001/rvec'
      dsetname=trim(adjustl(ci))
      call h5dopen_f(group_id,trim(dsetname),dataset_id,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(hdf5_read_array_d): h5dopen_f returned ",I6)')ierr
        goto 10
      endif
      call h5dget_space_f(dataset_id,dataspace_id, ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(hdf5_read_array_i4): h5dget_space_f returned ",I6)')ierr
        goto 10
      endif
      call h5sselect_hyperslab_f(dataspace_id,H5S_SELECT_SET_F,h_offset,h_dims,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(hdf5_read_array_i4): h5sselect_hyperslab_f returned ",I6)')ierr
        goto 10
      endif
      ! create memory dataspace
      call h5screate_simple_f(ndims,h_dims, memspace_id, ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(hdf5_read_array_i4): h5screate_simple_f returned ",I6)')ierr
        goto 10
      endif
      call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,eigvec_,h_dims,ierr,memspace_id,dataspace_id)
      if (ierr.ne.0) then
        write(errmsg,'("Error(hdf5_read_array_d): h5dread_f returned ",I6)')ierr
        goto 10
      endif
      ! write to output
      do j=1, inblock2d%blocksize
        inblock2d%zcontent(j,i)=cmplx(eigvec_(1,j),eigvec_(2,j))
      end do
      ! close the memory space
      call h5sclose_f(memspace_id, ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(hdf5_read_array_d): h5sclose_f returned ",I6)')ierr
        goto 10
      endif
      ! close the dataspace
      call h5sclose_f(dataspace_id, ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(hdf5_read_array_d): h5sclose_f returned ",I6)')ierr
        goto 10
      endif
      call h5dclose_f(dataset_id,ierr)
      if (ierr.ne.0) then
        write(errmsg,'("Error(hdf5_read_array_d): h5dclose_f returned ",I6)')ierr
        goto 10
      endif
    end do
    call h5gclose_f(group_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5gclose_f returned ",I6)')ierr
      goto 10
    endif
    call h5fclose_f(h5_root_id,ierr)
    if (ierr.ne.0) then
      write(errmsg,'("Error(hdf5_read_array_d): h5fclose_f returned ",I6)')ierr
      goto 10
    endif
    deallocate(eigvec_)
    return
    10 continue
    print *, 'problem with quick read!'
  end subroutine 

end module mod_blocks
