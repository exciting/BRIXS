module mod_blocks
  implicit none

  type :: block1d
    integer(4) :: nblocks     ! number of blocks
    integer(4) :: blocksize   ! length of each block
    integer(4) :: il, iu      ! absolute first & last index within block
    integer(4) :: offset      ! offset
    integer(4) :: id  ! ID (in this case number) of the block
    complex(8), allocatable :: zcontent(:)
    real(8), allocatable :: dcontent(:)
  end type block1d
  
  type :: block2d
    integer(4) :: nblocks     ! number of block
    integer(4) :: blocksize   ! size of each block is blocksize x blocksize
    integer(4) :: il, iu, jl, ju ! absolute index ranges within the block  
    integer(4), dimension(2) :: offset    ! offset
    integer(4), dimension(2) :: id    ! 2-D ID of block
    complex(8), allocatable :: zcontent(:,:)
    real(8), allocatable :: dcontent(:,:)
  end type block2d
  
  contains
  ! Methodenbereich
    !-----------------------------------------------------------------------------
    subroutine transform2matrix_b(inkoulims,insmap,inbl1d,outbl2d)
      implicit none
      integer(4), intent(in) :: inkoulims(:,:)
      integer(4), intent(in) :: insmap(:,:)
      type(block1d), intent(in) :: inbl1d
      type(block2d), intent(out) :: outbl2d

      ! local variables
      integer(4) :: lu,uu,lo,uo,nu,no,nk0,nkmax,i
      integer(4) :: dim_koulims(2),dim_smap(2),hamsiz,i1,i2
  
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
      ! blocksize cannot be defined because it is no longer a square matrix
      ! offset=0 since for each k-point there is only one block
      ! store the original block id in outbl2d%id(1)
      outbl2d%nblocks=inbl1d%nblocks
      outbl2d%il=1
      outbl2d%iu=nu
      outbl2d%il=1
      outbl2d%ju=no
      outbl2d%offset(:)=0
      outbl2d%id(1)=inbl1d%id
      outbl2d%id(2)=0
      
      ! allocate the output matrix
      if (allocated(outbl2d%zcontent)) deallocate(outbl2d%zcontent)
      allocate(outbl2d%zcontent(nu, no))

      ! loop over all transitions
      do i=inbl1d%il, inbl1d%iu
        i1=insmap(1,i)-lu+1
        i2=insmap(2,i)-lo+1
        outbl2d%zcontent(i1,i2)=inbl1d%zcontent(i-inbl1d%offset)
      end do
    end subroutine 
    !-----------------------------------------------------------------------------
    subroutine transform2vector_b(inkoulims,insmap,inbl2d,outbl1d)
      implicit none
      integer(4), intent(in) :: inkoulims(:,:)
      integer(4), intent(in) :: insmap(:,:)
      type(block2d), intent(in) :: inbl2d
      type(block1d), intent(out) :: outbl1d

      ! local variables
      integer(4) :: lu,uu,lo,uo,nu,no,nk0,i
      integer(4) :: dim_koulims(2),dim_smap(2),nkmax,hamsiz,i1,i2
  
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

      !generate the block output
      outbl1d%nblocks=inbl2d%nblocks
      outbl1d%blocksize=nu*no
      outbl1d%il=(inbl2d%id(1)-1)*outbl1d%blocksize+1
      outbl1d%iu=(inbl2d%id(1))*outbl1d%blocksize
      outbl1d%offset=(inbl2d%id(1)-1)*outbl1d%blocksize
      outbl1d%id=inbl2d%id(1)

      ! allocate the output matrix
      if (allocated(outbl1d%zcontent)) deallocate(outbl1d%zcontent)
      allocate(outbl1d%zcontent(outbl1d%blocksize))

      ! loop over all transitions
      do i=outbl1d%il, outbl1d%iu
        i1=insmap(1,i)-lu+1
        i2=insmap(2,i)-lo+1
        outbl1d%zcontent(i-outbl1d%offset)=inbl2d%zcontent(i1,i2)
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
    use mod_hdf5, only: hdf5_get_dims, hdf5_read_block
    implicit none
    type(block2d), intent(inout) :: inblock2d
    character(len=1024), intent(in) :: fname
    !local variables
    complex(8), allocatable :: eigvec_(:)
    integer(4) ::  i, dims_(2), offset_(1)
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
      call hdf5_get_dims(fname,path,ci,dims_)
      ! Allocate intermediate eigenvector array
      if (allocated(eigvec_)) deallocate(eigvec_)
      allocate(eigvec_(inblock2d%blocksize))
      ! Get data
      offset_(1)=inblock2d%offset(1)
      call hdf5_read_block(fname,path,dsetname,eigvec_(1),shape(eigvec_),offset_)
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
      
    ! here we assume that each k-pt is one block
    write(cik, '(I4.4)') inblock1d%id
    k=inblock1d%id
    !determine size of matrix in hdf5 file
    path=trim(adjustl('/pmat/'//trim(adjustl(cik))))
    call hdf5_get_dims(hdf5_file,path,trim(adjustl(dsetname)),dimensions)
    !allocate intermediate transition matrix for each k-point
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
    deallocate(pmat_)
  end subroutine 
!-----------------------------------------------------------------------------
  subroutine generate_tprime_block(inblock2d,pol,koulims,fname)
    use mod_hdf5
    implicit none
    real(8), intent(in) :: pol(3)
    integer(4), intent(in) :: koulims(:,:)
    character(1024), intent(in) :: fname
    type(block2d), intent(inout) :: inblock2d
    ! local variables
    integer(4) :: lu, uu, lo, uo, nu, no, nkmax, i,j
    integer(4) :: dimensions(4), inter(2)
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
    ! allocate output array
    if (allocated(inblock2d%zcontent)) deallocate(inblock2d%zcontent)
    allocate(inblock2d%zcontent(no,nu))
    
    dsetname='pmat'
    write(cik, '(I4.4)') inblock2d%id(1)
    !determine size of matrix in hdf5 file
    path='/pmat/'//trim(adjustl(cik))
    call hdf5_get_dims(fname,path,trim(adjustl(dsetname)),dimensions)
    !allocate intermediate transition matrix for each k-point
    if (allocated(pmat_)) deallocate(pmat_)
    allocate(pmat_(dimensions(2),dimensions(3),dimensions(4)))
    call hdf5_read(fname,path,dsetname,pmat_(1,1,1),shape(pmat_))
    do i=1,no
      do j=1,nu
        inblock2d%zcontent(i,j)=pol(1)*pmat_(1,i+lo-1,j+lu-1)+&
          &pol(2)*pmat_(2,i+lo-1,j+lu-1)+&
          &pol(3)*pmat_(3,i+lo-1,j+lu-1)
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
      block1%il=inblock2d%il
      block1%iu=inblock2d%iu
      block1%jl=(k-1)*blsz+1
      block1%ju=k*blsz
      block1%offset(1)=inblock2d%offset(1)
      block1%offset(2)=(k-1)*blsz
      block1%id(1)=inblock2d%id(1)
      block1%id(2)=k

      block2%nblocks=inblock2d%nblocks
      block2%blocksize=blsz
      block2%il=inblock2d%jl
      block2%iu=inblock2d%ju
      block2%jl=(k-1)*blsz+1
      block2%ju=k*blsz
      block2%offset(1)=inblock2d%offset(2)
      block2%offset(2)=(k-1)*blsz
      block2%id(1)=inblock2d%id(2)
      block2%id(2)=k
      ! block3 is a 1D block of the eigenvectors
      block3%nblocks=inblock2d%nblocks
      block3%blocksize=blsz
      block3%il=(k-1)*blsz+1
      block3%iu=k*blsz
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
      alpha=1.0d0
      beta=0.0d0

      call zgemm('n','c',blsz,blsz,blsz,alpha,inter,blsz,block2%zcontent(1:blsz,1:blsz),&
          &        blsz,beta,inter2,blsz)
      alpha=1.0d0
      beta=1.0d0
      call zgemm('N','N',blsz,blsz,blsz,alpha,block1%zcontent(1:blsz,1:blsz),blsz,inter2,&
        &        blsz,beta,inblock2d%zcontent(1:blsz,1:blsz),blsz)
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
      bl2d_%il=inbl%il
      bl2d_%iu=inbl%iu
      bl2d_%jl=(k-1)*inbl%blocksize+1
      bl2d_%ju=k*inbl%blocksize
      bl2d_%offset(1)=inbl%offset
      bl2d_%offset(2)=(k-1)*inbl%blocksize
      bl2d_%id(1)=inbl%id
      bl2d_%id(2)=k

      bl1d_%nblocks=inbl%nblocks
      bl1d_%blocksize=inbl%blocksize
      bl1d_%il=(k-1)*inbl%blocksize+1
      bl1d_%iu=k*inbl%blocksize
      bl1d_%offset=(k-1)*inbl%blocksize
      bl1d_%id=k
      ! generate core chi block
      call generate_chi_block(bl2d_,omega,broad,c_file)
      ! generate t vector block
      call generate_tblock(bl1d_,object%koulims,object%smap,object%ismap,pol,p_file)
      ! generate B vector block

      call zgemm('n','n',inbl%blocksize,1,inbl%blocksize,alpha,bl2d_%zcontent,inbl%blocksize,bl1d_%zcontent &
        &  ,inbl%blocksize,beta,inbl%zcontent,inbl%blocksize)
      deallocate(bl2d_%zcontent,bl1d_%zcontent)
    end do ! loop over blocks
    if (allocated(bl2d_%zcontent)) deallocate(bl2d_%zcontent)
    if (allocated(bl1d_%zcontent)) deallocate(bl1d_%zcontent)
    
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
    type(block2d) :: matB_b, tprime_b, matA_b
    integer(4) :: interdim(2), nkmax, nu, no
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
    vecB_b%blocksize=core%no*core%nu
    vecB_b%il=inbl%il
    vecB_b%iu=inbl%iu
    vecB_b%offset=inbl%offset
    vecB_b%id=inbl%id
    ! set up block for t' matrix
    tprime_b%nblocks=inbl%nblocks
    tprime_b%id(1)=inbl%id
    tprime_b%id(2)=0
    ! set up block for A matrix
    matA_b%nblocks=inbl%nblocks
    matA_b%id(1)=inbl%id
    matA_b%id(2)=0
    ! generate block of B vector
    call generate_Bvector_b(vecB_b,omega,broad,core,p_file,c_file, pol)
    ! generate block of B matrix
    call transform2matrix_b(core%koulims,core%smap,vecB_b,matB_b)
    ! generate block of tprimegenerate_Bvector_b
    call generate_tprime_block(tprime_b,pol,koulims_comb,p_file)
    ! allocate block of A matrix
    nu=optical%koulims(2,1)-optical%koulims(1,1)+1
    no=optical%koulims(4,1)-optical%koulims(3,1)+1
    if (allocated(matA_b%zcontent)) deallocate(matA_b%zcontent)
    allocate(matA_b%zcontent(nu,no))
    ! generate block of A matrix
    call matprod(matB_b%zcontent,tprime_b%zcontent,matA_b%zcontent)
    ! generate block of A vector
    call transform2vector_b(optical%koulims,optical%smap,matA_b,inbl)
    deallocate(koulims_comb) 
    deallocate(vecB_b%zcontent,matB_b%zcontent,matA_b%zcontent)
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


   
end module mod_blocks
