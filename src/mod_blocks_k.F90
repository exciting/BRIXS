module mod_blocks_k
  use mod_blocks
  use mod_mpi
  implicit none
  
  contains

    !-----------------------------------------------------------------------------
    subroutine transform_matrix2matrix_k(inkoulims, ik,insmap, inbl2d, out3d)
      implicit none
      integer(4), intent(in) :: inkoulims(:,:)
      integer(4), intent(in) :: insmap(:,:), ik
      type(block2d), intent(inout) :: inbl2d
      complex(8), allocatable, intent(out) :: out3d(:,:,:)

      ! local variables
      complex(8), allocatable :: inter4d(:,:,:,:)
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
      ! set k-ranges that are required for transform_matrix2matrix
      inbl2d%k1l=ik
      inbl2d%k1u=ik
      ! allocate the output matrix
      if (allocated(out3d)) deallocate(out3d)
      allocate(out3d(nu, no, inbl2d%blocksize(2)))
      
      ! generate intermediate 4D array
      call transform_matrix2matrix(inkoulims, insmap, inbl2d, inter4d)

      ! generate required slice
      out3d(:,:,:)= inter4d(:, :, 1, :)

      ! deallocate intermediate block
      deallocate(inter4d)
    end subroutine 
    
    !-----------------------------------------------------------------------------
    subroutine transform_matrix2vector_k(input, ik, in3d, outbl2d)
      use mod_io, only: io
      implicit none
      type(io), intent(in) :: input
      integer, intent(in) :: ik
      complex(8), allocatable, intent(in) :: in3d(:,:,:)
      type(block2d), intent(inout) :: outbl2d

      ! local variables
      integer(4) :: inshape(3), koulimsshape(2), nkmax
      complex(8), allocatable :: inter4d(:,:,:,:)
  
      ! get shape of in3d and koulins
      inshape = shape(in3d)
      koulimsshape = shape(input%koulims)
      nkmax=koulimsshape(2)
      ! allocate intermediate 4d block
      allocate(inter4d(inshape(1), inshape(2), 1, inshape(3)))
      ! fill the date from in3d to inter4d
      inter4d(:,:,:,:) = 0.0d0
      inter4d(:,:,1,:)= in3d(:,:,:)
      
      ! set k-ranges that are required for transform_matrix2vector
      outbl2d%k1l = ik
      outbl2d%k1u = ik
      !generate the block output
      call transform_matrix2vector(input, inter4d, outbl2d)

      deallocate(inter4d)
    end subroutine 
  
    !-----------------------------------------------------------------------------
    subroutine generate_t_k(inblock1d,ik,koulims,smap,ismap,pol,file_id)
      use hdf5, only: hid_t
      use mod_phdf5
      implicit none
      type(block1d), intent(inout) :: inblock1d
      integer(4), intent(in) :: ik, koulims(:,:)
      integer(4), intent(in) :: smap(:,:)
      integer(4), intent(in) :: ismap(:,:,:)
      real(8), intent(in) :: pol(3)
      integer(hid_t), intent(in) :: file_id
      !internal variables
      integer(4), dimension(2) :: dim_koulims, dim_smap
      integer(4) :: lu, uu, lo, uo, nu, no, nk0
      integer(4) :: i,j,dimensions(4), dimsg_(3), offset_(3)
      character(len=1024) :: path, dsetname, cik
      complex(8), allocatable :: pmat_(:,:,:)
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

      write(cik, '(I8.8)') ik
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
          inblock1d%zcontent(ismap(j,i,ik)-inblock1d%offset)=conjg(pmat_(1,i+lo-1,j+lu-1))*pol(1)+&
            & conjg(pmat_(2,i+lo-1,j+lu-1))*pol(2)+&
            & conjg(pmat_(3,i+lo-1,j+lu-1))*pol(3)
        end do
      end do
      deallocate(pmat_)
    end subroutine 

    !-----------------------------------------------------------------------------
    subroutine generate_tprime_k(in2d,ik,pol,koulims,file_id)
      use mod_phdf5
      use hdf5, only: hid_t
      implicit none
      real(8), intent(in) :: pol(3)
      integer(4), intent(in) :: ik, koulims(:)
      integer(hid_t), intent(in) :: file_id
      type(block2d), intent(inout) :: in2d
      ! local variables
      integer(4) :: lu, uu, lo, uo, nu, no, i,j
      integer(4) :: dimensions(4)
      integer :: dimsg_(3), offset_(3)
      integer(hid_t) :: dataset_id
      character(1024) :: cik, path, dsetname
      complex(8), allocatable :: pmat_(:,:,:)
      !complex(8) :: pmat_(3,2,35)
      !determine sizes
      lu=koulims(1)
      uu=koulims(2)
      lo=koulims(3)
      uo=koulims(4)
      nu=uu-lu+1
      no=uo-lo+1
      ! allocate output array
      if (allocated(in2d%zcontent)) deallocate(in2d%zcontent)
      allocate(in2d%zcontent(no,nu))
      in2d%blocksize=(/no, nu/)
      dsetname='pmat'
        
      write(cik, '(I8.8)') ik
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
          in2d%zcontent(i,j)=pol(1)*pmat_(1,i+lo-1,j+lu-1)+&
            & pol(2)*pmat_(2,i+lo-1,j+lu-1)+&
            & pol(3)*pmat_(3,i+lo-1,j+lu-1)
        end do
      end do
      deallocate(pmat_)
    end subroutine 

    !-----------------------------------------------------------------------------
    subroutine gen_t1_k(t1_b, k, core, core_id, pmat_id, inputparam)
      use  hdf5, only: hid_t
      use mod_io, only: io, input

      implicit none

      type(block1d), intent(inout) :: t1_b
      integer, intent(in) :: k
      type(io) :: core
      integer(hid_t) :: core_id, pmat_id
      type(input) :: inputparam
      ! local variables
      type(block1d) :: t_b
      complex(8) :: alpha, beta
      integer :: id_
      type(block2d) :: eigvec_b
      
      id_=t1_b%id

      ! set-up for the blocks of core eigenvectors
      eigvec_b%blocksize=(/ core%globalk, t1_b%blocksize /)
      eigvec_b%global=(/ core%global, inputparam%nstatc /)
      eigvec_b%il=(k-1)*core%globalk+1
      eigvec_b%iu=k*core%globalk
      eigvec_b%jl=t1_b%il
      eigvec_b%ju=t1_b%iu
      eigvec_b%offset(1)=(k-1)*core%globalk
      eigvec_b%offset(2)=t1_b%offset
      eigvec_b%id=(/ k, id_ /)
      ! set-up for the blocks of t
      t_b%blocksize=core%globalk
      t_b%global=core%global
      t_b%il=(k-1)*core%globalk+1
      t_b%iu=k*core%globalk
      t_b%offset=(k-1)*core%globalk
      t_b%id=k
      
      ! generate block of X
      call get_eigvecs(eigvec_b, core_id)
        ! generate block of t
      call generate_t_k(t_b, k, core%koulims, core%smap, core%ismap, &
       & inputparam%pol, pmat_id)
      ! matrix-vector multiplication
      alpha=cmplx(1.0d0, 0.0d0)
      beta=cmplx(0.0d0, 0.0d0)
      call zgemm('T','N', eigvec_b%blocksize(2), 1, eigvec_b%blocksize(1), &
       & alpha, eigvec_b%zcontent, eigvec_b%blocksize(1), t_b%zcontent,    &
       & t_b%blocksize, beta, t1_b%zcontent, t1_b%blocksize)
    end subroutine

    !-----------------------------------------------------------------------------
    subroutine gen_t2_k(outblock2d, ik, tprime_b, core, optical, core_id, &
        optical_id, inputparam)
      use  hdf5, only: hid_t
      use mod_io, only: io, input
      
      implicit none
      
      type(block2d), intent(inout) :: outblock2d
      integer(4), intent(in) :: ik
      type(io), intent(in) :: core, optical
      integer(hid_t), intent(in) :: core_id, optical_id
      type(block2d), intent(in) :: tprime_b
      type(input), intent(in) :: inputparam
      !local variables
      type(block2d) :: prod_, eigvec_
      integer(4) :: nblocks_, blk_, blk2_
      complex(8) :: alpha, beta
      
      nblocks_=inputparam%nblocks
      blk_=outblock2d%id(1)
      blk2_=outblock2d%id(2)
      ! prepare output array
      if (allocated(outblock2d%zcontent)) deallocate(outblock2d%zcontent)
      allocate(outblock2d%zcontent(outblock2d%blocksize(1), outblock2d%blocksize(2)))
      outblock2d%zcontent(:,:)=cmplx(0.0d0, 0.0d0)
          
      ! set-up block for prod vector
      prod_%nblocks=nblocks_
      prod_%blocksize=(/ optical%globalk, nofblock(blk2_, inputparam%nstatc, nblocks_)  /)
      prod_%global=(/ optical%global, inputparam%nstatc /)
      prod_%il=(ik-1)*optical%globalk+1
      prod_%iu=ik*optical%globalk
      prod_%jl=firstofblock(blk2_, inputparam%nstatc, nblocks_)
      prod_%ju=lastofblock(blk2_, inputparam%nstatc, nblocks_)
      prod_%offset=(/ prod_%il-1, prod_%jl-1 /)
      prod_%id=(/ ik, blk2_ /)

      ! set up block of optical eigenvectors
      eigvec_%nblocks=nblocks_
      eigvec_%blocksize=(/ optical%globalk, outblock2d%blocksize(1) /)
      eigvec_%global=(/ optical%global, inputparam%nstato /)
      eigvec_%il=(ik-1)*optical%globalk+1
      eigvec_%iu=ik*optical%globalk
      eigvec_%jl=outblock2d%il
      eigvec_%ju=outblock2d%iu
      eigvec_%offset(1)=eigvec_%il-1
      eigvec_%offset(2)=outblock2d%offset(1)
      eigvec_%id=(/ ik, blk_ /)
      ! generate block of eigenvectors
      call get_eigvecs(eigvec_, optical_id)
      ! generate block of intermediate product
      call gen_prod_k(prod_, tprime_b, ik, core, optical, core_id)

      !matrix-matrix multiplication
      alpha=cmplx(1.0d0, 0.0d0)
      beta=cmplx(1.0d0, 0.0d0)
      call zgemm('T', 'N', eigvec_%blocksize(2), prod_%blocksize(2), eigvec_%blocksize(1), alpha, eigvec_%zcontent, &
      & eigvec_%blocksize(1), prod_%zcontent, prod_%blocksize(1), beta, outblock2d%zcontent, outblock2d%blocksize(1))
    end subroutine
  
    !-----------------------------------------------------------------------------
    subroutine gen_prod_k(inbl, tprime_b, ik, core, optical, core_id)
      use mod_io, only: io, input
      use hdf5, only: hid_t
      implicit none
      type(block2d), intent(inout) :: inbl
      type(block2d), intent(in) :: tprime_b
      integer, intent(in) :: ik
      type(io), intent(in) :: core, optical
      integer(hid_t), intent(in) :: core_id
      !local variables
      type(block2d) :: eigvec
      complex(8), allocatable :: eigvec_matrix(:,:,:)
      complex(8) :: alpha, beta
      complex(8), allocatable :: prod_prime(:,:,:), prod_matrix(:,:,:)
      integer(4) :: id_(2), blsz_
      integer(4) :: lambda

      ! set temporary id and size
      ! core eigenvector block does not have
      ! the same size as the inblock
      id_=inbl%id
      blsz_=core%no*core%nu
      ! set up block for core eigenstates
      eigvec%nblocks=inbl%nblocks
      eigvec%blocksize=(/ blsz_, inbl%blocksize(2) /)
      eigvec%global=core%global
      eigvec%il=(inbl%id(1)-1)*blsz_+1
      eigvec%iu=inbl%id(1)*blsz_
      eigvec%jl=inbl%jl
      eigvec%ju=inbl%ju
      eigvec%nk=1
      eigvec%offset=(/ eigvec%il-1, inbl%offset(2) /)
      eigvec%id=inbl%id
      
      ! get block of core eigenstates
      call get_eigvecs(eigvec, core_id)
      
      ! generate block of B matrix
      call transform_matrix2matrix_k(core%koulims,ik,core%smap,eigvec, &
        & eigvec_matrix)
      
      ! in case optical & core calculations have different numbers of empty states, the matrices have to be adjusted
      if (allocated(prod_prime)) deallocate(prod_prime)
      allocate(prod_prime(core%nu, optical%no, inbl%blocksize(2)))
      prod_prime(:,:,:)=cmplx(0.0d0, 0.0d0)
      allocate(prod_matrix(optical%nu, optical%no, inbl%blocksize(2)))
      prod_matrix(:,:,:)=cmplx(0.0d0, 0.0d0)
      ! loop over block ofcore excitons
      do lambda=1, inbl%blocksize(2)
        alpha=1.0d0
        beta=1.0d0
        ! complex conjugate of eigenstates needed
        eigvec_matrix(:,:,lambda)=conjg(eigvec_matrix(:,:,lambda))
        call zgemm('N', 'N', core%nu, optical%no, core%no, alpha, eigvec_matrix(:,:,lambda), & 
          &core%nu, tprime_b%zcontent(:,:), core%no, beta, prod_prime(:,:,lambda), core%nu)
      end do
      if (optical%nu .gt. core%nu) then
        ! more empty states in optical calculation than in core one
        prod_matrix(:,:,:)=0.0d0
        prod_matrix(1:core%nu,:,:)=prod_prime(:,:,:)
      elseif (optical%nu .lt. core%nu) then
        ! less empty states in optical calculation than in core one
        prod_matrix(:,:,:)=prod_prime(1:optical%nu,:,:)
      else
        prod_matrix(:,:,:)=prod_prime(:,:,:)
      end if
      ! generate block of product vector
      call transform_matrix2vector_k(optical, ik, prod_matrix, inbl)
      
      deallocate(prod_prime, prod_matrix, eigvec_matrix)

    end subroutine
end module
