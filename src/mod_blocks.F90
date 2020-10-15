!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: mod_blocks
!
!> @author
!> Christian Vorwerk, Humboldt Universität zu Berlin.
!
! DESCRIPTION: 
!> This modules contains subroutines that perform all block-wise operations.
!> All subroutines are based on the objects "block1d", "block2d", and "block3d",
!> that contain the blockmatrix of the corresponding dimension together with
!> some metadata.
!
! REVISION HISTORY:
! 09 07 2020 - Added documentation
!------------------------------------------------------------------------------
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
    complex(8), allocatable :: zcontent(:) ! complex-valued content
    real(8), allocatable :: dcontent(:)    ! real-valued content
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
    complex(8), allocatable :: zcontent(:,:) ! complex-valued content
    real(8), allocatable :: dcontent(:,:)    ! real-valued content
  end type block2d

  type :: block3d
    integer(4) :: nblocks                 ! number of blocks
    integer(4) :: nk                      ! number of k-pts / block
    integer(4), dimension(3) :: blocksize ! non-square 3D blocksize
    integer(4) :: kl, ku                  ! absolute k-index within block
    integer(4) :: id                      ! store the id of the block it was generated from
    complex(8), allocatable :: zcontent(:,:,:) ! complex-valued content
  end type block3d
  contains
  ! Methodenbereich

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> Transforms a 2D blockmatrx in transition space into a 4D blockmatrix in 
    !> single-particle space, where the first index of the 2D matrix is 
    !> transformed into the 3D single-particle indices. The second index of the 
    !> 2D matrix is preserved.
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in] inkoulims   
    !> @param[in] insmap
    !> @param[in] inbl2d
    !> @return outbl4d    
    !---------------------------------------------------------------------------  
    subroutine transform_matrix2matrix(inkoulims,insmap,inbl2d,outbl4d)
      implicit none
      integer(4), intent(in) :: inkoulims(:,:)
      integer(4), intent(in) :: insmap(:,:)
      type(block2d), intent(in) :: inbl2d
      complex(8), allocatable, intent(out) :: outbl4d(:,:,:,:)

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
      if (allocated(outbl4d)) deallocate(outbl4d)
      allocate(outbl4d(nu, no, inbl2d%nk,inbl2d%blocksize(2)))
      ! loop over all excitons
      do lambda=1, inbl2d%blocksize(2)
        ! loop over all transitions
        do i=inbl2d%il, inbl2d%iu
          i1=insmap(1,i)-lu+1
          i2=insmap(2,i)-lo+1
          i3=insmap(3,i)-nk0+1
          outbl4d(i1,i2,i3-inbl2d%k1l+1,lambda)=inbl2d%zcontent(i-inbl2d%offset(1),lambda)
        end do !i
      end do !lambda
    end subroutine 

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> Transforms a 4D block matrix, where the first 3 indizes are the
    !> single-particle indices of a transition, into a 2D block matrix in
    !> transition space.
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in] input   
    !> @param[in] in4d
    !> @return outbl2d    
    !---------------------------------------------------------------------------  
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

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Reads block of BSE eigenvalues from HDF5 file. 
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] inblock1d  
    !> @param[in] file_id
    !---------------------------------------------------------------------------  
    subroutine get_evals(inblock1d,file_id)
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
    
    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Reads block of IP eigenvalues, i.e. energy differences between conduction
    !> and valence states, from HDF5 file. The array is stored in the inblock1d. 
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] inblock1d  
    !> @param[in] file_id
    !---------------------------------------------------------------------------  
    subroutine get_evalsIP(inblock1d,file_id)
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

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Read block of BSE eigenvectors from HDF5 file and store information in
    !> 2D block matrix, i.e. in transition-space representation. 
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] inblock2d  
    !> @param[in] file_id
    !---------------------------------------------------------------------------  
    subroutine get_eigvecs(inblock2d,file_id)
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

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Generates IP eigenstates in transition space. The function uses the
    !> ensortidx vector. 
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] inblock2d  
    !> @param[in] io_in
    !---------------------------------------------------------------------------  
    subroutine get_eigvecsIP(inblock2d, io_in)
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

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Generates block of the momentum-matrix elements of the initial
    !> excitation \f$ \langle c \mathbf{k} | \mathbf{e}_1 \cdot \mathbf{p} | \mu \mathbf{k} \f$.
    !> The output is stored in a vector in transition space. 
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] inblock1d  
    !> @param[in] koulims
    !> @param[in] smap
    !> @param[in] ismap
    !> @param[in] pol
    !> @param[in] file_id
    !---------------------------------------------------------------------------  
    subroutine generate_t(inblock1d,koulims,smap,ismap,pol,file_id)
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
        write(cik, '(I8.8)') k
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
    
    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Generates block of the momentum-matrix elements of the de-excitation 
    !> \f$ \langle \mu \mathbf{k} | \mathbf{e}_2 \cdot \mathbf{p} | v \mathbf{k} \f$.
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] in3d  
    !> @param[in] pol
    !> @param[in] koulims
    !> @param[in] file_id
    !---------------------------------------------------------------------------  
    subroutine generate_tprime(in3d,pol,koulims,file_id)
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
        write(cik, '(I8.8)') k
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

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Generates block of the momentum-matrix elements of the de-excitation 
    !> \f$ \langle \mu \mathbf{k} | \mathbf{e}_2 \cdot \mathbf{p} | v \mathbf{k} \f$.
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] in2d  
    !> @param[in] inputparam
    !> @param[in] core
    !> @param[in] optical
    !> @param[in] core_id
    !> @param[in] pmat_id
    !---------------------------------------------------------------------------  
    subroutine generate_product(in2d, inputparam, core, optical, core_id, pmat_id)
      use mod_io, only: io, input
      use hdf5, only: hid_t
      implicit none
      type(block2d), intent(inout) :: in2d
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
      id_=in2d%id
      blsz_=core%no*core%nu*in2d%nk
      global_=core%no*core%nu*nkmax
      ! set up block for core eigenstates
      eigvec%nblocks=in2d%nblocks
      eigvec%blocksize=(/ blsz_, in2d%blocksize(2) /)
      eigvec%global=global_
      eigvec%nk=in2d%nk
      eigvec%il=(in2d%id(1)-1)*blsz_+1
      eigvec%iu=in2d%id(1)*blsz_
      eigvec%jl=in2d%jl
      eigvec%ju=in2d%ju
      eigvec%k1l=in2d%k1l
      eigvec%k1u=in2d%k1u
      eigvec%offset=(/ eigvec%il-1, in2d%offset(2) /)
      eigvec%id=in2d%id
      ! set up block for t' matrix
      tprime_b%nblocks=in2d%nblocks
      tprime_b%nk=in2d%nk
      tprime_b%kl=in2d%k1l
      tprime_b%ku=in2d%k1u
      tprime_b%id=in2d%id(1)
      ! get block of core eigenstates
      if (inputparam%ip_c) then
        call get_eigvecsIP(eigvec, core)
      else
        call get_eigvecs(eigvec, core_id)
      end if
      
      ! generate block of B matrix
      call transform_matrix2matrix(core%koulims,core%smap,eigvec,eigvec_matrix)
      
      ! generate block of tprime
      call generate_tprime(tprime_b,inputparam%pol,koulims_comb,pmat_id)
      
      nk_=in2d%nk
      ! in case optical & core calculations have different numbers of empty states, the matrices have to be adjusted
      if (allocated(prod_prime)) deallocate(prod_prime)
      allocate(prod_prime(core%nu, optical%no, in2d%nk, in2d%blocksize(2)))
      prod_prime(:,:,:,:)=cmplx(0.0d0, 0.0d0)
      allocate(prod_matrix(optical%nu, optical%no, in2d%nk, in2d%blocksize(2)))
      prod_matrix(:,:,:,:)=cmplx(0.0d0, 0.0d0)

      ! loop over block ofcore excitons
      do lambda=1, in2d%blocksize(2)
        ! loop over k-points in input block
        do ik=1, in2d%nk
          alpha=1.0d0
          beta=0.0d0
          ! complex conjugate of eigenstates needed
          eigvec_matrix(:,:,ik,lambda)=conjg(eigvec_matrix(:,:,ik,lambda))
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
      call transform_matrix2vector(optical,prod_matrix,in2d)
      
      deallocate(koulims_comb) 
      deallocate(prod_prime, prod_matrix, eigvec_matrix)
    end subroutine
  
    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Writes 1D block matrix into HDF5 file. If parallel HDF5 is used, the
    !> blocks can be written parallely. 
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] in1d  
    !> @param[in] dataset_id 
    !---------------------------------------------------------------------------  
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

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Reads 1D block matrix into HDF5 file. If parallel HDF5 is used, the
    !> blocks can be read parallely. 
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] in1d  
    !> @param[in] dataset_id 
    !---------------------------------------------------------------------------  
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

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Writes 2D block matrix into HDF5 file. If parallel HDF5 is used, the
    !> blocks can be written parallely. 
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] in2d  
    !> @param[in] dataset_id 
    !---------------------------------------------------------------------------  
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

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Reads 2D block matrix into HDF5 file. If parallel HDF5 is used, the
    !> blocks can be read parallely. 
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[inout] in2d  
    !> @param[in] dataset_id 
    !---------------------------------------------------------------------------  
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
 
    !--------------------------------------------------------------------------
    function nofblock(block, globalsize, nblocks)
      implicit none
      integer :: nofblock
      integer, intent(in) :: block, globalsize, nblocks

      nofblock= globalsize / nblocks
      if (mod(globalsize,nblocks)> block-1) nofblock= nofblock+1
    
    end function nofblock

    !--------------------------------------------------------------------------
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

    !--------------------------------------------------------------------------
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
