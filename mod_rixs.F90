module mod_rixs
  use mod_hdf5
  use mod_io
  use mod_matmul
  implicit none
  
  public transform2vector
  public transform2matrix
  public generate_t
  public generate_tprime
  public generate_chi
  public generate_Bvector
  contains
  ! Methodenbereich
  !-----------------------------------------------------------------------------
  subroutine transform2matrix(inkoulims,insmap,inarray,outarray)
    implicit none
    integer(4), intent(in) :: inkoulims(:,:)
    integer(4), intent(in) :: insmap(:,:)
    complex(8), intent(in) :: inarray(:)
    complex(8), allocatable, intent(out) :: outarray(:,:,:)

    ! local variables
    integer(4) :: lu,uu,lo,uo,nu,no,nk0,nkmax,i
    integer(4) :: dim_koulims(2),dim_smap(2),hamsiz,i1,i2,i3
  
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
    if (allocated(outarray)) deallocate(outarray)
    allocate(outarray(nu,no,nkmax))

    ! loop over all transitions
    do i=1, hamsiz
      i1=insmap(1,i)-lu+1
      i2=insmap(2,i)-lo+1
      i3=insmap(3,i)-nk0+1
      outarray(i1,i2,i3)=inarray(i)
    end do
  end subroutine 

  !-----------------------------------------------------------------------------
  subroutine transform2vector(inkoulims,insmap,inarray,outarray)
    implicit none
    integer(4), intent(in) :: inkoulims(:,:)
    integer(4), intent(in) :: insmap(:,:)
    complex(8), intent(in) :: inarray(:,:,:)
    complex(8), allocatable, intent(out) :: outarray(:)

    ! local variables
    integer(4) :: lu,uu,lo,uo,nu,no,nk0,i
    integer(4) :: dim_koulims(2),dim_smap(2),nkmax,hamsiz,i1,i2,i3
  
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
    if (allocated(outarray)) deallocate(outarray)
    allocate(outarray(hamsiz))

    ! loop over all transitions
    do i=1, hamsiz
      i1=insmap(1,i)-lu+1
      i2=insmap(2,i)-lo+1
      i3=insmap(3,i)-nk0+1
      outarray(i)=inarray(i1,i2,i3)
    end do
  end subroutine 
  !-----------------------------------------------------------------------------
  subroutine generate_t(koulims,smap,ismap,pol,hdf5_file,pmat)
    implicit none

    integer(4), intent(in) :: koulims(:,:)
    integer(4), intent(in) :: smap(:,:)
    integer(4), intent(in) :: ismap(:,:,:)
    real(8), intent(in) :: pol(3)
    character(*), intent(in) :: hdf5_file
    complex(8), allocatable, intent(out) :: pmat(:)
    !internal variables
    integer(4), dimension(2) :: dim_koulims, dim_smap
    integer(4) :: lu, uu, lo, uo, nu, no, nk0, nkmax, hamsiz
    integer(4) :: k,i,j,dimensions(4)
    character(len=1024) :: path, dsetname, cik
    complex(8), allocatable :: pmat_(:,:,:)
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
    nkmax=dim_koulims(2)
    hamsiz=dim_smap(2)
    ! allocate pmat
    if (allocated(pmat)) deallocate(pmat)
    allocate(pmat(hamsiz))
    dsetname='pmat'
    ! loop over k-points
    do k=1,nkmax
      write(cik, '(I4.4)') k
      !determine size of matrix in hdf5 file
      path='/pmat/'//trim(adjustl(cik))
      call hdf5_get_dims(hdf5_file,path,trim(adjustl(dsetname)),dimensions)
      !allocate intermediate transition matrix for each k-point
      if (allocated(pmat_)) deallocate(pmat_)
      allocate(pmat_(dimensions(2),dimensions(3),dimensions(4)))
      call hdf5_read(hdf5_file,path,dsetname,pmat_(1,1,1),shape(pmat_))
      ! write  transition matrix into file for the states included 
      ! in the BSE calculation
      do i=1, no
        do j=1, nu
          pmat(ismap(j,i,k))=conjg(pmat_(1,i+lo-1,j+lu-1))*pol(1)+&
            & conjg(pmat_(2,i+lo-1,j+lu-1))*pol(2)+&
            & conjg(pmat_(3,i+lo-1,j+lu-1))*pol(3)
        end do
      end do
    end do
      deallocate(pmat_)
  end subroutine 
  !-----------------------------------------------------------------------------
  subroutine generate_tprime(pol,koulims,fname,tprime)
    implicit none
    real(8), intent(in) :: pol(3)
    integer(4), intent(in) :: koulims(:,:)
    character(1024), intent(in) :: fname
    complex(8), allocatable, intent(out) :: tprime(:,:,:)
    ! local variables
    integer(4) :: lu, uu, lo, uo, nu, no, nkmax, i,j
    integer(4) :: k, dimensions(4), inter(2)
    character(1024) :: cik, path, dsetname
    complex(8), allocatable :: pmat_(:,:,:)
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
    allocate(tprime(no,nu,nkmax))
    
    dsetname='pmat'
    
    do k=1,nkmax
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
          tprime(i,j,k)=pol(1)*pmat_(1,i+lo-1,j+lu-1)+&
            &pol(2)*pmat_(2,i+lo-1,j+lu-1)+&
            &pol(3)*pmat_(3,i+lo-1,j+lu-1)
        end do
      end do
    end do
  end subroutine 
  !-----------------------------------------------------------------------------
  subroutine generate_chi(omega,broad,object,fname,chi)
    implicit none
    real(8), intent(in) :: omega(:)
    real(8), intent(in) :: broad
    type(io), intent(inout) :: object
    character(1024), intent(in) :: fname
    complex(8), allocatable, intent(out) :: chi(:,:,:)
    ! internal variables
    integer(4) :: hamsiz,interdim(2),w,i
    complex(8), allocatable :: inter(:,:), eigvecs_H(:,:), inter2(:,:)
    complex(8) :: alpha, beta
    interdim=shape(object%smap)
    hamsiz=interdim(2)

    !get eigenvalues and eigenvectors
    call get_evals(object,fname)
    call get_eigvecs(object,fname)
    ! allocate and fill intermediate matrix
    if (allocated(inter)) deallocate(inter)
    allocate(inter(hamsiz,hamsiz))
    inter(:,:)=0.0d0
    
    ! allocate chi
    if (allocated(chi)) deallocate(chi)
    allocate(chi(hamsiz,hamsiz,size(omega)))
    
    ! pepare matrix-matrix multiplication
    allocate(inter2(hamsiz,hamsiz))
    alpha=1.0d0
    beta=0.0d0

    ! allocate array for hermetian conjugate of eigvecs
    if (allocated(eigvecs_H)) deallocate(eigvecs_H)
    allocate(eigvecs_H(hamsiz,hamsiz))
    eigvecs_H(:,:)=conjg(object%eigvecs(:,:))
    eigvecs_H=transpose(eigvecs_H)
    do w=1,size(omega)
      do i=1,size(object%evals)
        inter(i,i)=-1.0d0/(object%evals(i)*27.211d0-omega(w)+cmplx(0.0d0,broad))
      end do
      !chi(:,:,w)=matmul(object%eigvecs,matmul(inter,eigvecs_H))
      ! generate intermediate matrix
      call zgemm('N','C',hamsiz,hamsiz,hamsiz,alpha,inter,hamsiz,object%eigvecs,&
        &        hamsiz,beta,inter2,hamsiz)
      call zgemm('N','N',hamsiz,hamsiz,hamsiz,alpha,object%eigvecs,hamsiz,inter2,&
        &        hamsiz,beta,chi(:,:,w),hamsiz)
    end do
    deallocate(inter, inter2, eigvecs_H)
  end subroutine   
  !-----------------------------------------------------------------------------
  subroutine generate_Bvector(omega,broad,object,p_file,c_file, pol, B)
    implicit none
    real(8), intent(in) :: omega(:), broad, pol(3)
    type(io) :: object
    character(1024), intent(in) :: p_file, c_file
    complex(8),allocatable, intent(out) :: B(:,:)
    !internal variables
    complex(8), allocatable :: chi(:,:,:), pmat(:)
    integer(4) :: w
    !generate core chi
    call generate_chi(omega,broad,object,c_file,chi)
    ! generate t_vector
    call generate_t(object%koulims,object%smap,object%ismap,pol,p_file,pmat)
    ! allocate output B vector
    allocate(B(size(pmat),size(omega)))
    do w=1, size(omega)
      B(:,w)=matmul(chi(:,:,w),pmat(:))
    end do
    deallocate(chi,pmat)
  end subroutine 
  !-----------------------------------------------------------------------------
  subroutine generate_oscstr(object,A,oscstr)
    implicit none
    type(io), intent(in) :: object
    complex(8), intent(in) :: A(:,:)
    complex(8), intent(out) :: oscstr(:,:)
    ! local variables
    integer :: hamsiz, lambda, omega, dims(2)
    complex(8) :: alpha1, alpha2, beta, alpha
    complex(8), allocatable :: inter(:,:), inter2(:)
    complex(8), allocatable :: oscstr2(:,:)
    ! set paramteres
    hamsiz=object%hamsize
    alpha1=-1.0d0
    alpha2=1.0d0
    beta=0.0d0
    allocate(inter(hamsiz,hamsiz))
    allocate(inter2(hamsiz))
    dims=shape(A)
    ! loop over core frequencies
    do lambda=1,dims(1)
      ! generate eigvecxeigvec matrix
      call zgemm('N','C',hamsiz,hamsiz,1,alpha1,object%eigvecs(:,lambda),hamsiz,object%eigvecs(:,lambda),&
        &        hamsiz,beta,inter,hamsiz)
      ! loop over excitons
      do omega=1,dims(2)
        ! generate eigvecxeigvecxA vector
        call zgemm('N','N',hamsiz,1,hamsiz,alpha2,inter,hamsiz,A(:,omega),hamsiz,beta,inter2,hamsiz)
        ! generate oscillator strength
        call zgemm('C','N',1,1,hamsiz,alpha2,A(:,omega),hamsiz,inter2,hamsiz,beta,oscstr(lambda,omega),hamsiz)
      end do
    end do
    ! test other way of calculating
    allocate(oscstr2(dims(1),dims(2)))
    do omega=1,dims(2)
      alpha=1.0
      beta=0.0
      call zgemm('C','N',hamsiz,1,hamsiz,alpha,object%eigvecs,hamsiz,A(:,omega),hamsiz,beta,oscstr2(:,omega),hamsiz)
    end do
    do omega=1,dims(2)
      do lambda=1,dims(1)
        print *, 'oscstr(', lambda, ',', omega,')=', abs(oscstr(lambda,omega))-abs(oscstr2(lambda,omega))**2
      end do
    end do
    deallocate(inter,inter2)
  end subroutine
end module 
