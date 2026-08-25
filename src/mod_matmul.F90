!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: mod_matmul
!
!> @author
!> Christian Vorwerk, Humboldt Universität zu Berlin.
!
! DESCRIPTION: 
! This module contains subroutines for convenient matrix-matrix multiplications.
!
! REVISION HISTORY:
! 09 07 2020 - Added documentation
!------------------------------------------------------------------------------
module mod_matmul
  implicit none
  
  interface matprod
    module procedure matprod_matmat, &
        &            matprod_matvec, &
        &            matprod_vecvec
  end interface

  contains
  ! Methodenbereich
  !-----------------------------------------------------------------------------
  
    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Matrix-matrix multiplication of two 2D matrices.
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in]  a   
    !> @param[in]  b 
    !> @param[out] c    
    !---------------------------------------------------------------------------  
    subroutine matprod_matmat(a,b,c)
      implicit none
      complex(8), intent(in) :: a(:,:), b(:,:)
      complex(8), intent(out) :: c(:,:)
      !local variables
      integer, dimension(2) :: dim1, dim2, dim3
      complex(8) :: alpha, beta
      ! get dimensions
      dim1=shape(a)
      dim2=shape(b)
      dim3=shape(c)

      ! set alpha and beta
      alpha=1.0d0
      beta=0.0d0
      call zgemm('N','N',dim1(1),dim2(2),dim1(2),alpha,a,dim1(1),b, &
        & dim2(1),beta,c,dim3(1))

    end subroutine matprod_matmat

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Matrix-vector multiplication of 2D matrix with 1D vector.
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in]  a   
    !> @param[in]  b 
    !> @param[out] c    
    !---------------------------------------------------------------------------  
    subroutine matprod_matvec(a,b,c)
      implicit none
      complex(8), intent(in) :: a(:,:), b(:)
      complex(8), intent(out) :: c(:)
      ! local variables
      integer :: dim1(2), dim3, dim2, M, N, K
      complex(8) :: alpha, beta
      ! zgemm expects rank-2 B(LDB,*) and C(LDC,*); view the vectors as
      ! single-column matrices
      complex(8) :: b_(size(b),1), c_(size(c),1)
      !get dimensions
      dim1=shape(a)
      dim2=size(b)
      dim3=size(c)
      
      !set alpha & beta
      alpha=1.0d0
      beta=0.0d0
      ! set M,N and K
      M=dim1(1)
      N=1
      K=dim1(2)
      b_(:,1)=b
      call zgemm('N','N', M, N, K, alpha, a, dim1(1), b_, dim2, beta, &
       & c_, dim3)
      c=c_(:,1)
    end subroutine matprod_matvec

  
    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> @brief
    !> Scalar product of two 1D vectors.
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in]  a   
    !> @param[in]  b 
    !> @param[out] c    
    !---------------------------------------------------------------------------  
    subroutine matprod_vecvec(a,b,c)
      implicit none
      complex(8), intent(in) :: a(:), b(:)
      complex(8), intent(out) :: c
      !local variables
      complex(8) :: alpha, beta
      ! zgemm expects rank-2 A(LDA,*), B(LDB,*) and C(LDC,*); view the
      ! vectors/scalar as single-column/1x1 matrices
      complex(8) :: a_(size(a),1), b_(size(b),1), c_(1,1)

      alpha=1.0d0
      beta=0.0d0
      a_(:,1)=a
      b_(:,1)=b
      call zgemm('C','N',1,1,size(a),alpha,a_,size(a),b_,size(b),beta,c_,1)
      c=c_(1,1)
    end subroutine matprod_vecvec
end module
