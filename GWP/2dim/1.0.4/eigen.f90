      subroutine eigen(N,A)
      implicit real*8(a-h,o-z)
      INTEGER*4,intent(in) :: N
      complex*16, intent(in) :: A(n,n)
      complex*16 :: mat(n,n) 

      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
      INTEGER          INFO, LWORK
!     .. Local Arrays ..
!     RWORK dimension should be at least MAX(1,3*N-2)
      DOUBLE PRECISION W( N ), RWORK( 3*N-2 )
      COMPLEX*16       WORK( LWMAX )
      LDA = N
!     .. External Subroutines ..
!      EXTERNAL         ZHEEV
!      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
!
!     .. Intrinsic Functions ..
!      INTRINSIC        INT, MIN
!
!     .. Executable Statements ..
!      WRITE(*,*)'ZHEEV Example Program Results'
!
!     Query the optimal workspace.

      mat = A 

      LWORK = -1
      CALL ZHEEV( 'Vectors', 'Lower', N, mat, LDA, W, WORK, LWORK, RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

!     Solve eigenproblem.

      CALL ZHEEV( 'Vectors', 'Lower', N, mat, LDA, W, WORK, LWORK, RWORK, INFO)

!     Check for convergence.

      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
     
      ak = maxval(w)/minval(w)
      if (ak > 1d7) then
        stop
      else
        print *,'condition number', maxval(w)/minval(w)  
      endif 

      return 
      end subroutine

!     End of ZHEEV Example.
