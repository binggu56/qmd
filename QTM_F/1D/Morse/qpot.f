      subroutine qpot(m1,x,w,Ntraj,qp,qfx)
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Ntraj

      integer INFO
      integer, parameter :: nb=4
      real*8,intent(IN) :: m1
      real*8,intent(IN) :: w(Ntraj),x(Ntraj)
      real*8, intent(OUT) :: qfx(Ntraj),qp(Ntraj)
      real*8 :: rx(ntraj)
      real*8 :: f(nb),s(nb,nb),c(nb,1)

! define c matrix = -1/2*df/dx
      c(1,1) =-0.5d0
      c(2,1) = 0d0
      c(3,1) = 0d0

      do i=1,ntraj
        c(2,1) = c(2,1)+w(i)*2d0*x(i)
        c(3,1) = c(3,1)+w(i)*3d0*x(i)**2
      enddo

      c(2,1) = -0.5d0*c(2,1)
      c(3,1) = -0.5d0*c(3,1)
      c(4,1) = 0d0

!     f*transpose(f) matrix
      s=0d0
      do i=1,Ntraj
      f=(/x(i),x(i)**2,x(i)**3,1d0/)
        do m=1,nb
          do n=1,nb
            s(m,n)=w(i)*f(m)*f(n)+s(m,n)
          end do
        end do
      end do

! calculate matrix c(t)
        call DPOSV('U',nb,1,s,nb,c,nb,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, "matrix fails"
        stop
        end if

! the momentum operator r=cf
      ams = -1d0/2d0/m1
      do i=1,Ntraj
        rx(i) = c(1,1)*x(i)+c(2,1)*x(i)**2+c(3,1)*x(i)**3+c(4,1)

! calculate quantum potential
        dr = c(1,1)+2d0*c(2,1)*x(i)+c(3,1)*3d0*x(i)**2
        qp(i)=ams*(rx(i)**2+c(1,1)+c(2,1)*2d0*x(i)+c(3,1)*3d0*x(i)**2)
        qfx(i)=-ams*(2d0*rx(i)*dr+2d0*c(2,1)+c(3,1)*6d0*x(i))

      end do

      return
      end subroutine
