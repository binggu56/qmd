      subroutine quantum_potential(m1,x,w,Ntraj,qp,qfx,c)
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Ntraj

      integer INFO
      integer, parameter :: nb=3
      real*8,intent(IN) :: m1
      real*8,intent(IN) :: w(Ntraj),x(Ntraj)
      real*8, intent(OUT) :: qfx(Ntraj),qp(Ntraj),c(nb,1)
      real*8 :: f(nb),rx(Ntraj)
      real*8 ::  s(nb,nb),dr(ntraj)

! define c matrix = -1/2*df/dx
      c(1,1) =-0.5d0
      c(2,1) = 0d0
c      c(3,1) = 0d0

      do i=1,ntraj
        c(2,1) = c(2,1)+w(i)*2d0*x(i)
c        c(3,1) = c(3,1)+w(i)*3d0*x(i)**2
      enddo

      c(2,1) = -0.5d0*c(2,1)
c      c(3,1) = -0.5d0*c(3,1)
c      c(4,1) = 0d0

!     f*transpose(f) matrix
      s=0d0
      do i=1,Ntraj
      f=(/x(i),x(i)**2,1d0/)
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
        rx(i) = c(1,1)*x(i)+c(2,1)*x(i)**2+c(3,1)

! calculate quantum potential
        dr(i) = c(1,1)+2d0*c(2,1)*x(i)
        qp(i)=ams*(rx(i)**2+dr(i))
        qfx(i)=-ams*(2d0*rx(i)*dr(i)+2d0*c(2,1))
      end do

      return
      end subroutine
