      subroutine qpot(Nb,ddv,c,am,p,x,w,u,du) 
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Nb
      real*8,intent(IN) :: w(Nb),x(Nb),am,p(nb),ddv(nb)
      complex*16 c(nb)
      real*8, intent(OUT) :: du(nb),u(Nb)

      an = 0d0
      x1 = 0d0
      x2 = 0d0
      do i=1,nb
        an=an+w(i)
        x1=x1+w(i)*x(i)
        x2=x2+w(i)*x(i)**2
      enddo

      b1=-0.5d0*an*an/(an*x2-x1*x1)
      b2 = 0.5d0*x1*an/(an*x2-x1*x1)
      do i=1,Nb
        r=b1*x(i)+b2
        u(i) = -(r*r+b1)/2d0/am
        du(i) = -b1*r/am
      enddo
      


        !common/pot_par/m1,m2
! define c matrix 
!        c(1,1)=-0.5d0
!C        c(1,2)=0d0
!        c(2,1)=0d0
!C       c(2,2)=-0.5d0
!
!! quantum force
!        s=0d0
!        do i=1,Ntraj
!        f=(/x(i),1d0/)
!        do m=1,2
!        do n=1,2
!        s(m,n)=w(i)*f(m)*f(n)+s(m,n)
!        end do
!        end do
!        end do
!
!        !do m=1,3
!        !print *, s(m,1),s(m,2),s(m,3)
!        !enddo
!        !STOP
!      print *,'qpot'
!      print *,s(1,1),s(1,2)
!      print *,c(1,1),c(2,1)
!! calculate matrix c(t)
!        call DPOSV('U',2,1,s,2,c,2,INFO)
!        if(INFO/=0) then
!        print *, "info=",info
!        print *, "QPOT: matrix fails"
!        stop
!        end if
!
!! the momentum operator r=cf
!        do i=1,Ntraj
!        rx(i)=c(1,1)*x(i)+c(2,1)
!
!! calculate quantum potential
!        qp(i)=-rx(i)**2/(2d0*m1)
!     &     -c(1,1)/(2d0*m1)
!        qfx(i)=rx(i)*c(1,1)/m1
!
!        end do

        return
        end subroutine
