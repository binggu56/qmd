      subroutine qpot(nb,ddv,c,am,p,q,w,u,du)
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: nb
      real*8,intent(IN) :: w(nb),q(nb),am,p(nb),ddv(nb)
      complex*16,intent(in) :: c(nb)
      real*8, intent(OUT) :: u(nb),du(nb)
      complex*16 s,z0,z,psi(np)
      common/wave/al,q0,p0,a0
      common/grid/xmin,xmax,dx,np

      u = 0d0

!      do i=1,nb
!        s = (0d0,0d0)
!        do j=1,nb
!          do k=1,nb
!            dq = q(j) - q(k)
!            z = dq/2d0-im*(p(j)+p(k))/2d0/al
!            z0 = exp(-al*z*conjg(z)+im/2d0*(p(j)+p(k))*dq)          
!            s = s+conjg(c(j))*c(k)*(z0*z+(q(k)-q(i))*z0)
!          enddo
!        enddo
!        print *,s
!        du(i) = real(s)*ddv(i)
!      enddo


!--------numerical integration---------------
      do i=1,np
        xi = xmin+dx*(i-1)
        psi(i) = (0d0,0d0)
        do j=1,nb
          psi(i) = psi(i)+c(j)*exp(-al/2d0*(xi-q(j))**2+
     +             im*p(j)*(xi-q(j)))*
     +            dsqrt(dsqrt(al/pi))
        enddo
      enddo

      do k=1,nb
        d = 0d0
        do i=1,np
          xi = xmin+dx*(i-1)
          d = d+abs(psi(i))**2*(xi-q(k))*ddv(k)*dx
        enddo
        du(k) = d
      enddo

      return
      end subroutine
