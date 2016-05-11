
      double complex function corr(nb,c,q,p,s)

      use cdat, only : im,al,pi

      implicit real*8 (a-h,o-z)

      real*8, intent(in) :: q(nb),p(nb),s(nb)
      complex*16, intent(in) :: c(nb)

      complex*16 :: mat(nb,nb),z,z0,c0(nb)
      common/wave/q0,p0,a0,sigma

      a = al
      alfa = 2d0*a0
      as=a+alfa
      av=a*alfa/as
      an=dsqrt(2d0*dsqrt(av/as))
      an0=dsqrt(dsqrt(a/pi))
      an2=dsqrt(dsqrt(alfa/pi))
!---- compute c0_k = <psi0|g_k>
      do j=1,nb
        dq=q(j)-q0
        dp=p(j)-p0
        d2=-0.5d0*av*dq*dq
        d0=-0.5d0/as*dp*dp
        d1=(alfa*p(j)+a*p0)/as*dq
        c0(j)=exp(d2+d0+im*d1)*an
      enddo

      corr = (0d0,0d0)
      do i=1,nb
        corr = corr+conjg(c0(i))*exp(im*s(i))*c(i)
      enddo

      end function corr 
