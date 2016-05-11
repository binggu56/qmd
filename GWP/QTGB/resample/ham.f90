      subroutine pot(x,v,ddv)
      implicit real*8(a-h,o-z)

      ipot = 2

      if (ipot == 0) then
        ep = 0.0
        v = x**2/2d0 + ep*x**4/4.0
        ddv = 1d0 + 3.0*x**2*ep

      elseif(ipot == 1) then

        eta = 1.3544d0
        v = 1d0/16d0/eta*x**4-0.5d0*x**2
        ddv = 3d0/4d0/eta*x**2-1d0
      elseif(ipot == 2) then

        D = 16d0
!        a = 1.3624d0
        a = 1d0 

        v = D/cosh(a*x)**2
        ddv = 6d0*D*a**2*sinh(a*x)**2/cosh(a*x)**4-2d0*D*a**2/cosh(a*x)**2
      endif

      return
      end subroutine

      subroutine ham(c,nb,q,p,s,dc)

      use cdat 
      implicit real (kind = 8) (a-h,o-z) 

      integer*4, intent(in) :: nb
      complex*16, dimension(nb), intent(in) :: c
      real*8, intent(in), dimension(nb) :: q,p,s

      complex*16,intent(out) :: dc(nb)

      real*8 :: v,dv,ddv

      complex*16, dimension(nb,nb) :: h,mat
      real*8 :: aux(nb,nb)

      complex*16 :: z,z0,d0,d1,d2,b(nb)
      real*8 b1(nb,2)

      call multi(nb,q,p,s,mat)

      anrm = dsqrt(dsqrt(al/pi))
!      eta = 1.3544d0

      do j=1,nb
        do k=1,nb
          dq = q(k)-q(j)
          z = dq/2d0
          z0 = exp(-al*dq**2/4d0)
          d0 = al/(4d0*am)*(1d0-al/2d0*dq**2)
          y = (q(j)+q(k))/2d0

! --- local harmonic approximation 
          call pot(y,v,ddv)
          d1 = v + ddv/(4d0*al)

! --- analytical potential integration for double-well 
!          d1 = 1d0/(16d0*eta)*(y**4 + 3d0*y**2/al + 3d0/(4d0*al**2)) -
!          0.5d0*(y**2+1d0/(2d0*al))

          d2 = im*p(k)*al/(2d0*am)*dq
          h(j,k) = z0*(d0+d1+d2)
        enddo
      enddo

      b = matmul(h,c)

      do i=1,nb
        b1(i,1) = real(b(i))
        b1(i,2) = imag(b(i))
      enddo

!      call inverse(mat,mativs,nb)

!-------save the overlap matrix-----------
      do i=1,nb
      do j=1,nb
        aux(i,j) = real(mat(i,j))
      enddo
      enddo

      call dposv('U', nb, 2, aux, nb, b1, nb, INFO)
      if(INFO .ne. 0) then
        write(*,*) 'dC/dt: Matrix fails, INFO =',INFO
        stop
      endif
!------compute increment dc------------
      do j=1,nb
        dc(j) = -im*(b1(j,1)+im*b1(j,2))
      enddo

      return
      end subroutine

