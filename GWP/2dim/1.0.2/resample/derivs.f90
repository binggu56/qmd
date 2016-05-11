
! --- potential and force fields

      subroutine derivs(x,Ntraj,v,dv,ddv)
      
      use cdat, only : ipot 

      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Ntraj
      real*8,intent(IN) :: x(Ntraj)
      real*8,intent(OUT) :: v(Ntraj),dv(Ntraj),ddv(Ntraj)

! --- eckart barrier 
!      v0 = 16d0
!      a = 0.8d0
!
!      do i=1,ntraj
!        v(i) = v0/cosh(x(i)/a)**2
!        dv(i) = -2d0*v0*tanh(x(i)/a)/a/(cosh(x(i)/a))**2
!        ddv(i) = 4d0*v0*tanh(x(i)/a)**2/a**2/cosh(x(i)/a)**2-2d0*v0* &
!                 (1d0-tanh(x(i)/a)**2)/a**2/cosh(x(i)/a)**2
!      enddo
      if (ipot == 1) then 
! ----- quatic oscilator
        ak = 1.d0
        ak2 = 0.01d0 

        do i=1,Ntraj
          dv(i) = ak*x(i) + 4d0*ak2*x(i)**3 
          v(i)  = ak*x(i)**2/2d0 + ak2*x(i)**4
          ddv(i) = ak + 12d0*ak2*x(i)**2 
!          ddv(i) = 0d0 
        enddo 

      elseif(ipot == 2) then 

!------morse potential with extra bound potential at the right side 
      de = 0.176d0
      x0 = 1.4d0
      a = 1.02d0
      beta = 100d0 
      bound = 5d0 
      b = 0d0 

      do i=1,Ntraj
        d = (1d0-exp(-a*(x(i)-x0)))
        v(i) = de*d**2
        dv(i) = 2*de*d*a*exp(-a*(x(i)-x0))+b*exp(-beta*(x(i)-bound)**2)
        ddv(i)=2*de*(-d*exp(-a*((x(i)-x0)))*a**2+(exp(-a*(x(i)-x0)))**2*a**2)       
      enddo
      
      endif 
!---------doulbe well---------------------------
!      eta = 1.3544d0
!      do i=1,Ntraj
!        v(i) = 1d0/16d0/eta*x(i)**4-0.5d0*x(i)**2
!        dv(i) = 1d0/4d0/eta*x(i)**3-x(i)
!        ddv(i) = 3d0/4d0/eta*x(i)**2-1d0
!      enddo

      return
      end subroutine derivs 
