! ----------------------------------------------------------------------
!
!     this subroutine computes the local energy and potential energy
!     of a configuration.
!
! ----------------------------------------------------------------------

      subroutine local(ndim,ntraj,x, env, dv)

      implicit double precision (a-h, o-z)

      real (kind = 8) x(ndim,ntraj) 

      real*8, intent(out) :: env(ntraj),dv(ndim,ntraj)
      
      ak1 = 5.0 
      ak2 = 15.0 
      b = 0.5d0 
      v0 = 16d0 
      g = 1.3624d0 
      
      do i=1,ntraj 

      ak = 0.5*(ak1+ak2)+0.5*(ak2-ak1)*tanh(b*x(2,i))

      env(i) = 0.5*ak*x(1,i)**2 + v0/cosh(g*x(2,i))**2  
      
      dv(1,i) = ak*x(1,i) 
      dv(2,i) = x(1,i)**2*(ak2-ak1)*(1d0-tanh(b*x(2,i))**2)*b/4d0- & 
                2d0*v0*g*sinh(g*x(2,i))/cosh(g*x(2,i))**3
      enddo 
      
              
      return
      end subroutine

