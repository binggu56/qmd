
! --- initialize parameters 

      subroutine init_pars()
      
      use cdat 

      implicit real*8 (a-h,o-z)  

      open(5,file='IN')
      read(5,*) nbx,nby 
      read(5,*) kmax,dt,kout
!      read(5,*) xmin,xmax 
      read(5,*) amx,amy 
      read(5,*) ax0,ay0,ax,ay,beta 
      read(5,*) qx0,qy0 
      read(5,*) px0,py0 
      read(5,*) idum1
      read(5,*) order 
      close(5)

!      ipot = 1 ! 1. HO 2. Morse  

      end subroutine
