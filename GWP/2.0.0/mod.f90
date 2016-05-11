      module cdat

      integer*4 :: ipot 

      real*8,public,parameter :: PI=4.d0*atan(1.0)
      
      integer,parameter :: nmax = 16 

      complex*16,public,parameter::im=(0d0,1d0)

      real*8 :: al,beta 
 
      integer*4 ::  Ntraj,kmax,kout,idum1,order
      real*8 :: dt,am,a0,q0,p0
      real*8 :: xmin,xmax
      save
      end module cdat

