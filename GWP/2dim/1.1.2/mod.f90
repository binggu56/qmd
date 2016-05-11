      module cdat

      integer*4 :: ipot 

      real*8,public,parameter :: PI=4.d0*atan(1.0)
      
      integer,parameter :: nmax = 16 

      complex*16,public,parameter::im=(0d0,1d0)

      real*8 :: al,beta 
 
      integer*4 ::  Ntraj,kmax,kout,idum1,order,nbx,nby
      real*8 :: dt,amx,amy,ax0,ay0,ax,ay,qx0,qy0,px0,py0 
      real*8 :: xmin,xmax,ymin,ymax
      save
      end module cdat

