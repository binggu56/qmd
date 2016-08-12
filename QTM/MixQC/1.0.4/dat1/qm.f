
!     QSATS version 1.0 (3 March 2011)

!     file name: eloc.f

! ----------------------------------------------------------------------

!     this computes the total energy and the expectation value of the
!     potential energy from the snapshots recorded by QSATS.

! ----------------------------------------------------------------------

      module domain 

      integer (kind = 4), parameter :: id = 2 ! classical DoF defining the domain function 
      integer (kind = 4), parameter :: nd = 4 ! number of domains  
      
!      integer (kind = 4), parameter :: 

      real (kind = 8), dimension(nd+1) ::  ydom
      real (kind = 8) ymin,ymax,dy 
      
      save 
      end module 
      
      module cdat 


      integer (kind = 4), parameter :: nq = 1 ! number of QDoF  
      real (kind = 8),parameter ::  half = 0.5d0, one = 1.0d0 
      real (kind = 8),parameter ::  pi = 4.0*atan(1d0)  

      
      contains 
       
 
      double precision function step(x,y,z) 
      implicit double precision (a-h,o-z) 
      if(z > x .and. z < y) then 
         step = 1d0 
      else
        step = 0d0 
      endif 
      return 
      end function     
 
      end module 
! ----------------------------
! --- main program 
! ------------------------------
      program eloc
      use cdat
      use domain  

      implicit double precision (a-h, o-z)

!------------trajectory initial allocate------------------------
      real*8, dimension(:,:),allocatable :: x,p,du,ap,ar,rp,fr

      real*8, dimension(:,:),allocatable :: dv 
      real*8, dimension(:),allocatable :: am,p0,w,alpha,x0,pe,ke,u,cf
      real*8 :: gasdev
      
      integer (kind = 4), dimension(:), allocatable :: idum 

      real (kind = 8), dimension(:), allocatable :: enk,env 
      real (kind = 8), dimension(:), allocatable :: q  

!------------------------------------------------------
!      character*4 :: atom(NATOMS)
!      real*8 :: dvl(NATOM3)
!      integer*4  :: ntraj,ndim,kmax

c --- this common block is used to enable interpolation in the potential
c     energy lookup table in the subroutine local below.


      open(100, file = 'en.dat'  )
      open(101, file = 'xoutput' )
      open(102, file = 'traj.dat')
      open(103, file = 'den.dat' ) ! density overlap
      open(104, file = 'xy.dat'  ) 
      open(105, file = 'xy0.dat'  ) 
      
      open(110, file = 'ratio.dat')
      open(111, file = 'ndom.dat' )

!------read parameters from IN----------------------      
      open(10,file='IN')

      read(10,*) Ntraj
      read(10,*) ndim 
      read(10,*) kmax,dt
!      read(10,*) a0

      close(10)

!      Ndim = NATOM3

c --- initialization.

!      call tstamp

      allocate(dv(ndim,ntraj))
!-----------------------------------------------------
      allocate(env(ntraj),enk(Ntraj),x(Ndim,Ntraj),
     +         p(Ndim,Ntraj),du(Ndim,Ntraj),w(Ntraj))

      allocate(p0(Ndim),alpha(Ndim),am(Ndim),x0(Ndim),cf(Ndim))        

      allocate( idum(ndim) ) 
      allocate( q(ndim) ) 
     
!      binvrs=one/bin


c ------ compute the local energy and the potential energy.
      
      dt2 = dt/2d0
      t   = 0d0
      pow = 6d0
      am0 = 1d0 

      am(1) = 1d0 
      am(2) = 10d0 

      alpha(1) = 1.118 
      alpha(2) = 11.18 
      
      x0(1) = 0d0 
      x0(2) = -3d0 
      
      p0(1) = 0d0 
      p0(2) = 18d0 


! --- domain setup 
      
      ymin = -4.0d0 
      ymax = 6.0d0 
      dy = (ymax-ymin)/nd 
      do i=2,nd 
        ydom(i) = ymin+(i-1)*dy 
      enddo 
      ydom(1) = -100d0 
      ydom(nd+1) = 100d0

      write(6,1016) ydom 
1016  format('DOMAIN BOUNDARY',20(f9.3,1x)) 

      call seed(idum,Ndim)

      write(*,1010) Ntraj,Ndim,kmax,dt,am(1),am(2)
1010  format('Initial Conditions'//,
     +       'Ntraj   = ' , i6/, 
     +       'Ndim    = ' , i6/,
     +       'kmax    = ' , i6/,
     +       'dt      = ' , f7.4/,
     +       'Mass    = ' , f7.4,f7.4/)

      write(*,1011) alpha(1),alpha(2),x0(1),x0(2),p0(1),p0(2)
1011  format('Initial Wavepacket'//,
     +       'GWP width   = ' , 2(f7.4,1x)/,
     +       'GWP center  = ' , 2(f7.4,1x)/,
     +       'Momentum    = ' , 2(f7.4,1x)/)


! --- initial Lagrangian grid points (evolve with time)
      do i=1,Ntraj
        do j=1,Ndim
1100      x(j,i)=gasdev(idum(j))
          x(j,i)=x(j,i)/dsqrt(4d0*alpha(j))+x0(j)
!          if((x(j,i)-x0(j))**2 .gt. pow/2d0/alpha(j)) goto 1100
        enddo
      enddo
      
      do i=1,ntraj 
        write(105,1000) x(1,i),x(2,i) 
      enddo 

! --- initial momentum for QTs and weights

      do i=1,Ntraj
        p(1,i) = p0(1) 
      enddo 
      do i=1,ntraj 
            p(2,i) = p0(2)+2d0*alpha(2)*(x(2,i)-x0(2)) 
!            rp(j,i) = -2d0*alpha(j)*(x(j,i)-x0(j))
      enddo

      do i=1,ntraj 
        w(i) = 1d0/dble(ntraj)
      enddo

!---- expectation value of x(1,ntraj)--------
      av = 0d0
      do i=1,ntraj
        av = av+w(i)*x(1,i)
      enddo
      write(*,6008) av
6008  format('Expectation value of x(1,ntraj) = ', f10.6)





! --- half timestep for p 

      call fit(am,x,p,w,ndim,ntraj,enq,du)
      call local(ndim,ntraj,x,env,dv) 
      do i=1,ntraj         
          do j=1,ndim 
            p(j,i) = p(j,i)-(dv(j,i)+du(j,i))*dt2 
          enddo 
      enddo 

! --- trajectories propagation	
      do 10 kt=1,kmax

        t=t+dt



! ----- half-step increments of momenta & full step increment of positions
        do i=1,Ntraj

          
!          call derivs(q,vloc,dv)
!          call long_force(q,vlong,dvl)
!          write(*,*) 'test2'
!          env(i) = vloc 

          do j=1,Ndim
!              p(j,i)=p(j,i)+(-dv(j)-du(j,i))*dt2
              x(j,i)=x(j,i)+p(j,i)*dt/am(j)
          enddo   
        enddo 
    
!          write(*,*) 'test3',x(1,2),x(2,3)

! ----- half-step increments of momenta


          call local(ndim,ntraj,x,env,dv)
!          call long_force(q,vlong,dvl)
          call fit(am,x,p,w,ndim,ntraj,enq,du)
       
        do i=1,Ntraj
          do j=1,Ndim
            p(j,i)=p(j,i)+(-dv(j,i)-du(j,i))*dt
          enddo
        enddo

!-------update potential, kinetic, and total energy
          do i=1, Ntraj
            enk(i) = 0d0
            do j=1,Ndim
              enk(i) = p(j,i)**2/(2d0*am(j)) + enk(i) 
            enddo
          enddo

          write(101,1000) t,(x(1,i),i=1,20)
          write(102,1000) t,(x(j,4),j=1,20)

        call aver(Ntraj,w,enk,enk_ave)
        call aver(Ntraj,w,env,env_ave)
!        call aver(Ntraj,w,u,qu)
      
        ent = env_ave+enk_ave+enq 

!        do i=1,ntraj
!          write(105,1000) dble(i),ev(i)
!        enddo

        write(100,1000) t,enk_ave,env_ave,enq,ent 
        call flush(100)

! --- density overlap 
      d0 = 0d0 
      d1 = 0d0 
      ak1 = 5d0 
      ak2 = 15d0 

      alfa = dsqrt(am(1)*ak1)/2d0 
      beta = dsqrt(am(1)*ak2)/2d0 

      do i=1,ntraj
        if(x(2,i) > 0d0) then  
        d0 = d0 + w(i)*exp(-alfa*x(1,i)**2)
        d1 = d1 + w(i)*exp(-beta*x(1,i)**2)
       endif 
      enddo 
      d0 = d0*dsqrt(alfa/pi) 
      d1 = d1*dsqrt(beta/pi) 
      

      write(103,1000) t,d0,d1 

! --- store temporary data
!      open(106,file='temp.dat')
!      if (mod(kt,1000) .eq. 0) then
!        do i=1,Ntraj
!          do j=1,Ndim
!            write(106,1000) x(j,i),p(j,i)
!            call flush(106)
!          enddo
!        enddo
!      endif

!      close(106)

10    enddo


      do i=1,ntraj 
        write(104,1000) x(1,i),x(2,i) 
      enddo 

      write(*,1020) ent 
1020  format('Total Energy =', f10.5/ ,
     +       'MISSION COMPLETE.')

      
!          call local(q,vloc,dv)

!       convert to eV 
!          vloc = vloc*27.211d0
c ------ convert to kelvin per atom.

!         tloc=tloc/(3.1668513d-6*dble(NATOMS))
!          vloc=vloc/(3.1668513d-6*dble(NATOMS))

c ------ accumulate the results.

!         vtavg(k)=vtavg(k)+vloc
!         vtavg2(k)=vtavg2(k)+(vloc)**2
!         etavg(k)=etavg(k)+tloc+vloc
!         etavg2(k)=etavg2(k)+(tloc+vloc)**2

!350      continue

!      end do

!      goto 300

c --- account for overshooting.

!600   loop=loop-1

!      write (6, 6600) loop
!6600  format ('number of snapshots = ', i6/)

c --- compute the averages and standard deviations.

!      do k=1, NREPS, 11
!
!      vtavg(k)=vtavg(k)/dble(loop)
!      vtavg2(k)=vtavg2(k)/dble(loop)
!      etavg(k)=etavg(k)/dble(loop)
!      etavg2(k)=etavg2(k)/dble(loop)
!
!      vsd=sqrt(vtavg2(k)-vtavg(k)**2)
!      esd=sqrt(etavg2(k)-etavg(k)**2)
!
!      write (6, 6610) k, 'VAVG = ', vtavg(k)
!6610  format ('replica ', i3, 1x, a7, f10.5, ' Kelvin')
!          write (100,1000) q(1),vloc
!       end do
!1001    format('potential at equilibrim configuration [K/atom] =',
!     &         f12.7)
!
!      write (6, 6610) k, 'V SD = ', vsd
!
!      write (6, 6610) k, 'EAVG = ', etavg(k)
!
!      write (6, 6610) k, 'E SD = ', esd
!
!      end do
1000  format(20(e14.7,1x))
      end program


c ----------------------------------------------------------------------

c     quit is a subroutine used to terminate execution if there is
c     an error.

c     it is needed here because the subroutine that reads the parameters
c     (subroutine input) may call it.

c ----------------------------------------------------------------------

      subroutine quit

      write (6, *) 'termination via subroutine quit'

      stop
      end subroutine
! ---------------------------------------------------------------
!     random number seeds
! ---------------------------------------------------------------
      subroutine seed(idum,Ndim)
      implicit real*8(a-h,o-z)
      integer*4, intent(IN) :: Ndim
      integer*4, intent(OUT) :: idum(Ndim)

      do i=1,Ndim
        idum(i) = 5 + i
      enddo

      return
      end subroutine
! ---------------------------------------------------------------
!     average over trajectories
! ---------------------------------------------------------------
      subroutine aver(Ntraj,w,x,y)
      implicit real*8(a-h,o-z)
      real*8 :: x(Ntraj), w(Ntraj)

      y = 0d0

      do i=1,Ntraj
        y = y+x(i)*w(i)
      enddo
      
      return
      end subroutine
! ----------------------------------------------------------------
