      program main 
      implicit real*8(a-h,o-z)
      real*8 :: ki
     	real*8,allocatable :: ke(:),V(:),x(:),y(:),px(:),py(:),w(:),c(:)
      real*8,allocatable :: qfx(:),qfy(:),qp(:),fx(:),rx(:),      & 
                           fy(:),s(:),dpx(:),dpy(:),s0(:),dr(:),ddr(:), & 
                           ap(:),ar(:),ddp(:)
     	real*8 :: f(3),err(2) 

      real :: gasdev
      complex*16 :: cor,im
      common/smaple/xmin,xmax,dx

      open(100,file='energy.out')
      open(102,file='xav.out')
	open(103,file='pes.out')
	open(101,file='p.out') 
	open(104,file='ke.out')
      open(105,file='weight.out')
      open(106,file='traj.out')
      open(107,file='cor.out')
      open(119,file='err.out')
      open(120,file='px.out')

      open(5,file='IN')
      read(5,*) Ntraj
      read(5,*) kmax,dt
      read(5,*) am
      read(5,*) idum1
      read(5,*) ax
      read(5,*) qx0,p0 
      read(5,*) fc
      close(5)

      allocate(ke(Ntraj),V(Ntraj),x(Ntraj))
      allocate(y(Ntraj),px(Ntraj),py(Ntraj),s(Ntraj),w(Ntraj),qfx(Ntraj)) 
      allocate(qfy(Ntraj),qp(Ntraj),fx(Ntraj),fy(Ntraj),dpx(Ntraj),dpy(Ntraj),s0(Ntraj))
      allocate(c(ntraj),rx(ntraj),dr(ntraj),ddr(ntraj),ap(ntraj),ar(ntraj),ddp(ntraj))

      dt2 = dt/2d0
      t = 0d0
      s = 0d0
      pi = 4.d0*atan(1.d0)
      pow = 6d0
      
      im = (0d0,1d0)
!      xmin = qx0-dsqrt(pow/ax)
!      xmax = qx0+dsqrt(pow/ax)

      ww = 0d0

      var = sqrt(1.0/2.0/ax) 

      do i=1,Ntraj
1100    x(i) = gasdev(idum1)
!        if(abs(x(i)) > 4.0) goto 1100 

        x(i) = x(i)/sqrt(4d0*ax)+qx0
        
!        y(i)=y(i)/sqrt(4d0*ay)+qy0
!        x(i) = xmin+dx*(i-1)
!        c(i) = -ax*(x(i)-qx0)**2+log(dsqrt(2d0*ax/pi))/2d0
!        w(i) = exp(2d0*c(i))*dx
        w(i) = 1d0/Ntraj
        rx(i) = -2d0*ax*(x(i)-qx0) 
        ww = ww+w(i)
      enddo 
      
      write(*,*) 'Weights = ', ww
       
      do i=1,Ntraj
        px(i) = p0 
!        py(i) = 2d0*ay*(y(i)-qy0)
!        s0(i) = ax*(x(i)-qx0)**2+ay*(y(i)-qy0)**2
      enddo


!        print out the initial conditions        
        write(*,*) 'Initial Conditions'
        print *,'ax =',ax
        write(*,*) 'Number of trajectories =', Ntraj
        write(*,*) 'Time steps =',kmax,dt
        write(*,*) 'Mass =',am 
       
      s = 0d0
! 	initial value for action function s
	do i=1,Ntraj
		s(i) = px(i)*(x(i)-x0)
	enddo


! --- print intial trajectories and potential
	xmin = 0.1 
	dx = 0.01 
    do i=1,300 
      x_grid = xmin + (i-1) * dx
      call derivs(x_grid, v0,dv) 
      write(103,10000) x_grid, v0, dv 
    enddo

      write(*,*) 'friction coeff =',fc
      
!     begin the time step loop
        do 10 k=1,kmax

          t=t+dt   ! increase t by dt

          call fitp(am,x,px,rx,w,Ntraj,dpx,ddp,dr,ddr,qp,ap,ar,err)

          do i=1,Ntraj
            x(i) = x(i) + px(i) / am * dt

            call derivs(x(i),v0,dv) 

            v(i) = v0 ! potential energy of i-th trajectory 

            px(i) = px(i) + ( - dv - fc * px(i) + 1d0/2d0/am * (2d0 * ar(i)*dr(i) + ddr(i)))*dt
!            c(i) = c(i)+(-1d0/m1*(rx(i)*px(i)+dpx(i)/2d0))*dt
            rx(i) = rx(i) - 1d0/2d0/am * (2.0 * dpx(i) * ar(i) + ddp(i))*dt
          enddo
  
          do i=1,Ntraj
            s(i) = s(i)+(- px(i)**2/2d0/am - v(i)-qp(i))*dt
          enddo
         

!
!	do  i=1,Ntraj
!          px(i)=px(i)+(-fx(i)-qfx(i))*dt2-px(i)*dpx(i)/m1*dt2
!          py(i)=py(i)+(-fy(i)-qfy(i))*dt2-py(i)*dpy(i)/m2*dt2
!        enddo
! update weights for guassian trajectories
!	do i=1,Ntraj
!	  w(i) = exp(-2d0*(s(i)-s0(i)))*1d0/Ntraj
!	enddo
        
         write(119,10000) t,(err(i),i=1,2)   

! update potential, kinetic, and total energy for each trajectory	
         Ek = 0d0 
         do i=1,Ntraj
           Ek =  Ek + px(i)**2/(2d0*am) * w(i) 
         enddo
! print out a random trajectory	
        
         write(106,10000) t,(x(i),i=1,20)

! calculate the total energy, the sum of all the trajectories	
          call aver(Ntraj,w,V,VAve)
          call aver(Ntraj,w,qp,UAVe)
        if(UAve < 0) then 
			print *,'Quantum potential = ',UAve, ' Cannot be a negative number'
			stop 
		endif 

	  ww = 0d0
	  do i=1,Ntraj
	    ww = ww + w(i)
	  enddo
	  write(105,10000) t,ww

      
! --- autocorrelation 
      cor = (0d0,0d0)
      do i=1,ntraj
        cor = cor+w(i)*exp(2d0*s(i)*im)
      enddo
        write(107,10000) t,cor,abs(cor)

        Etot = (Ek + VAve + UAve)
        write(100,10000) t, Ek, VAve, UAve, Etot

10	end do

      x1 = 0d0 
      x2 = 0d0 

      do i=1,ntraj 
        x1 = x1 + x(i)*w(i) 
        x2 = x2 + x(i)**2*w(i) 
      enddo 

      write(*,*) 'Expectation of x, var',x1,x2 - x1**2 
      write(*,*) 'Total Energy = ', Etot, 'Hartree.' 

      do i=1,ntraj
        write(120,10000) x(i),px(i)
      enddo

10000 format(100(e14.7,1x))
1001  format(A10,20(f10.6,1x))

      close(101)
      close(102)
      close(103)
      close(104)
      close(105)
      close(106)
      
      return
      end program main

        
!----------------------------------------------
! average over trajectories y = sum(w(i)**x(i))
!----------------------------------------------
        subroutine aver(Ntraj,w,x,y)
        implicit real*8(a-h,o-z)
        real*8 :: x(Ntraj),w(Ntraj)
        y = 0d0

        do i=1,Ntraj 
          y = y+x(i)*w(i) 
        enddo
      
        return
        end subroutine

!     -------------------------------------
!     potential and derivatives
!     -------------------------------------
      subroutine derivs(x,v0,dv)
      implicit real*8(a-h,o-z)
      real*8,   intent(in) :: x
      character*20 :: PES 
     
      PES = 'Morse' 
      
 
      a = 1.02
      x0 = 1.4 
      De = 0.176
    
        d = (1.0-exp(-a*(x-x0)))
    
        v0 = De * (d**2) 
        dv = 2.0 * De * d * a * exp(-a*(x-x0))
    
!        ddv = 2.0 * De * (-d*exp(-a*((x-x0)))*a**2 + (exp(-a*(x-x0)))**2*a**2)
          

      return 
      end subroutine
      
      
