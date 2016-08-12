	module pH2 
	real*8, parameter :: a2bohr = 0.52917721092, hartree_wavenumber = 219474.63 
	real*8, parameter :: re = 3.47005, De = 24.2288 
	real*8, parameter :: r_ref = 4.60, Vmin = -24.2288 
	
	real*8, parameter :: C6 = 5.820364d04, C8 = 2.87052154d05, C10 = 1.80757343d06 
	real*8, parameter, dimension(0:10) :: b = (/-6.631d-02, 1.346d-01, -3.300d-02, 6d0, -1.4d01,  & 
				-1.193d02, 2.290d02, 1.110d03, -1.850d03, -3.5d03, 6.0d03/)
	
	contains 
	
	double precision function damp(r,n)
		implicit real*8(a-h,o-z)
		den = 1.10 
		 
		damp = (1.0 - exp(-3.30 * den * r/n - 0.423 * (den * r)**2/sqrt(dble(n))))**(n-1) 
		 
	end function 
	 
	end module 
	 
	


	program main 
	
	use pH2, only : hartree_wavenumber 

	implicit real*8(a-h,o-z)

	real*8,allocatable :: ke(:),V(:),x(:),y(:),px(:),py(:),w(:),c(:)
	real*8,allocatable :: qfx(:),qfy(:),qp(:),fx(:),rx(:),      & 
					   fy(:),s(:),dpx(:),dpy(:),s0(:),dr(:),ddr(:), & 
					   ap(:),ar(:),ddp(:)
	real*8 :: f(3),err(2) 

	real :: gasdev
	complex*16 :: cor,im
	common/smaple/xmin,xmax,dx

	open(100,file='energy.dat')
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

    ww = 0d0
      do i=1,Ntraj
        x(i) = gasdev(idum1)
        x(i)=x(i)/sqrt(4d0*ax)+qx0
!        y(i)=y(i)/sqrt(4d0*ay)+qy0
!        x(i) = xmin+dx*(i-1)
!        c(i) = -ax*(x(i)-qx0)**2+log(dsqrt(2d0*ax/pi))/2d0
!        w(i) = exp(2d0*c(i))*dx
        w(i) = 1d0/Ntraj
        rx(i) = - 2d0*ax*(x(i)-qx0) 
        ww = ww+w(i)
      enddo 
      
      write(*,*) 'Weights = ', ww
       
       px = p0 


! ---- print out the initial conditions        
	write(*,*) 'Initial Conditions'
	print *,'ax =',ax
	write(*,*) 'Number of trajectories =', Ntraj
	write(*,*) 'Time steps =',kmax,dt
	write(*,100) am, p0 
100	format('Mass = ',f14.7, 'Initial mementum',f14.7) 
       
	s = 0d0
	! --- initial value for action function s
	do i=1,Ntraj
		s(i) = px(i)*(x(i)-x0)
	enddo


! --- print intial trajectories and potential
	xmin = 4 
	dx = 0.02 
    do i=1,300 
		x_grid = xmin + (i-1) * dx
		write(103,10000) x_grid,vpot(x_grid) 
    enddo

      write(*,*) 'friction coeff =',fc

    print *,'initial trajectory', minval(x),maxval(x) 
 
          call derivs(x,ntraj,fx,v)
          call fitp(am,x,px,rx,w,Ntraj,dpx,ddp,dr,ddr,qp,ap,ar,err)    

    print *,'fx',fx(1:10)
    print *,'dr',dr(1:10)
    print *,'dp',dpx(1:10) 

!     begin the time step loop
        do 10 k=1,kmax

          t = t + dt   ! increase t by dt

          do i=1,Ntraj

            px(i) = px(i)+ ( fx(i) + 1d0/2d0/am * (2d0 * rx(i)*dr(i) + ddr(i)))*dt2 - fc * px(i) * dt2 
            rx(i) = rx(i) - 1d0/2d0/am * (2.0 * dpx(i) * rx(i) + ddp(i))*dt2

          enddo
 
        print *,'px',px(1:10) 

          do i=1,ntraj  
            x(i) = x(i) + px(i) / am * dt
          enddo 

        print *,'x',x(1:10) 

          call derivs(x,ntraj,fx,v)
          call fitp(am,x,px,rx,w,Ntraj,dpx,ddp,dr,ddr,qp,ap,ar,err)    
 
          do i=1,Ntraj
          
            px(i) = px(i)+ ( fx(i) - fc * px(i) + 1d0/2d0/am * (2d0 * rx(i)*dr(i) + ddr(i)))*dt2
            rx(i) = rx(i) - 1d0/2d0/am * (2.0 * dpx(i) * rx(i) + ddp(i))*dt2

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

	print *, 'p',px(1:10) 
      
! --- autocorrelation 
      cor = (0d0,0d0)
      do i=1,ntraj
        cor = cor+w(i)*exp(2d0*s(i)*im)
      enddo
        write(107,10000) t,cor,abs(cor)

		Ek = Ek * hartree_wavenumber
		VAve = VAve * hartree_wavenumber
		UAve = UAVe * hartree_wavenumber
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
      write(*,*) 'Total Energy = ', Etot  , 'cm-1.' 

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
      subroutine derivs(x,Ntraj,fx,v)
      use pH2 
      implicit real*8(a-h,o-z)
      integer*4,intent(in) :: ntraj
      real*8,   intent(in) :: x(Ntraj)
      real*8,   intent(out) :: v(ntraj),fx(ntraj)
      real*8 :: dv(ntraj)
      character*20 :: PES 
     
      PES = 'pH2' 
      
      if(PES == 'Morse') then 
 
      a0 = -0.180516
      a1 = 0.196123
      a2 = 0.95861
      a3 = 1.33249
      a4 = 0.0322998

      do i=1,ntraj
        z = 1d0-exp(-a2*(x(i)-a3))
        v(i)= a0 + a1*z**2 + a4/x(i)**4
        dv(i) = 2d0*a1*z*a2*x(i)*exp(-a2*(x(i)-a3)) - 4d0*a4/x(i)**5    
        fx(i) = -dv(i) 
      enddo 

      elseif(PES == 'pH2') then 
      
        dr = 1.0d-4 
            
        do i = 1, Ntraj 

           r = x(i) + re/a2bohr 

           v(i) = vpot(r) 
           dv(i) = (vpot(r+dr)-vpot(r))/dr    ! finite difference for force      
         enddo 
         
         fx = -dv 
             
      endif 
    
      return 
      end subroutine
      
      
	double precision function vpot(r) 
	use pH2 
	implicit real*8(a-h,o-z) 
	
	r = r * a2bohr

        beta_inf = log(2.0 * De / u_LR(re)) 
           s = 0.0 
          do j = 0,10 
            s = s + b(j) * y_ref(r,1)**j
          enddo

          
          beta = y_ref(r,6) * beta_inf + (1.0-y_ref(r,6))  * s  
    
    vpot = De*(1.0 - u_LR(r)/u_LR(re) * exp(-beta * y_eq(r,6)))**2
    
	vpot = vpot + Vmin 
    
    vpot = vpot / hartree_wavenumber 
          
	end  
	

    double precision function y_eq(r,n)
    use pH2, only : re 
    implicit real*8(a-h,o-z)
     
    y_eq = (r**n - re**n)/(r**n + re**n) 
 
    end function 
    
    double precision function y_ref(r,n)
    use pH2, only : r_ref  
    implicit real*8(a-h,o-z)
     
    y_ref = (r**n - r_ref**n)/(r**n + r_ref**n) 
 
    end function     
    
    double precision function u_LR(r)
    use pH2, only : C6, C8, C10, damp 
    implicit real*8(a-h,o-z)
	
	u_LR = damp(r,6) * C6/r**6 + damp(r,8) * C8/r**8 + damp(r,10) * C10 / r**10 
      
    end function 
    
