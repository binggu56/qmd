! one dimensional fourier transform
	program fourier
	implicit none 
	integer*4 Nt,i,j
	real*8 dt,Emax,x,x0,p0,a0,Pi,t,gy
	real*8 ymin,y
	complex*16,allocatable :: data(:),out(:)
	complex*16 :: im
	real*8,allocatable :: E(:)
!	complex*16,allocatable :: data(:)

	open(100,file='OUT')
	open(5,file='INf')
	read(5,*) Nt
	read(5,*) dt
	close(5)
!	Nt = 200
!	dt = 0.01d0
	Emax = 400d0
	im = (0d0,1d0)
	Pi = 4.0d0*atan(1d0)

	x0 = -5d0
	p0 = 12d0
	a0 = 10d0
!	gy = 2d0*pi/(Nt*dt)
	gy = 1d0/(nt*dt)
	ymin = -10d0

	allocate(data(Nt),out(Nt),E(Nt))
!	open(6,file='corrf')
!	do i=1,Nt
!	read(6,10000) data(i)
!	enddo
!	close(6)

	out = (0d0,0d0)
!
	do i=1,Nt
	y = gy*(i-1) + ymin
	do j=1,Nt
	x = 0d0 +(j-1)*dt
!	data(j) = exp(-a0*(x-x0)**2+im*p0*(x-x0))
	data(j) = exp(-a0*x**2)
	out(i) = 1d0/2d0/Pi*data(j)*exp(im*y*x)*dt + out(i)
	enddo
	write(100,10000) y,out(i)
	enddo
	
10000	format(1000(e14.7,1x))
 	end program fourier
	
	
	
	
	
