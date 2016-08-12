! one dimensional fourier transform
	subroutine fourier(data,Nt,dt,out)
	implicit none 
	integer*4 Nt,i,j
	real*8 dt,Emax,x,x0,p0,a0,Pi,t,gy
	real*8 ymin,y,en,emin,de
	complex*16 :: data(Nt)
	complex*16,intent(OUT) :: out(Nt)
	complex*16 :: im
	
	common/four/emin,de
	
	Emax = 400d0
	im = (0d0,1d0)
	Pi = 4.0d0*atan(1d0)

!	de = 1d0/(Nt*dt)
	De = 1.d-8
	emin = 0d0

	out = (0d0,0d0)

	do i=1,Nt
	en = de*(i-1) + emin
	do j=1,Nt
	t = 0d0 +j*dt
	out(i) = 1d0/2d0/Pi*data(j)*exp(im*en*t)*dt + out(i)
	enddo
	enddo
	
	return
 	end subroutine
	
	
	
	
	
