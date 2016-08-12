! one dimensional fourier transform, to use this subroutine 
! 1.define {emin,de} and make them common in the main program
! 2.define an output array for {out}
	program main 

	implicit real*8(a-h,o-z) 
	integer*4 Nt,i,j
	real*8 dt,Emax,x,x0,p0,a0,Pi,t(nt)
	real*8 ymin,y,en,emin,de
	complex*16 :: data(Nt)
	complex*16,intent(OUT) :: out(Nt)
        complex*16 :: im
        
        open(100,file='fft.dat',status='new')

	open(10,file='corr')
        do i=1,nt 
          read(10,*) t(i), data(i) 
        enddo 

	Emax = 0.1d0
	im = (0d0,1d0)
	Pi = 4.0d0*atan(1d0)

        de = 1d0/(Nt*dt)
        emin = -0.18d0
        out = (0d0,0d0)

	i = 1
	do while (Emin+(i-1)*de < Emax)
        
        en = Emin + (i-1)*de
!	if(en > Emax) goto 100
        do j=1,Nt
          t = 0d0 +j*dt
!	data(j) = exp(-a0*(x-x0)**2+im*p0*(x-x0))
          out(i) = data(j)*exp(im*en*t)*dt + out(i)
        enddo

        i = i+1

        write(100,1000) y,out(i)
        end do
        
1000    format(1000(e14.7,1x))

        return 
        end program 


        
        
        
