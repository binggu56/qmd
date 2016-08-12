	program ftr
	implicit doubleprecision (a-h,o-z)
        integer*4,parameter :: nf=20000
	real*8, dimension(nf) :: w, t
        complex*16, dimension(nf) :: ac,ft 
	complex*16 im,te
	real*4 a1,a2,a3,a4,a5
        character*32 fil1,fil2,fil3
	
        im=(0.d0,1.d0)

	open(20,file='ftr_in',status='old')
	read(20,*) nm,beta
	read(20,*) fil1

        open(25,file='cor1', status='old')
        open(26,file='fft.dat')
!       energy  units
        eu=27.2d0/1836.15
        emx=4d0
        write(*,*) 'energy unit',eu

	do 1 i=1,nm
	  read(25,1000) x1,ac(i),xx
	 t(i)=x1
1	continue

        print *,'ac',ac(1) 

	ac(1)=ac(1)/2.d0
	ac(nm)=ac(nm)/2.d0

!-------- damp the signal ----------------
!        beta = 0d0 
!        beta=-log(beta/x2)/t(nm)**2

!       do  i=1,nm
!         ac(i)=ac(i)*exp(-beta*t(i)**2)
!        enddo

	pi=4.d0*atan(1.d0)
	dt=abs(t(2)-t(1))
	write(*,*) dt,nm
	h=0.5d0/nm/dt*pi
	nw=emx/h
	write(*,*) h,nw

        de = 0.05
        nw = 100
	do 3 i=1,nw
	  w(i)=de*i
3	continue

        print *,im,ac(1),dt,t(1),w(1)
        
	do 2 k=1,nw
	  ft(k)=(0.d0,0.d0)
	  do  i=1,nm
	    ft(k) = ac(i)*exp(im*w(k)*t(i)) + ft(k)
          enddo
!	  fin=dreal(ft(k))*dt
	  fin=abs(ft(k))*dt
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          print *,'ft',ft(k)
          write(26,1000) w(k),ft(k)*dt
2	  continue

1001    format(20(e14.5,1x))
1000    format(20(e14.7,1x))
	stop
	end
