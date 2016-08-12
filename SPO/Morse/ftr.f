	program ftr
	implicit doubleprecision (a-h,o-z)
	dimension w(20000),t(20000)
	complex*16 ac(20000),ft(20000),ci,te
	real*4 a1,a2,a3,a4,a5
	character*32 fil1,fil2,fil3
	ci=(0.d0,1.d0)
	open(20,file='ftr_in',status='old')
	read(20,*) nm,beta
	read(20,*) fil1,fil2
        open(25,file=fil1,status='old')
        open(26,file=fil2)
c       energy  units
	eu=27.2d0/1836.15
	emx=50d0
	write(*,*) 'energy unit',eu

	do 1 i=1,nm
	  read(25,1000) x1,ac(i),x2
	 t(i)=x1
1	continue
	ac(1)=ac(1)/2.d0
	ac(nm)=ac(nm)/2.d0
c-------- damp the signal ----------------
        beta=-log(beta/x2)/t(nm)**2
        do 21 i=1,nm
21       ac(i)=ac(i)*exp(-beta*t(i)**2)

	pi=4.d0*atan(1.d0)
	dt=abs(t(2)-t(1))
	write(*,*) dt,nm
	h=0.5d0/nm/dt*pi
	nw=emx/h
	write(*,*) h,nw

	do 3 i=1,nw
	  w(i)=h*i
3	continue

	do 2 k=1,nw
	  ft(k)=(0.d0,0.d0)
	  do 4 i=1,nm
	    te=exp(ci*w(k)*t(i))
4	    ft(k)=ac(i)*te+ft(k)
c	  fin=dreal(ft(k))*dt
	  fin=abs(ft(k))*dt
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          write(26,1000) w(k),fin,ft(k)*dt
2	  continue

1001    format(20(e14.5,1x))
1000    format(20(e14.7,1x))
	stop
	end
