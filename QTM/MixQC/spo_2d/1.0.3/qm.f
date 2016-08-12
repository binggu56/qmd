! ------------------------
      program qm
      use mass 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc Coupled Morse oscillators ccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	implicit real*8 (a-h,o-z)
	integer*4 ni,nj
	parameter (ni=128, nj=512)
	real*8 ak2(ni,nj),norm,akx(ni,nj),aky(ni,nj)
	complex*16 v(ni,nj)
	complex*16 im,psi(ni,nj),hpsi(ni,nj)
	complex*16 psi0(ni,nj),c1,c2,po,ki
	real*8 :: corr,wx,wy,hbar
      integer*4 ipot
	common /params/ d,x0,y0,di,wi
	common /para/ hbar
	common /para1/wx,wy,corr
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
	common /ini/ qx0,qy0,px0,py0
	common /wav/ ax,ay
      common /pes/ ipot
      common /four/emin,emax,DE
C
C
	im=(0.d0,1.d0)
	pi=4.d0*atan(1.d0)
	N2=ni*nj
	d=160d0
	X0=0d0
	y0=0d0
c	hbar=0.06d0
c	z0=1.04435d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	open(15,file='cor1')
	open(16,file='cor2')
	open(17,file='wf0')
	open(18,file='wft')
	open(19,file='xav')
	open(20,file='en.dat')
      open(21,file='cwf')
      open(22,file='cor')
      open(23, file = 'norm.dat') 
      open(120,file='ps.dat')
	open(5,file='IN')
c time propagation parameters
	read(5,*) dt,kmax,nout
	read(5,*) mx,my 
c initial wave packet parameters
	read(5,*) ax,ay
	read(5,*) qx0,qy0
	read(5,*) px0,py0
c potential surface parameters
        read(5,*) DI,WI
	read(5,*) xmin,ymin
	read(5,*) xmax,ymax
	read(5,*) wx,wy
	read(5,*) corr
      read(5,*) ipot
	close(5)
c inverse mass prefactor 1/2/m
c hydrogen
c deuterium m=2*mh
c	am=0.5d0
c	alfa=alfa*sqrt(2.d0)
cc tritium=3*mh
c	am=1.d0/3.d0
cc	alfa=alfa*sqrt(3.d0)
cc
C	ax=wx*sqrt(mx)/(2d0)
C	ay=wy*sqrt(my)/(2d0)
	dx=(xmax-xmin)/ni
	dy=(ymax-ymin)/nj
	ds=dx*dy
      
      hbar = 1d0 

      emin = 0.001
      emax = 4d0
      de = 0.001d0
c--------------------------------------------------
c 	print out initial conditions
c--------------------------------------------------
        write(*,*)'initial state'
	print *,'constants for initial wavefunction'
	print *,'ax =',ax
	print *,'ay =',ay
	print *,'initial center of  wavepacket'
	print *,'qx0 =',qx0
	print *,'qy0 =',qy0
	print *,'initial momentum'
	print *,'px0 =',px0
	print *,'py0 =',py0
	print *,'grid size'
	print *,'xmin =',xmin
	print *,'ymin =',ymin
	print *,'xmax =',xmax
	print *,'ymax =',ymax
	print *,'dx =',dx
	print *,'dy =',dy
	print *,'mass'
	print *,'mx =',mx
	print *,'my =',my
c----------------------------------------
c 	define vector ak2
c----------------------------------------
	consI=2.d0*hbar*pi/dx/ni
	consJ=2.d0*hbar*pi/dy/nj
	ni2=ni/2
	nj2=nj/2
	do j=1,nj
	 akj=(j-1)*consj
	 if (j.gt.nj2) akj=-(nj+1-j)*consj
	 do  i=1,ni
	  aki=(i-1)*consi
	  if (i.gt.ni2) aki=-(ni+1-i)*consi
	akx(i,j)=(aki*aki/(2d0*mx))
	aky(i,j)=(akj*akj/(2d0*my))
	ak2(i,j)=aki*aki/(2d0*mx)+akj*akj/(2d0*my)
        enddo
        enddo
c------------ set up the initial wavefunction -----------------------
      an0=hbar*pi/2d0/dsqrt(ax*ay)
      an0=1d0/dsqrt(an0)
	
	ch_n=0d0
	 do j=1,nj
         y=ymin+dy*(j-1)
          do i=1,ni
           x=xmin+dx*(i-1)
           psi(i,j)= exp(-ax*(x-qx0)**2/hbar-ay*(y-qy0)**2/hbar
     &             + im*px0*(x-qx0)/hbar+im*py0*(y-qy0)/hbar)
          psi(i,j) = an0*psi(i,j)

	   !if(abs(psi(i,j)).gt.1d-3) 
	 write(17,*) x,y,abs(psi(i,j))
         ch_n=ch_n+conjg(psi(i,j))*psi(i,j)*ds
	enddo
	write (17,*) ' '
	enddo

      write(*,*) 'NORM=',ch_n
      
      psi0 = psi
cccccccccccccccccc propagate psi ccccccccccccccccccccccccccccccccc
	call ham(ni,nj,N2,ak2,akx,aky,psi,v,hpsi,en0,enx,eny,po,ki)
	call correl(ni,nj,psi0,psi,c1,c2)
	t=0.d0
	write(15,1000) t,c1,abs(c1)
	write(16,1000) 2d0*t,c2,abs(c2)
      print *,'proceeding........'
	call split(kmax,nout,ni,nj,N2,dt,t,v,ak2,akx,aky,psi0,psi)
c------------ print the final  wavefunction -----------------------
c	norm=0d0
	do j=1,nj
         y=ymin+(j-1)*dy
	 do i=1,ni
	  x=xmin+dx*(i-1)
	  if(abs(psi(i,j)).gt.1d-3) then 
c	norm=norm+conjg(psi(i,j))*psi(i,j)*ds
	    write(18,*) x,y,abs(psi(i,j))
          endif 
	enddo 
	write (18,*) ' '
	enddo
c	print *,'Norm_Final =',norm
	
c--------------------------------------------------------------------

1000  format(1000(e14.7,1x))
	stop
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine ham(ni,nj,N2,ak2,akx,aky,psi,v,hpsi,
     &   en0,enx,eny,po,ki)
	implicit real*8(a-h,o-z)
	integer*4,intent(IN) :: ni,nj
	complex*16 psi(ni,nj),hpsi(ni,nj),tpsi(ni,nj),vpsi(ni,nj),
     & aux(ni,nj),auy(ni,nj),vxpsi(ni,nj),vypsi(ni,nj)
	real*8,intent(IN):: ak2(ni,nj),akx(ni,nj),aky(ni,nj)
	complex*16 v(ni,nj),c,enkx,enky,po,ki
        real*8,intent(OUT)::en0,enx,eny
	real*8 :: evx,evy
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	call potpsi(ni,nj,v,psi,vpsi)
	call aver(ni,nj,psi,vpsi,po)
c	en0=dreal(po)
c	call aver(ni,nj,psi,vxpsi,envx)
c	call aver(ni,nj,psi,vypsi,envy)
c	write(*,*) 'potential energy',C
	
	call kinpsi(ni,nj,N2,ak2,psi,tpsi)
	call aver(ni,nj,psi,tpsi,ki)
	call kinpsi(ni,nj,N2,akx,psi,aux)
	call aver(ni,nj,psi,aux,enkx)
	call kinpsi(ni,nj,N2,aky,psi,auy)
	call aver(ni,nj,psi,auy,enky)
	call xex2(ni,nj,psi,evx,evy)
c	write(*,*) 'kinetic energy', C
	en0=dreal(ki+po)
	enx=dreal(enkx)+evx
	eny=dreal(enky)+evy
C	write(*,*) 'TOTAL E', en0
	do  j=1,nj
	 do  i=1,ni
	  hpsi(i,j)=tpsi(i,j)+vpsi(i,j)
	enddo
	enddo
	return
	end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine potpsi(ni,nj,v,psi,vpsi)
	implicit real*8(a-h,o-z)
	integer*4,intent(IN) :: ni,nj
      integer*4 ipot
	real*8 corr,wx,wy,De,Xe,k2
      complex*16,intent(IN)::psi(ni,nj)
      complex*16::v(ni,nj),vx(ni,nj),vy(ni,nj)
	complex*16,intent(OUT)::vpsi(ni,nj)
C
	common /params/ d,x0,y0,di,wi
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
	common /para1/wx,wy,corr
      common /pes/ ipot
c------------- set up potetial -------------------------------------------
c--------------------------------------------
      k2=0.3696726511d0
      De=160d0
      Xe=1.4d0
      epi = 0.d0


      DO j=1,nj
         y=ymin+(j-1)*dy
      DO i=1,ni
         x=xmin+(i-1)*dx
C Morse
      IF(IPOT .EQ. 1) THEN
      v(i,j)=De*(1d0-exp(-(x-xe)))**2+y**2/2d0
      vpsi(i,j)=psi(i,j)*v(i,j)
C Eckart
      ELSEIF(IPOT .EQ. 2) THEN
      v(i,j)=0.063740d0-0.060553d0/(1d0+exp(-1.5d0*x))-
     & 0.0981148d0/(cosh(0.75d0*x))**2+k2*(y-1.40d0)**2/2d0
      vpsi(i,j)=psi(i,j)*v(i,j)
C Henon-Hails
      ELSEIF(IPOT .eq. 3) THEN
        v(i,j)=(x**2 + y**2)/2d0 + epi*x*y

! --- reactive scattering 
      elseif(ipot == 4) then 
      
      ak1 = 5.0 
      ak2 = 15.0 
      b = 0.5d0 
      v0 = 16d0 
      g = 1.3624d0 
      
      ak = 0.5*(ak1+ak2)+0.5*(ak2-ak1)*tanh(b*y)
      v(i,j) = 0.5*ak*x**2 + v0/cosh(g*y)**2
c        vx(i,j)=wx**2*x**2/2d0
c        vy(i,j)=wy**2*y**2/2d0
c          if(vr.gt.5d0*d) vr=5d0*d
c          vi=0d0
c          r=dsqrt(rx*rx+ry*ry)
c          if(r.gt.wi) vi=vi+di*(r-wi)**2
c          v(i,j)=vr-im*vi
        vpsi(i,j)=psi(i,j)*v(i,j)
c        vxpsi(i,j)=psi(i,j)*vx(i,j)
c        vypsi(i,j)=psi(i,j)*vy(i,j)
      ENDIF
	enddo
	enddo
        return
        end subroutine


        

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine kinpsi(ni,nj,N2,ak2,cx1,aux)
        implicit real*8(a-h,o-z)
c	parameter(nij=512)
	integer*4,intent(IN) :: ni,nj
        dimension ak2(ni,nj)
        complex*16 cx1(ni,nj),aux(ni,nj)
        dimension nn(2)
        nn(1)=ni
        nn(2)=nj
        do 10 j=1,nj
         do 10 i=1,ni
10      aux(i,j)=cx1(i,j)
        call fourn(aux,nn,2,1)
        do 11 j=1,nj
         do 11 i=1,ni
11      aux(i,j)=aux(i,j)*ak2(i,j)
        isign=-1
        call fourn(aux,nn,2,-1)
        do 12 j=1,nj
         do 12 i=1,ni
12      aux(i,j)=aux(i,j)/N2
        return
        end

c-------------------------------------------------------
c-------------------------------------------------------
      subroutine split(kmax,nout,ni,nj,N2,dt,t,v,ak2,akx,
     &  aky,psi0,psi)
      use mass 
      use sci 

      implicit real*8(a-h,o-z)
      integer*4,intent(IN)::ni,nj
      dimension ak2(ni,nj),akx(ni,nj),aky(ni,nj)
      dimension nn(2)
      complex*16:: psi(ni,nj),v(ni,nj),po,ki,cor(kmax)
      complex*16 psi0(ni,nj),hpsi(ni,nj),c1,c2,out(kmax)

      real (kind = 8), dimension(ni,nj) ::  dphi,phi 
      real (kind = 8) fd4(5) 

	common /grid/ xmin,ymin,xmax,ymax,dx,dy
      common/correlation/cor_d
      common/four/emin,emax,de
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open(23,file = 'den.dat' )
      open(24,file = 'norm.dat')
 
      dt2=dt*0.5
	n2=ni*nj
      nn(1)=ni
      nn(2)=nj


      do 11 k=1,kmax
	  t=t+dt
      
! --- nonlinear operator exp(-i*dt*L) 
! --- subtract quantum potential of CDOF 
        
        call cqpot(ni,nj,psi,dt)

         call fourn(psi,nn,2,1)
         call diff(psi,ak2,dt2,ni,nj)
         call fourn(psi,nn,2,-1)
         do 12 j=1,nj
          do 12 i=1,ni
12       psi(i,j)=psi(i,j)/N2
 
         call phase(psi,v,dt2,ni,nj)
         call fourn(psi,nn,2,1)
         call diff(psi,ak2,dt2,ni,nj)
         call fourn(psi,nn,2,-1)
         do 13 j=1,nj
          do 13 i=1,ni
13       psi(i,j)=psi(i,j)/N2
      call phase(psi,v,dt2,ni,nj)

! --- observables             
      
      call correl(ni,nj,psi0,psi,c1,c2)
      cor(k) = c1

      call ham(ni,nj,N2,ak2,akx,aky,psi,v,hpsi,en0,enx,eny,po,ki)
!	call xex(ni,nj,psi,xav,yav,xx,yy,ex)
!	call cor(ni,nj,psi0,psi,auto,au_x,au_y)
!      write(22,1000) t,cor_d
      write(20,1000) t,po,ki,en0 
	write(15,1000) t,c1,abs(c1)
!	write(16,1000) 2d0*t,c2,abs(c2)
!	write(19,1000) t,xav,yav,xx,yy,auto,au_x,au_y
cccccc 	print out the wavefunction for every time step
C	do s=1,100
C	if(k==s*kmax/100) then
C	do j=1,nj
C		y=ymin+(j-1)*dy
C	do i=1,ni
C		x=xmin+(i-1)*dx
C	write(21,*) x,y,abs(psi(i,j))
C	enddo
C	write(21,*) ' '
C	enddo
C	write(21,*) ' '
C	write(21,*) ' '
C	endif
C	enddo

! --- normalization 
      a1 = 0d0 
	do j=1,nj
	  y=ymin+(j-1)*dy
	do i=1,ni
	  x=xmin+(i-1)*dx
	  a1 = a1 + abs(psi(i,j))**2*dx*dy 
      enddo 
      enddo 
      write(100,1000) t,a1  

      ak1 = 15d0 
!      ak2 = 15d0 
      pi = 4d0*atan(1.0) 
      alfa = dsqrt(ak1)/2d0
      anrm = dsqrt(alfa/pi) 

      d0 = 0d0 
      d1 = 0d0 
      do j=1,nj
        y=ymin+(j-1)*dy
      if( y > 0d0) then 
      do i=1,ni
        x=xmin+(i-1)*dx
        d0 = d0+abs(psi(i,j))**2*exp(-alfa*x**2)*dy*dx  
        d1 = d1+abs(psi(i,j))**2*x**2*exp(-alfa*x**2)*dx*dy 
      enddo 
      endif 
      enddo 

      d0 = d0*anrm 
      d1 = d1*dsqrt(alfa/pi)*(2d0*alfa) 

      write(23,1000) t,d0,d1 
      

11	enddo

      close(23) 

      call fourier(cor,kmax,dt,out)
      
      i=1
      do while (Emin+(i-1)*De<Emax .and. i<=kmax)
        en = Emin + (i-1)*de
        write(120,1000) en, out(i) 
        i = i+1
      enddo

1000  format(1000(e14.7,1x))
15	return
	end

!-----fourier transform----------------------
      subroutine fourier(data,Nt,dt,out)
      implicit real*8(a-h,o-z)
      integer*4, intent(in) :: nt
      complex*16 :: im,data(Nt),out(Nt)
      common/four/Emin,Emax,de

      im = (0d0,1d0)
      Pi = 4.0d0*atan(1d0)
!      alfa = 1d0/2d0**4
      alfa = 0d0

      out = (0d0,0d0)

        i = 1
        do while (Emin+(i-1)*De<Emax.and.i<=Nt)
        en = Emin + (i-1)*de
        do j=1,Nt
          t = 0d0 +j*dt
          out(i) = 1d0/2d0/Pi*2d0*real(data(j)*exp(im*en*t))*
     +             exp(-alfa*t**2)*dt + out(i)
        enddo
        i=i+1
        end do

        return
        end subroutine
!------------------------------------------------
c-------------------------------------------------------
	subroutine diff(cwf,ak2,ts,ni,nj)
	implicit real*8(a-h,o-z)
	integer*4,intent(IN)::ni,nj
	complex*16 nim,cwf(ni,nj)
	dimension ak2(ni,nj)
	real*8 :: hbar
	common /para/ hbar
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	nim=(0.d0,-1.d0)
	do 11 j=1,nj
	 do 11 i=1,ni
11	cwf(i,j)=cwf(i,j)*exp(nim*ts*ak2(i,j)/hbar)
	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine phase(cwf,v,ts,ni,nj)
	implicit real*8(a-h,o-z)
	complex*16 nim,cwf(ni,nj)
	complex*16 v(ni,nj)
	real*8 :: hbar
	common /para/ hbar
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	nim=(0.d0,-1.d0)
	do 11 j=1,nj
	 do 11 i=1,ni
	if (ts.ge.0) cwf(i,j)=cwf(i,j)*exp(nim*ts*v(i,j)/hbar)
11	if (ts.lt.0) cwf(i,j)=cwf(i,j)*exp(nim*ts*conjg(v(i,j))/hbar)
	return
	end

c-------------------------------------------------------
c-------------------------------------------------------
	subroutine correl(ni,nj,psi0,psi,c1,c2)
	implicit real*8(a-h,o-z)
      integer*4,intent(IN)::ni,nj
	complex*16 psi0(ni,nj),psi(ni,nj),c1,c2
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
ccc---------- compute correlation functions ------------------------
	c1=(0d0,0d0)
	c2=(0d0,0d0)
	
      do 10 j=1,nj
	 y=ymin+(j-1)*dy
	  do 10 i=1,ni
	    x=xmin+(i-1)*dx
          c1=c1+conjg(psi0(i,j))*psi(i,j)*dx*dy
          c2=c2+psi(i,j)**2*dx*dy
10	continue
C        c1=c1*dx*dy
C        c2=c2*dx*dy
	return
	end subroutine

c-------------------------------------------------------
c-------------------------------------------------------
        subroutine aver(ni,nj,psi0,psi,c)
        implicit real*8(a-h,o-z)
c        parameter(ij=512)
        integer*4,intent(IN)::ni,nj
        complex*16 c,psi0(ni,nj),psi(ni,nj)
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
	c=(0.d0,0.d0)
	ct=0.d0
        do 11 j=1,nj
         do 11 i=1,ni
	ct=ct+abs(psi0(i,j))**2
11      c=c+psi(i,j)*conjg(psi0(i,j))*dx*dy
        return
        end

      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL*8 data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL*8 tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END


      subroutine cqpot(ni,nj,psi,dt) 
      use mass
      use sci  

      implicit real (kind = 8) (a-h, o-z) 
      real (kind = 8), dimension(ni,nj) :: psi,phi,dphi 
      real (kind = 8) fd4(5) 
 
      common/grid/xmin,ymin,xmax,ymax,dx,dy       
      
      fd4 = (/-1.0,16.0,-30.0,16.0,-1.0/) ! fourth-order finite-difference coefficients
      dt2 = dt/2d0 
      del = 1d-8 

      do j=1,nj 
      do i=1,ni       
        phi(i,j) = abs(psi(i,j))
      enddo 
      enddo 

      dphi = 0d0 
      do j=3,nj-2 
        do i=1,ni  
        dphi(i,j) = (fd4(1)*phi(i,j+2)+fd4(2)*phi(i,j+1)+
     &       fd4(3)*phi(i,j)+fd4(4)*phi(i,j-1)+
     &       fd4(5)*phi(i,j-2))/12d0/dy**2
        enddo 
      enddo 

! ---- boundary points, treat out-of-boundary points as 0 
      j = nj-1 
      do i=1,ni 
        dphi(i,j) = (fd4(2)*phi(i,j+1)+
     &       fd4(3)*phi(i,j)+fd4(4)*phi(i,j-1)+
     &       fd4(5)*phi(i,j-2))/12d0/dy**2
      enddo 
      
      j = nj
      do i=1,ni   
        dphi(i,j) = (fd4(3)*phi(i,j)+fd4(4)*phi(i,j-1)+
     &       fd4(5)*phi(i,j-2))/12d0/dy**2
      enddo 
      
      j = 2 
      do i=1,ni 
        dphi(i,j) = (fd4(1)*phi(i,j+2)+fd4(2)*phi(i,j+1)+
     &       fd4(3)*phi(i,j)+fd4(4)*phi(i,j-1))/12d0/dy**2
      enddo 
      
      j = 1
      do i=1,ni  
        dphi(i,j) = (fd4(1)*phi(i,j+2)+fd4(2)*phi(i,j+1)+
     &       fd4(3)*phi(i,j))/12d0/dy**2
      enddo


!      do i=1,ni  
!        dphi(i,1) = (phi(i,3)+phi(i,1)-2d0*phi(i,2))/dy**2
!      enddo 
!      do i=1,ni  
!        dphi(i,nj) = (phi(i,nj)+phi(i,nj-2)-2d0*phi(i,nj-1))/dy**2
!      enddo 

      do j=1,nj 
        do i=1,ni 
          psi(i,j) = psi(i,j) - im*dt*psi(i,j)*dphi(i,j)/
     &               (2d0*my*(phi(i,j)+del))
        enddo 
      enddo 

      return 
      end subroutine 
