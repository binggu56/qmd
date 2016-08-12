        program main
c
ccccccccc real H and the convergence factor ccccccccccccccc
ccccccc FFT and split propgation 
ccccccc Number of the grid points has to be a power of 2 ccccccc
c
	implicit real*8(a-h,o-z)
	parameter(nx=2048)
        REAL*8 AK2(nx),x(nx)
	complex*16 im,c1,c2,c4,c3,v(nx),HPSI(nx),VPSI(nx)
     & ,psia(nx),psib(nx),psi0(nx),psi(nx)
!     	complex*16 test
C	
	integer*4 :: nn(1) 
	real*8 xmin,xmax,dx,prob
	real*8 mx,k,di,wi,vr,vi
	common/grid/xmin,xmax,dx
	common/four/emin,de
	complex*16 :: sab,eta_a,eta_b
	complex*16,allocatable :: dc(:),f_c(:),c(:)
	
C	common/prob/prob
cc I2	parameter(a10=0.98267d0,a11=0.05716d0,a12=5.04d0,ami2=1.165d4
c     &	,z1=.2047033605d0,z2=4.82456958d0)
cc hyd	parameter(a10=1.0144d0,a11=0.1744d0,a12=5.0d0,ami2=911.4d0)
C	ami2=1d0
C**********************************************************
C
	open(44,file='pot')
        open(54,file='wf0')
        open(55,file='wft')
        open(53,file='corr')
        open(56,file='rho')
        open(57,file='prob')
	open(59,file='fourc')
	open(60,file='corrf')
	OPEN(101,FILE='smat')
	open(102,file='corr1')
!	OPEN(103,file='test')
        OPEN(15,FILE='IN',STATUS='OLD')
c grid and  wavepacket
        READ(15,*) N1,xmin,xmax
	read(15,*) a0
	read(15,*) x0,p0
c potential 
	read(15,*) di,wi
c number of time steps and time increment 
        READ(15,*) Nt,dt,nout,ipot
	read(15,*) mx
!	read(15,*) Depth
C
	write(*,*) 'pot',di,ai
	print *,xmin,xmax
        write(*,*) 'wp',N1,a0,x0,p0
        write(*,*) 'time',Nt,dt,nout,ipot
	
	allocate(dc(Nt),c(Nt),f_c(Nt))

C       CONSTANTS TO SAVE COMPUTAION TIME
!        mx = (m1*(2d0*m2)/(m1+2d0*m2))*1836d0
	im=(0.d0,1.d0)
        AM=.5d0/mx
	pi=4d0*atan(1d0)
	an=dsqrt(dsqrt(2d0*a0/pi))
	al = a0

cccccccccccccccccccccccccc define the grid cccccccccccccccccccccccccccccc
C        xmin=-xmax
cc WARN for h2
cc	xmin=3.5d0
cc	if(ipot.eq.0) then
cc	 xmin=z2-xmax
cc	 xmax=z2+xmax
cc	endif
        xn=xmax-xmin
        dx=xn/(n1-1)
        write(*,*) 'dx=',dx
        do 10 i=1,N1
          x(i)=xmin+dx*(i-1)
C	print *,x(i)
10      continue
!	write(*,*) 'xfin',x(n1),xmax
!	write(*,*) 'bug,al,x0,p0',a0,x0,p0

c	do 777 ip0=1,30
c	 e0=2d0+ip0*1d0
c	 p0=dsqrt(2d0*e0)
	 csave=0d0
cccccccccccccccccccccccccc define the potential ccccccccccccccccccccccccc
        do 11 i=1,N1
	xi = xmin+(i-1)*dx

C 	ECKART BARRIER
	if(ipot .eq. 1) then
	        D = 16d0
        	a = 1.3624d0
		vi = 0d0
		Vr = D/cosh(a*xi)**2
		if ( abs(xi) .gt. wi) then 
			vi=vi+di*(xi-wi)**2
		endif
		v(i) = vr -im*vi
	endif
cccccccc harmonic ccccccccccccccc
c approx morse
	 if(ipot.eq.0) v(i)=z1*(x(i)-z2)**2
cccccc parabolic  cccccc
	  if(ipot.eq.-1) v(i)=-d*d*x(i)*x(i)/2d0
cccccc  double well cccccccccccccc
	if(ipot.eq.2) v(i)=x(i)**2*(aw*x(i)**2-bw)
cccccccc Morse oscillator cccccccccccccccccccccccc
	if(ipot.eq.3) then
	  t5=exp(-a10*(x(i)-a12))
	  v(i)=a11*(1.d0-t5)**2
	  if(abs(v(i)).gt.5d0*a11) v(i)=5.d0*a11
	endif
cccccccccccccccc metastable well ccccccccccccccccc
	if(ipot.eq.4) then
	 x2=x(i)*x(i)
	 v(i)=x2*(1.d0-0.1d0*x(i)+0.1d0*x2)
	endif
11	continue
cccccccccccccccc  initial wavepacket cccccccccccc

	do i=1,N1
	  psi(i)=an*exp(-a0*(x(i)-x0)**2+im*p0*(x(i)-x0))
	  psi0(i)=psi(i)
	  psib(i)=an*exp(-a0*(x(i)+x0)**2+im*p0*(x(i)+x0))
	  write(44,1000) x(i),v(i)
	enddo

	do i=1,n1
	 if(abs(psi(i)).gt.1d-3) write(54,*) x(i),abs(psi(i))
	enddo

        PI=4.d0*atan(1.d0)
        TWPI=2.d0*PI
        CONS=TWPI/(DX*N1)
        N2=N1/2+1
C       THIS DO LOOP DEFINES THE VECTOR AKX2
        DO 1 N=1,N1
        AK=(N-1)*CONS
        IF(N.GT.N2) AK=-(N1+1-N)*CONS
 1      AK2(N)=AK*AK*AM
C
	c0=(0.d0,0.d0)
	do 13 j=1,n1
13	  c0=c0+abs(psi(j))**2*v(j)*dx
	
        WRITE(*,*) 'POT ENERGY=',c0
        call avec(N1,dx,x,psi0,psib,psi,c1,c2,c0,as,cm1,cm2)
	write(56,1000) 0d0,c0,as,as/cm1-(c0/cm1)**2,cm1
	write(57,1000) 0.d0,cm1,cm2,abs(c1),abs(c2),c0,as
	write(53,1000) 0.d0,c1,c2
!	write(60,1000) 0.d0,c2
C       THE MAIN TIME INTEGRATION LOOP
	dt2=dt/2.d0
	c = 0d0

        DO 2 J=1,NT
         do jj=1,nout
          call split(n1,dt,dt2,v,ak2,psi)
         enddo

         call avec(N1,dx,x,psi0,psib,psi,c1,c2,c0,as,cm1,cm2)
         t=dt*j*nout
	
	do i=1,N1
	  c(j) = c(j) + psi(i)*conjg(psib(i))*dx
	enddo
	write(57,1000) t,cm1,cm2,abs(c1),abs(c2),c0,as
	write(53,1000) t,c1,c2
	write(56,1000) t,c0,as,as/cm1-(c0/cm1)**2,cm1
C REACTION PROBALITY
!        prob=0d0
!         do i=1,n1
C          x=xmin+dx*(i-1)
!        if( x(i) < -5.d0) prob=prob+abs(psi(i))**2*dx
C       endif
!        enddo
!        WRITE(103,1000) t,c12
C
       if(c0.gt.csave) csave=c0
2      CONTINUE
	
!	do i=1,N1
	
	do i=1,Nt
	write(102,1000) i*1.d0,c(i)
	enddo
	call fourier(c,Nt,dt,f_c)
	do i=1,Nt
	en = emin + (i-1)*de
        write(59,1000) en,f_c(i)
        enddo
	
	do i=2,Nt

	en = emin + (i-1)*de
	k = dsqrt(2d0*mx*en)
	const = an*dsqrt(mx/(2d0*k*a0))
! fourier transform of initial wavepacket and the random wavepacket in the product region
	eta_a=exp(-(p0-k)**2/(4d0*a0)-im*x0*k)*const
        eta_b=exp(-(p0-k)**2/(4d0*a0)+im*x0*k)*const
	sab = f_c(i)/(conjg(eta_b)*eta_a)
	write(101,1000) en,eta_a,eta_b,sab,abs(sab)**2
	enddo
	
!	nn(1) = Nt
!	call fourn(dc,nn,1,1)
!	do i=1,Nt
!         y = (i-1)/Nt/dt
!        IF(i .GT. Nt/2) y = -(N1+1-i)/Nt/dt
!	write(59,1000) y,dc(i)
!	ENDDO

c	write(57,1000) e0,csave,c0

cc777	continue

	do i=2,n1
	 if(abs(psi(i))**2.gt.1d-4) then 
	  u1=abs(psi(i-1))**2
	  u2=abs(psi(i))**2
	  u3=abs(psi(i+1))**2
	  der=(u3-u1)/2.d0/dx/u2
	  der2=(u3+u1-2.d0*u2)/dx/dx/u2
	  write(55,1000) x(i),abs(psi(i)),der2
	 endif
	enddo


1000    format(1000(e14.7,1x))

        stop
ccccccccccccccc do Chebyshev operator time propagation ccccccccccccc
c101        DO 3 J=1,NT
c         do jj=1,nout
c          call cheb(n1,dt,v,ak2,psi)
c         enddo
c         CALL AVEC(N1,psi0,psia,psib,psi,c0,c1,c2,c3,c4)
c         t=dt*j*nout
c         write(200,*) t,c0
c         write(201,*) t,c1
c         write(202,*) t,c2
c         write(203,*) t,c3
c         write(204,*) t,c4
c  3     CONTINUE
c        stop
ccccccccccccccc do Chebyshev operator order sequence ccccccccccccccc
c102          call seq(n1,nt,v,ak2,psi0,psia,psib,psi)
c        close(200)
c        close(201)
c        close(202)
c        close(203)
c        close(204)
c
c
c
C        STOP
        END program
!----------------------------------------------------
!        SUBROUTINE cmat(N,PSI,c12)
!        implicit none
!        integer*4,intent(IN) :: N
!        complex*16,intent(IN) :: psi(N)
!        complex*16,intent(OUT) :: c12
!        complex*16 :: im,psib(N)
!        integer*4 :: i
!        real*8 xmin,xmax,dx,xb,beta,pb,x
!        common/grid/xmin,xmax,dx
!        im = (0d0,1d0)
!        xb = -6d0
!        pb = 6d0
!        do i=1,N
!        x = xmin+(i-1)*dx
!        psib(i) = exp(-beta*(x-xb)**2+im*pb*(x-xb))
!        C12 = conjg(psib(i))*psi(i) + c12
!        enddo
!       return
!        END SUBROUTINE

C***************************************************************
      SUBROUTINE AVE(N,PSI,RPSI,ANS)
C      CALCULATION OF THE COMPLEX INNER PRODUCT
C      PSI AND RPSI ARE INPUT,
C      INTEGRAL OF CONJG(PSI)*RPSI*DX RETURNED
C
	implicit real*8 (a-h,o-z)
      COMPLEX*16 ans,PSI(1),RPSI(1)

      ANS=(0.d0,0.d0)
      DO 1 K=1,N
      ANS=ANS+CONJG(PSI(K))*RPSI(K)
 1    CONTINUE
      RETURN
      END
C***************************************************************
C***************************************************************
      	subroutine avec(N,dx,x,psi0,psib,psi,a1,a2,a0,as,cm1,cm2)
	implicit real*8 (a-h,o-z)
	parameter(nj=2048)
	real*8 x(nj)
      	complex*16 psi0(nj),psi(nj),psib(nj),a1,a2
      A0=0.d0
	as=0.d0
      cm1=0.d0
      cm2=0.d0
	a1=(0d0,0d0)
	a2=(0d0,0d0)
      DO 1 K=1,N
      if(x(k).gt.0d0) A0=A0+abs(psi(k))**2
        As=As+abs(psi(k))**2*x(k)
cc        as=as+abs(psi(k))**2*x(k)**2
	a1=a1+conjg(psi0(k))*psi(k)
	a2=a2+conjg(psib(k))*psi(k)
	cm1=cm1+abs(psi(k))**2
	cm2=cm2+abs(psib(k))**2*abs(psi(k))**2
 1    CONTINUE
	a0=a0*dx
	as=as*dx
	a1=a1*dx
	a2=a2*dx
	cm1=cm1*dx
	cm2=cm2*dx
	if(abs(a2).lt.1d-16) a2=(0d0,0d0)
	if(abs(a1).lt.1d-16) a1=(0d0,0d0)
	if(abs(cm1).lt.1d-16) cm1=0d0
	if(abs(cm2).lt.1d-16) cm2=0d0
      RETURN
      END
c************************************************************
        subroutine split(ni,dt,dt2,v,ak2,cwf)
        implicit real*8(a-h,o-z)
        parameter (nj=2048)
        dimension ak2(ni)
        dimension nn(1)
        complex*16 cwf(nj),aux(nj),psi0(nj),C,df(nj),v(ni),
     &             psia(nj),psib(nj),fr(nj),im,z
	real*8 prob,x
	common/grid/xmin,xmax,dx
        im=(0.d0,1.d0)
	nn(1)=ni
        call fourn(cwf,nn,1,1)
         call diff(cwf,ak2,dt2,ni)
         call fourn(cwf,nn,1,-1)
          do 12 i=1,ni
12       cwf(i)=cwf(i)/ni

         call phase(cwf,v,dt2,ni)

         call fourn(cwf,nn,1,1)
         call diff(cwf,ak2,dt2,ni)
         call fourn(cwf,nn,1,-1)
          do 13 i=1,ni
13       cwf(i)=cwf(i)/ni

        call phase(cwf,v,dt2,ni)
c
        return
        end

        subroutine diff(cwf,ak2,ts,ni)
        implicit real*8(a-h,o-z)
        complex*16 nim,cwf(ni)
        dimension ak2(ni)
        nim=(0.d0,-1.d0)
         do 11 i=1,ni
11      cwf(i)=cwf(i)*exp(nim*ts*ak2(i))
        return
        end

        subroutine phase(cwf,v,ts,ni)
        implicit real*8(a-h,o-z)
        complex*16 nim,cwf(ni),v(ni)
        nim=(0.d0,-1.d0)
         do 11 i=1,ni
        if(ts.ge.0) cwf(i)=cwf(i)*exp(nim*ts*v(i))
        if(ts.lt.0) cwf(i)=cwf(i)*exp(nim*ts*conjg(v(i)))
11      continue


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
c
c*************************************************************
        subroutine psi_c(N,psi)
        implicit real*8 (a-h,o-z)
        parameter(nj=2048)
        complex*16 psi(nj)
ccccccccccccccccccc set up the wavepackets ccccccccccc
        sum=0.d0
        iseed=31456729
        do 16 i=1,N
         psi(i)=(2.*ran0(iseed)-1.)
         sum=sum+abs(psi(i))**2
16      continue
        an=sqrt(abs(sum))
        sum=0.d0
        do 17 i=1,N
          psi(i)=psi(i)/an
          sum=sum+abs(psi(i))**2
17      continue
        write(*,*) 'wavepacket normalizaiton',sum
        write(12,*) 'wavepacket normalizaiton',sum
        return
        end
C--------------------------------------------------
      real*8 function ran0(idum)
      integer idum,IA,IM,IQ,IR,MASK
      real*8 AM
      parameter(IA=16807,IM=2147483647,AM=1./IM,
     &          IQ=127773,IR=2836,MASK=123459876)
      integer k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if(idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      end
C--------------------------------------------------

        subroutine row_mat(n,a,v,eig)
        implicit real*8 (a-h,o-z)
        parameter(nj=2048)
        real*8 ans(nj),a(nj,nj),v(nj),eig(nj)
cccccccccccccccccccc row times a matrix cccccccccccccccccccccccccccccc
        do 10 i=1,n
          ans(i)=(0.d0,0.d0)
          do 11 j=1,n
11          ans(i)=ans(i)+v(j)*a(j,i)
10      eig(i)=ans(i)/v(i)
         return
        end

        subroutine mat_mat(n,a,b)
        implicit real*8 (a-h,o-z)
        parameter(nj=2048)
        real*8 ans,an(nj),a(nj,nj),b(nj,nj)
cccccccccccccccccccc matrix times column ccccccccccccccccccccccccccccc
        do 10 i=1,n
          do 10 j=1,n
            ans=(0.d0,0.d0)
            do 12 k=1,n
12            ans=ans+a(k,i)*b(k,j)
          if(i.eq.j) an(i)=sqrt(ans)
10      continue
        do 11 j=1,n
          do 11 i=1,n
            b(i,j)=b(i,j)/an(j)
11      a(i,j)= a(i,j)/an(j)
        return
        end
c
        subroutine mat_col(n,a,v,ans)
        implicit real*8 (a-h,o-z)
        parameter(nj=2048)
        real*8 a(nj,nj)
	complex*16 v(nj),ans(nj)

        do 10 i=1,n
          ans(i)=0.d0
          do 10 j=1,n
10      ans(i)=ans(i)+a(i,j)*v(j)

         return
        end


	subroutine cheb(N,T,V,AK2,PSI)
	implicit real*8 (a-h,o-z)
	parameter(nj=2048)
	real*8 ak2(nj),vt(nj),gam(nj),coef(4000)
	complex*16 rt,phase,v(nj),psi(nj),w1(nj),w2(nj),w3(nj)
	SAVE IFLAG,COEF,gam,VT,RT,phase
	data iflag/0/
      	IF (IFLAG.EQ.0) THEN
      	VMAX=dreal(V(1))
      	VMIN=dreal(V(1))
      	AKMAX=0.
      	DO 10 I=1,N
	 gam(i)=exp(dimag(v(i)))
      	 VMIN=DMIN1(dreal(V(I)),VMIN)
      	 VMAX=DMAX1(dreal(V(I)),VMAX)
      	 AKMAX=DMAX1(AK2(I),AKMAX)
10  	CONTINUE
cmy	vmax=vmax*2.d0
cmy	vmin=vmin/2.d0
	DELTAE=(VMAX-VMIN)+AKMAX
        DET=DELTAE*T/2.
        WRITE(12,2001)DET,DELTAE,VMAX,VMIN,AKMAX
2001    FORMAT(//,3X,5E16.5)
C
        IF (DET.LT.100.) FAC=1.5
        IF (DET.GT.100. .AND. DET.LT.400.) FAC=1.2
        IF (DET.GT.400. ) FAC=1.05
C
C      COMPUTE EXPANSION COEFFICIENTS
C
        ACC=1.E-16
  40    NTEST=DET*FAC
        WRITE(16,2002)DET,NTEST
 2002   FORMAT(2X,E16.5,2X,I5)
C
        CALL COF(DET,NTEST,COEF)
        J=NTEST
 20     CONTINUE
        IF (ABS(COEF(J)).LT.ACC .AND. J.GT. 1) THEN
         J=J-1
         GOTO 20
        END IF
        IF (J.EQ.NTEST) THEN
         FAC=FAC*1.1
         GOTO 40
        END IF
        NCHEB=J
        WRITE(*,2000)NCHEB,DET,COEF(NCHEB)
        WRITE(12,2000)NCHEB,DET,COEF(NCHEB)
 2000 FORMAT(/,3X,'NCHEB,DET,COEF(NCHEB) ',I4,2X,2E15.6,/)

C
C       SHIFT THE POTENTIAL.
C       RT*HAM = -I*HAMNORM
C
      RT = (0.,-2.)/DELTAE
      DO 50 I=1,N
      VT(I)=dreal(V(I))-VMIN-DELTAE/2.
  50  CONTINUE
C
      PHASE=EXP((0.,-1.)*(DELTAE/2.+VMIN)*T)
      IFLAG=1
      END IF
C     END OF INITIALIZATION
C
C       START OF RECURSION
C
      DO 60 I=1,N
      W1(I)=PSI(I)
      PSI(I)=COEF(1)*PSI(I)
 60   CONTINUE
      CALL HAM(N,AK2,VT,W1,W2)
      DO 70 I=1,N
      W2(I)=gam(i)*W2(I)*RT
      PSI(I)=PSI(I)+W2(I)*COEF(2)
 70   CONTINUE
C
C	CHEBYHEV RECURSION

      DO 80 J=3,NCHEB
      CALL HAM(N,AK2,VT,W2,W3)
      DO 75 I=1,N
      W3(I)=gam(i)*(W3(I)*RT*2.+gam(i)*W1(I))
      W1(I)=W2(I)
      W2(I)=W3(I)
      PSI(I)=PSI(I)+COEF(J)*W3(I)
 75   CONTINUE
 80   CONTINUE
      DO 90 I=1,N
 90   PSI(I)=PHASE*PSI(I)
C
C      PHASE SHIFT INCLUDED IN PSI
C
      RETURN
      END

        subroutine seq(N,nt,V,AK2,psi0,psia,psib,psi)
        implicit real*8 (a-h,o-z)
        parameter(nj=2048)
        real*8 ak2(nj),vt(nj),gam(nj),coef(4000)
        complex*16 v(nj),psi(nj),w1(nj),w2(nj),w3(nj)
     &          ,psi0(nj),psia(nj),psib(nj),c0,c1,c2,c3
        SAVE IFLAG,COEF,gam,VT,RT,phase
        data iflag/0/
        IF (IFLAG.EQ.0) THEN
        VMAX=dreal(V(1))
        VMIN=dreal(V(1))
        AKMAX=0.
        DO 10 I=1,N
	 gam(i)=exp(dimag(v(i)))
         VMIN=DMIN1(dreal(V(I)),VMIN)
         VMAX=DMAX1(dreal(V(I)),VMAX)
         AKMAX=DMAX1(AK2(I),AKMAX)
10      CONTINUE
cmy     vmax=vmax*2.d0
cmy     vmin=vmin/2.d0
        DELTAE=(VMAX-VMIN)+AKMAX
C
C       SHIFT THE POTENTIAL.
C       RT*HAM = HAMNORM
C
      RT = 2./DELTAE
      DO 50 I=1,N
      VT(I)=(V(I))-VMIN-DELTAE/2.
  50  CONTINUE
      IFLAG=1
      END IF
C     END OF INITIALIZATION
C
C       START OF RECURSION
C
      DO 60 I=1,N
      W1(I)=PSI(I)
 60   CONTINUE
      CALL HAM(N,AK2,VT,W1,W2)
      DO 70 I=1,N
      W2(I)=gam(i)*W2(I)*RT
 70   CONTINUE
c       CALL AVEC(N,psi0,psia,psib,w2,c0,c1,c2,c3,c4)
c       t=2*1.d0
c       write(200,*) t,c0
c       write(201,*) t,c1
c       write(202,*) t,c2
c       write(203,*) t,c3
c       write(204,*) t,c4

C
C       CHEBYHEV RECURSION

      DO 80 J=3,nt
      CALL HAM(N,AK2,VT,W2,W3)
      DO 75 I=1,N
       W3(I)=gam(i)*(W3(I)*RT*2.-gam(i)*W1(I))
       W1(I)=W2(I)
       W2(I)=W3(I)
 75   CONTINUE
c       CALL AVEC(N,psi0,psia,psib,w3,c0,c1,c2,c3,c4)
       t=j*1.d0
       write(200,*) t,c0
       write(201,*) t,c1
       write(202,*) t,c2
       write(203,*) t,c3
       write(204,*) t,c4

 80   CONTINUE
C
      RETURN
      END



	subroutine ham(ni,ak2,v,a,b)
	implicit real*8(a-h,o-z)
	parameter(nj=2048)
	real*8 ak2(nj),v(nj)
	complex*16 cwf(nj),a(nj),b(nj)
	dimension nn(1)
	nn(1)=ni

	do 11 i=1,ni
11	cwf(i)=a(i)
        call fourn(cwf,nn,1,1)
	do 12 i=1,ni
12	cwf(i)=cwf(i)*ak2(i)
        call fourn(cwf,nn,1,-1)
         do 13 i=1,ni
13      b(i)=cwf(i)/ni+v(i)*a(i)

	return
	end

C************************************************************
C
        SUBROUTINE COF(R,NCH,CF)
C       CALCULATES BESSEL COFICIENTS
C
C          WARNING: REAL STATEMENTS ARE MACHINE-DEPENDENT
C                   ALSO CHECK MMBSJN DATA DECLARATIONS
C
        implicit real*8 (a-h,o-z)
        REAL*8 B(3000),R1
        REAL*8 CF(2000)
        NR=0
        R1=R
C       R1=1.0E2
C       NCH=101
        CALL MMBSJ1(R1,NCH,B,IER)
C       WRITE(6,2000)B(1),B(NCH)
C2000   FORMAT(2X,'IN COF',2X,2E16.5,/)
C       IF(R1.GT.0)STOP 13
C
        DO 1 I=2,NCH
        CF(I)=2.*B(I)
1       CONTINUE
        CF(1)=B(1)
        RETURN
        END
c
c************************************************************
      SUBROUTINE MMBSJ1 (ARG,N,B,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IER
      REAL*8        ARG,B(3000)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L,LARGEX,MAGX,NCALC,NN,NBMX,M,NSTART,NEND
      REAL*8      TEST,TEMPA,TEMPB,TEMPC,EXPARG,P
      REAL*8                 RSIGN,SUM,TOVER,PLAST,POLD,PSAVE,PSAVEL
      REAL*8                 DSIG,RTEN,TMPA4,SMALLX,XINF
      DATA               DSIG/14.4d0/
C       DSIG CHANGED FROM 14.4 TO 7.4
      DATA               RTEN/38.d0/
C       RTEN CHANGED FROM 322E0 TO 38E0
      DATA               EXPARG/741.66748319914d0/
      DATA               XINF/.12650140831d+36/
C    CHANGE IN XNIF   FROM E 325 TO E 36
       DATA             LARGEX/500000/
C                                  FIRST EXECUTABLE STATEMENT

      IER = 0
      TEMPA =dABS(ARG)
      MAGX = TEMPA
ccc        print*,'arg n',arg,n,magx,largex

      IF(N.GT.0 .AND. MAGX.LE.LARGEX) GO TO 10
C                                  ERROR RETURN -- ARG,N IS OUT OF RANGE
      IER = 129
      B(1) = XINF
      IF(N.LT.2) GO TO 9000
      DO 5 L=2,N
         B(L) = XINF
    5 CONTINUE
      GO TO 9000
   10 RSIGN = 1.
      NCALC = N
C                                  USE 2-TERM ASCENDING SERIES FOR
C                                    SMALL ARG
      TMPA4 = TEMPA**4.
      SMALLX = .1**DSIG
      IF(TMPA4.GE.SMALLX) GO TO 20
C                                  TWO-TERM ASCENDING SERIES FOR
C                                    SMALL ARG
      TEMPA = 1.d0
      TEMPB = -.25*ARG*ARG*RSIGN
      B(1) = 1.d0+TEMPB
      IF(N.EQ.1) GO TO 9005
      DO 15 NN=2,N
         TEMPA = TEMPA*ARG/(dFLOAT(2*NN-2))
         B(NN) = TEMPA*(1.+TEMPB/(dFLOAT(NN)))
   15 CONTINUE
      GO TO 9005
C                                  INITIALIZE THE CALCULATION OF P*S
   20 NBMX = N-MAGX
      NN = MAGX+1
      PLAST = 1.d0
      P = (dFLOAT(2*NN))/TEMPA
C                                  CALCULATE GENERAL SIGNIFICANCE TEST
      TEST = 2.d0*1.d1**DSIG
      M = 0
      IF(NBMX.LT.3) GO TO 30
C                                  CALCULATE P*S UNTIL NN=N-1.
C                                    CHECK FOR POSSIBLE OVERFLOW.
      TOVER = 1.d1**(RTEN-DSIG)
      NSTART = MAGX+2
      NEND = N-1
      DO 25 NN=NSTART,NEND
         POLD = PLAST
         PLAST = P
         P = (dFLOAT(2*NN))*PLAST/TEMPA-RSIGN*POLD
         IF(P-TOVER) 25, 25, 35
   25 CONTINUE
      NN = NEND
C                                  CALCULATE SPECIAL SIGNIFICANCE TEST
C                                    FOR NBMX.GT.2.
C
      TEST =dMAX1(TEST,dSQRT(PLAST*1.d1**DSIG)*dSQRT(2.d0*P))
C
C                                  CALCULATE P*S UNTIL SIGNIFICANCE
C                                    TEST PASSES
   30 NN = NN+1
      POLD = PLAST
      PLAST = P
      P = (dFLOAT(2*NN))*PLAST/TEMPA-RSIGN*POLD
      IF(P.LT.TEST) GO TO 30
      IF(M.EQ.1) GO TO 55
C                                  FOR J*S, A STRONG VARIANT OF THE TEST
C                                    IS NECESSARY. CALCULATE IT, AND
C                                    CALCULATE P*S UNTIL THIS TEST IS
C                                    PASSED.
      M = 1
      TEMPB = P/PLAST
      TEMPC = (dFLOAT(NN+1))/TEMPA
      IF(TEMPB+1./TEMPB.GT.2.*TEMPC)TEMPB=TEMPC+dSQRT(TEMPC**2-1.)
      TEST = TEST/dSQRT(TEMPB-1./TEMPB)
      IF(P-TEST) 30, 55, 55
C                                  TO AVOID OVERFLOW, DIVIDE P*S BY
C                                    TOVER.  CALCULATE P*S UNTIL
C                                    ABS(P).GT.1.
   35 TOVER = 1.d1**RTEN
      P = P/TOVER
      PLAST = PLAST/TOVER
      PSAVE = P
      PSAVEL = PLAST
      NSTART = NN+1
   40 NN = NN+1
      POLD = PLAST
      PLAST = P
      P = (dFLOAT(2*NN))*PLAST/TEMPA-RSIGN*POLD
      IF(P.LE.1.) GO TO 40
      TEMPB = (dFLOAT(2*NN))/TEMPA
      TEMPC = .5d0*TEMPB
      TEMPB = PLAST/POLD

      IF(TEMPB+1.d0/TEMPB.GT.2.d0*TEMPC)TEMPB=TEMPC+dSQRT(TEMPC**2-1.d0)
C
C                                  CALCULATE BACKWARD TEST, AND FIND
C                                    NCALC, THE HIGHEST NN SUCH THAT THE
C                                    TEST IS PASSED.
      TEST = .5d0*POLD*PLAST*(1.d0-1.d0/TEMPB**2)/1.d1**DSIG
      P = PLAST*TOVER
      NN = NN-1
      NEND = MIN0(N,NN)
      DO 45 NCALC=NSTART,NEND
         POLD = PSAVEL
         PSAVEL = PSAVE
         PSAVE = (dFLOAT(2*NN))*PSAVEL/TEMPA-RSIGN*POLD
         IF(PSAVE*PSAVEL-TEST) 45, 45, 50
   45 CONTINUE
      NCALC = NEND+1
   50 NCALC = NCALC-1
C                                  THE SUM B(1)+2B(3)+2B(5)... IS USED
C                                    TO NORMALIZE. M, THE COEFFICIENT OF
C                                    B(NN), IS INITIALIZED TO 2 OR 0.
   55 NN = NN+1
      M = 2*NN-4*(NN/2)
C                                  INITIALIZE THE BACKWARD RECURSION AND
C                                    THE NORMALIZATION SUM
      TEMPB = 0.
      TEMPA = 1.d0/P
      SUM = (dFLOAT(M))*TEMPA
      NEND = NN-N
      IF(NEND) 80, 70, 60
C                                  RECUR BACKWARD VIA DIFFERENCE
C                                    EQUATION, CALCULATING (BUT NOT
C                                    STORING) B(NN), UNTIL NN=N.
   60 DO 65 L=1,NEND
         NN = NN-1
         TEMPC = TEMPB
         TEMPB = TEMPA
         TEMPA = ((dFLOAT(2*NN))*TEMPB)/ARG-RSIGN*TEMPC
         M = 2-M
         SUM = SUM+(dFLOAT(M))*TEMPA
   65 CONTINUE
C                                  STORE B(NN)
   70 B(NN) = TEMPA
      IF(N.GT.1) GO TO 75
C                                  N=1.  SINCE 2*TEMPA IS ADDED TO THE
C                                    SUM, TEMPA MUST BE SUBTRACTED
      SUM = SUM-TEMPA
      GO TO 110
C                                  CALCULATE AND STORE B(NN-1)
   75 NN = NN-1
      B(NN) = ((dFLOAT(2*NN))*TEMPA)/ARG-RSIGN*TEMPB
      IF(NN.EQ.1) GO TO 105
      M = 2-M
      SUM = SUM+(dFLOAT(M))*B(NN)
      GO TO 90
C                                  NN.LT.N, SO STORE B(NN) AND SET
C                                  HIGHER ORDERS TO ZERO
   80 B(NN) = TEMPA
      NEND = -NEND
      DO 85 L=1,NEND
         ITEMP = NN+L
         B(ITEMP) = 0.0
   85 CONTINUE
   90 NEND = NN-2
      IF(NEND.EQ.0) GO TO 100
C                                  CALCULATE VIA DIFFERENCE EQUATION AND
C                                    STORE B(NN), UNTIL NN=2
      DO 95 L=1,NEND
         NN = NN-1
         B(NN) = ((dFLOAT(2*NN))*B(NN+1))/ARG-RSIGN*B(NN+2)
         M = 2-M
         SUM = SUM+(dFLOAT(M))*B(NN)
   95 CONTINUE
C                                  CALCULATE B(1)
  100 B(1) = 2.*B(2)/ARG-RSIGN*B(3)
  105 SUM = SUM+B(1)
C                                  NORMALIZE--IF IZE=1, DIVIDE SUM BY
C                                    COSH(ARG). DIVIDE ALL B(NN) BY SUM.
  110 CONTINUE
      DO 115 NN=1,N
  115 B(NN) = B(NN)/SUM
      IF(NCALC.EQ.N) GO TO 9005
      IER = 129+NCALC
 9000 CONTINUE
C      CALL UERTS1(IER,'MMBSJN')
      CALL UERTS1
 9005 RETURN
      END
      SUBROUTINE UERTS1
      PRINT *,'ERROR IN IMSL MMBSJN'
      STOP
      END
	



