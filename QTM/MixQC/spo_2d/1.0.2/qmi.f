	program qm
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc Coupled Morse oscillators ccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	implicit real*8 (a-h,o-z)
	integer*4 ni,nj
	parameter (ni=256, nj=256)
	real*8 ak2(ni,nj)
	complex*16 v(ni,nj)
	complex*16 im,psi(ni,nj)
	complex*16 psi0(ni,nj),c1,c2
	common /params/ d,x0,z0,di,wi
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
C
C
	im=(0.d0,1.d0)
	pi=4.d0*atan(1.d0)
	N2=ni*nj
	d=160.d0
	X0=1.40083d0
c	z0=1.04435d0
	z0=1.0435d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c	open(15,file='cor1')
c	open(16,file='cor2')
	open(19,file='eni_qm')
	open(17,file='wf0')
	open(18,file='wft')
	open(5,file='INq',status='old')
c time propagation parameters
	read(5,*) dt,kmax,nout
c grid spacing
	read(5,*)dx,dy 
c initial wave packet parameters
	read(5,*) alfa,q0,p0,beta
c potential surface parameters
        read(5,*)DI,WI
	close(5)
c inverse mass prefactor 1/2/m
c hydrogen
	am=1.0
c deuterium m=2*mh
c	am=0.5d0
c	alfa=alfa*sqrt(2.d0)
cc tritium=3*mh
c	am=1.d0/3.d0
cc	alfa=alfa*sqrt(3.d0)
cc
	xmin=0.1d0
	xmax=xmin+ni*dx
	ymin=0.1d0
	ymax=ymin+nj*dy
	ds=dx*dy
        q0=q0+x0
        write(*,*)'initial state',alfa,q0,p0
        write(*,*) 'kinetic energy coupling',beta
c define vector ak2
	consI=2.d0*pi/dx/ni
	consJ=2.d0*pi/dy/nj
	ni2=ni/2
	nj2=nj/2
	do 10 j=1,nj
	 akj=(j-1)*consj
	 if (j.gt.nj2) akj=-(nj+1-j)*consj
	 do 20 i=1,ni
	  aki=(i-1)*consi
	  if (i.gt.ni2) aki=-(ni+1-i)*consi
20	ak2(i,j)=(aki*aki-2d0*beta*aki*akj+akj*akj)/2.d0/am
10	continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c------------ set up the initial wavefunction -----------------------
        an0=dsqrt(2d0*alfa/pi)
        ch_n=0d0
	do 222 j=1,nj
	 y=ymin+dy*(j-1)
	 do 222 i=1,ni
	  x=xmin+dx*(i-1)
          psi(i,j)=an0*exp(-alfa*((x-q0)**2+(y-q0)**2)
     &          +im*p0*(x-q0+y-q0))
          psi0(i,j)=psi(i,j)
	   if(abs(psi(i,j)).gt.1d-3) write(17,*) x,y,abs(psi(i,j))
         ch_n=ch_n+abs(psi(i,j))**2*ds
222     continue
        write(*,*) 'NORM=',ch_n
cccccccccccccccccc propagate psi ccccccccccccccccccccccccccccccccc
	call ham(ni,nj,N2,ak2,psi,v,en0)
	call correl(ni,nj,psi0,psi,c1,c2)
	t=0.d0
        write(19,*) t, en0
c	write(15,1000) t,c1,abs(c1)
c	write(16,1000) 2d0*t,c2,abs(c2)
	call split(kmax,nout,ni,nj,N2,dt,t,v,ak2,psi0,psi)
	call ham(ni,nj,N2,ak2,psi,v,en0)
        write(*,*) 'FINAL ENERGY',en0
        write(19,*) t, en0
c------------ print the final  wavefunction -----------------------
	do 23 j=1,nj
         y=ymin+(j-1)*dy
	 do 23 i=1,ni
	  x=xmin+dx*(i-1)
	  if(abs(psi(i,j)).gt.1d-3) write(18,*) x,y,abs(psi(i,j))
23	continue
c--------------------------------------------------------------------
1000    format(20(e14.7,1x))
	stop
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine ham(ni,nj,N2,ak2,psi,v,en0)
	implicit real*8(a-h,o-z)
	parameter(nij=256)
	complex*16 psi(nij,nij),tpsi(nij,nij),vpsi(nij,nij)
	dimension ak2(nij,nij)
	complex*16 v(nij,nij),c,anc
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
	call aver(ni,nj,psi,psi,anC)
	call potpsi(ni,nj,v,psi,vpsi)
	call aver(ni,nj,psi,vpsi,C)
c	write(*,*) 'potential energy',C/anc
	en0=dreal(c)
	call kinpsi(ni,nj,N2,ak2,psi,tpsi)
	call aver(ni,nj,psi,tpsi,C)
c	write(*,*) 'kinetic energy', C/anc
	en0=(en0+dreal(c))/abs(anc)
c	write(*,*) 'TOTAL E', en0
	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine potpsi(ni,nj,v,psi,vpsi)
	implicit real*8(a-h,o-z)
	parameter(nij=256)
	complex*16 psi(nij,nij),vpsi(nij,nij),v(nij,nij),im
	common /params/ d,x0,z0,di,wi
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
c------------- set up potetial -------------------------------------------
        im=(0d0,1d0)
        do 10 j=1,nj
         y=ymin+(j-1)*dy
         ry=exp(-z0*(y-x0))
         do 10 i=1,ni
          x=xmin+(i-1)*dx
          rx=exp(-z0*(x-x0))
          vr=d*(rx-1d0)**2+d*(ry-1d0)**2
          if(vr.gt.5d0*d) vr=5d0*d
          vi=0d0
          r=dsqrt(rx*rx+ry*ry)
          if(r.gt.wi) vi=vi+di*(r-wi)**2
          v(i,j)=vr-im*vi
          vpsi(i,j)=psi(i,j)*v(i,j)
10      continue
        return
        end
        

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine kinpsi(ni,nj,N2,ak2,cx1,aux)
        implicit real*8(a-h,o-z)
	parameter(nij=256)
        dimension ak2(nij,nij)
        complex*16 cx1(nij,nij),aux(nij,nij)
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
        subroutine split(kmax,nout,ni,nj,N2,dt,t,v,ak2,psi0,cwf)
        implicit real*8(a-h,o-z)
	parameter(ij=256)
        dimension ak2(ij,ij)
        dimension nn(2)
        complex*16 cwf(ij,ij),v(ij,ij)
        complex*16 psi0(ij,ij),hpsi(ij,ij),c1,c2
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        dt2=dt*0.5
	n2=ni*nj
        nn(1)=ni
        nn(2)=nj
        do 11 k=1,kmax
	 t=t+dt
         call fourn(cwf,nn,2,1)
         call diff(cwf,ak2,dt2,ni,nj)
         call fourn(cwf,nn,2,-1)
         do 12 j=1,nj
          do 12 i=1,ni
12       cwf(i,j)=cwf(i,j)/N2
 
         call phase(cwf,v,dt2,ni,nj)
 
         call fourn(cwf,nn,2,1)
         call diff(cwf,ak2,dt2,ni,nj)
         call fourn(cwf,nn,2,-1)
         do 13 j=1,nj
          do 13 i=1,ni
13       cwf(i,j)=cwf(i,j)/N2

        call phase(cwf,v,dt2,ni,nj)
c	call correl(ni,nj,psi0,cwf,c1,c2)
	call ham(ni,nj,N2,ak2,cwf,v,en0)
c	write(15,1000) t,c1,abs(c1)
c	write(16,1000) 2d0*t,c2,abs(c2)
        write(19,*) t,en0
11	continue

1000    format(20(e14.7,1x))
15	return
	end

c-------------------------------------------------------
c-------------------------------------------------------
	subroutine diff(cwf,ak2,ts,ni,nj)
	implicit real*8(a-h,o-z)
	parameter(nij=256) 	
	complex*16 nim,cwf(nij,nij)
	dimension ak2(nij,nij)
	nim=(0.d0,-1.d0)
	do 11 j=1,nj
	 do 11 i=1,ni
c11	cwf(i,j)=cwf(i,j)*exp(nim*ts*ak2(i,j))
11	cwf(i,j)=cwf(i,j)*exp(-ts*ak2(i,j))
	return
	end

	subroutine phase(cwf,v,ts,ni,nj)
	implicit real*8(a-h,o-z)
	complex*16 nim,cwf(ni,nj)
	complex*16 v(ni,nj)
	nim=(0.d0,-1.d0)
	do 11 j=1,nj
	 do 11 i=1,ni
c	if (ts.ge.0) cwf(i,j)=cwf(i,j)*exp(nim*ts*v(i,j))
c11	if (ts.lt.0) cwf(i,j)=cwf(i,j)*exp(nim*ts*conjg(v(i,j)))
11	cwf(i,j)=cwf(i,j)*exp(-ts*v(i,j))
	return
	end

c-------------------------------------------------------
c-------------------------------------------------------
	subroutine correl(ni,nj,psi0,psi,c1,c2)
	implicit real*8(a-h,o-z)
	parameter(ij=256)
	complex*16 psi0(ij,ij),psi(ij,ij),c1,c2
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
ccc---------- compute correlation functions ------------------------
	c1=(0d0,0d0)
	c2=(0d0,0d0)
	do 10 j=1,nj
	 y=ymin+(j-1)*dy
	 do 10 i=1,ni
	  x=xmin+(i-1)*dx
          c1=c1+conjg(psi0(i,j))*psi(i,j)
          c2=c2+psi(i,j)**2
10	continue
        c1=c1*dx*dy
        c2=c2*dx*dy
	return
	end

c-------------------------------------------------------
c-------------------------------------------------------
        subroutine aver(ni,nj,psi0,psi,c)
        implicit real*8(a-h,o-z)
        parameter(ij=256)
        complex*16 c,psi0(ij,ij),psi(ij,ij)
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
	c=(0.d0,0.d0)
	ct=0.d0
        do 11 j=1,nj
         do 11 i=1,ni
	ct=ct+abs(psi0(i,j))**2
11      c=c+psi(i,j)*conjg(psi0(i,j))
        c=c/dx/dy
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

