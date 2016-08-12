	program qm
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc 3D H+H2 Reaction in Jacobi coordinates  cccccccccccccccc
cccccccc use Gaussians instead of the eigenstates ccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	implicit real*8 (a-h,o-z)
	integer*4 ni,nj
	parameter (ni=256, nj=256,nth=175)
	real*8 ak2(ni,nj),rp(2,2),v(ni,nj,nth)
	real*8 V_a,V_b,ka,kb,ma,mb,mc
	complex*16 im,c,psi(ni,nj,nth),psib(ni,nj,nth)
	complex*16 eta_a,eta_b,sab
        complex*16,allocatable :: corr(:),cf(:)
	common /wave/ alfa,q0,y0,px0,py0,dt
	common /params/ di,wi
	common /grid/ xmin,ymin,xmax,ymax,dx,dy,amuh,amud,am(2)
	common/bend/z(nth),w(nth),rot(ni,nj),tr(nth,nth)
	common/mass/ma,mb,mc
	common/iw/iwc,iws
!	common/four/Emin,Emax,de
c---------- constants 
	im=(0.d0,1.d0)
	pi=4.d0*atan(1.d0)
	n2=ni*nj
	iwc=26
	iws=27
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	open(15,file='prob')
	open(16,file='nofe')
	open(17,file='corq')
        open(18,file='cor')
!        open(19,file='corf')
	open(20,file='wf0')
	open(21,file='wft')
!	open(22,file='smat')
        open(26,file='corf')
        open(27,file='smat')
	open(5,file='INq',status='old')
c time propagation parameters
	read(5,*) dt,kmax
	read(5,*) ma,mb,mc
c grid spacing
	read(5,*) xmin,ymin
	read(5,*) xmax,ymax
	read(5,*) cut 
c initial wave packet parameters
	read(5,*) alfa,beta,x0,y0,px0,py0 
c number of wavewpackets 
	read(5,*) nmom,dmom
	read(5,*) nk,iso
	close(5)

        allocate(cf(kmax),corr(kmax))
c-------- masses -------------------------------------------
        ame=1836.15d0
        amh=ame*1.007825d0
        amd=ame*1.007825d0*iso
cc        amo=amh*15.994915d0
	amo=amh
        am(1)=amo*(amh+amd)/(amo+amh+amd)
        am(2)=amh*amd/(amh+amd)
        amud=amd/(amh+amd)
        amuh=amh/(amh+amd)
	ma=1836.15d0*ma
	mb=1836.15d0*mb
	mc=1836.15d0*mc
        write(*,*) 'REDUCED MASSES', am(1),am(2)
!	xmin=0.1d0
!	xmax=xmin+ni*dx
!	ymin=0.1d0
!	ymax=ymin+nj*dy
	dx=(xmax-xmin)/ni
	dy=(ymax-ymin)/nj
	ds=dx*dy
        print *,'GRID SIZE'
        print *,'xmin=',xmin
        print *,'xmax=',xmax
        print *,'ymin=',ymin
        print *,'ymax=',ymax
	print *,'dx=',dx
	print *,'dy=',dy
        print *,'nx=',NI
        PRINT *,'ny=',NJ
        print *,'nk=',NK      

c-------- product wavepacket -----------------------
	rp(1,1)=-amo/(amo+amh)
	rp(1,2)=-amh/(amo+amh)*(amo+amh+amd)/(amh+amd)
	rp(2,1)=1d0
	rp(2,2)=-amd/(amh+amd)
	yp=1.85d0
	xp=-x0
	pyp=py0
	pxp=-px0
!	xf=xp*rp(2,2)-yp*rp(1,2)
!	yf=-xp*rp(2,1)+yp*rp(1,1)
	xf = 6.d0
	yf = 1.444d0
	pxf=pxp*rp(2,2)-pyp*rp(1,2)
	pyf=-pxp*rp(2,1)+pyp*rp(1,1)
        print *,'INITIAL WAVEPACKET'
        print *,'alpha=',alfa,'x0=',x0
        print *,'beta=',beta,'y0=',y0
        print *,'px0=',px0,'py0=',py0
	write(*,*) 'FINAL WAVEPACKET'
        print *,'xf=',xf
        print *,'yf=',yf
	print *,'pxp=',pxp
	print *,'pyp=',pyp
        print *,'pxf=',pxf
        print *,'pyf=',pyf
c----------------------------------------------------

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
20	ak2(i,j)=(aki*aki/am(1)+akj*akj/am(2))/2.d0
10	continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	DO 888 nn=1,nmom
	t=0d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	call psigaus(ni,nj,nk,psi,v,pi,alfa,beta,x0,y0,px0,py0)
	call psiprod(ni,nj,nk,psib,pi,alfa,beta,xf,yf,pxp,pyp,rp)
	call plot(20,ni,nj,nk,psi,cut)
cccccccccccccccccc propagate psi ccccccccccccccccccccccccccccccccc
	 call correl(ni,nj,nk,psib,psi,poh,pod,phd,c)
	 write(15,1000) t,poh,pod,phd
	 write(17,1000) t,c,abs(c)
	 call split(kmax,ni,nj,nk,N2,ak2,v,psib,psi,dt,corr)
220	continue
	iwr=20+nn
	call plot(iwr,ni,nj,nk,psi,cut)
	call correl(ni,nj,nk,psib,psi,poh,pod,phd,c)
	call tw(corr,alfa,beta,px0,pxp,x0,xf,kmax,dt,sab)
!	do i=1,kmax
!	En = Emin+(i-1)*de
!	write(22,1000) En,abs(sab(i))**2
!	enddo
	write(16,1000) -px0,poh,pod,phd
	px0=px0+dmom
888 	CONTINUE
1000    format(20(e14.7,1x))
	stop
	end
	
	subroutine plot(iwr,ni,nj,nk,psi,cut)
	implicit real*8(a-h,o-z)
c----- write down wavefunction --------------
	parameter(ij=256,nth=175)
	complex*16 psi(ij,ij,nth)
	common /grid/ xmin,ymin,xmax,ymax,dx,dy,amuh,amud,am(2)
	common/bend/z(nth),w(nth),rot(ij,ij),tr(nth,nth)
	do 10 j=1,nj
	 do 10 i=1,ni
	  do 10 k=1,nk
	   if(abs(psi(i,j,k))**2.gt.cut) then
	    x=xmin+dx*(i-1)
	    y=ymin+dy*(j-1)
	   write(iwr,*) x,y,z(k)
	  endif
10	continue
	return
	end

	subroutine psigaus(ni,nj,nk,psi,v,pi,alfa,beta,xi,yi,pxi,pyi)
	implicit real*8(a-h,o-z)
cccccccccccc harmonic oscillator wavefunctions ccccccccccccccccc
	real*8 yr(3)
	parameter (ij=256,nth=175)
	complex*16 im,psi(ij,ij,nth)
	real*8 v(ij,ij,nth)
	common /grid/ xmin,ymin,xmax,ymax,dx,dy,amuh,amud,am(2)
	common/bend/z(nth),w(nth),rot(ij,ij),tr(nth,nth)
c------- initialize potential  ----------------
	call lsth
c----------------------------------------------
	npx=(xi-xmin)/dx
	npy=(yi-ymin)/dy
	call quadrature(nk,z,w,tr)
	if(nk.eq.1) then
	  write(*,*) 'Collinear'
	  z(nk)=1d0
	endif
	write(*,*) npx,npy
	im=(0d0,1d0)
	an=dsqrt(2d0*dsqrt(alfa*beta)/pi/2d0)
	pot=0d0
	anrm=0d0
	do 11 j=1,nj
	 y=ymin+dy*(j-1)
	 qy=y-yi
	 do 11 i=1,ni
	  x=xmin+dx*(i-1)
	  qx=x-xi
	  rot(i,j)=1d0/x/x/am(1)+1d0/y/y/am(2)
	  do 11 k=1,nk
	  psi(i,j,k)=exp(-alfa*qx**2+im*pxi*qx)*
     &               exp(-beta*qy**2+im*pyi*qy)*
     &		     dsqrt(w(k))*an
c----- >>>>
	  yr(3)=y
	  yr(1)=dsqrt(x*x+(amuh*y)**2+2d0*z(k)*x*y*amuh)
	  yr(2)=dsqrt(x*x+(amud*y)**2-2d0*z(k)*x*y*amud)
	  call potls(yr,vpot)
	  v(i,j,k)=vpot
	  wt=dx*dy
	  pot=pot+abs(psi(i,j,k))**2*vpot*wt
	  anrm=anrm+abs(psi(i,j,k))**2*wt
10	 continue
11	continue	
	write(*,*) 'nrm',anrm
	write(*,*) 'pot',pot
	write(*,*) 'kin',beta/2d0/am(2)
	return
	end

	subroutine psiprod(ni,nj,nk,psi,pi,alfa,beta,xi,yi,pxp,pyp,rp)
	implicit real*8(a-h,o-z)
cccccccccccc harmonic oscillator wavefunctions ccccccccccccccccc
	real*8 rp(2,2)
	parameter (ij=256,nth=175)
	complex*16 im,psi(ij,ij,nth)
	common /grid/ xmin,ymin,xmax,ymax,dx,dy,amuh,amud,am(2)
	common/bend/z(nth),w(nth),rot(ij,ij),tr(nth,nth)
	im=(0d0,1d0)
	an=dsqrt(2d0*dsqrt(alfa*beta)/pi/2d0)
	anrm=0d0
	do 11 j=1,nj
	 y=ymin+dy*(j-1)
	 do 11 i=1,ni
	  x=xmin+dx*(i-1)
!	  qx=x*rp(1,1)+y*rp(1,2)-xi
!	  qy=x*rp(2,1)+y*rp(2,2)-yi
	qx=0.5d0*x+0.75d0*y-xi
	qy=x-0.5d0*y-yi
	  do 10 k=1,nk
	  psi(i,j,k)=exp(-alfa*qx**2+im*pxp*qx)*
     &               exp(-beta*qy**2+im*pyp*qy)*
     &		     an*dsqrt(w(k))
	  anrm=anrm+abs(psi(i,j,k))**2*dx*dy
10	continue
11	continue	
	write(*,*) 'nrmB',anrm
	return
	end



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


        subroutine split(kmax,ni,nj,nk,N2,ak2,v,psib,psi,dt,corr)
        implicit real*8(a-h,o-z)
	complex*16 im
	parameter(ij=256,nth=175,im=(0d0,1d0))
	complex*16 c,z1,z2,z4,zx(nth),cwf(ij,ij),fx(nth),zrot,zz
        dimension ak2(ij,ij)
        dimension nn(2)
        complex*16 psib(ij,ij,nth),psi(ij,ij,nth),corr(kmax)
	real*8 v(ij,ij,nth)
        common /grid/ xmin,ymin,xmax,ymax,dx,dy,amuh,amud,am(2)
	common/bend/z(nth),w(nth),rot(ij,ij),tr(nth,nth)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	z1=-im*dt
	z2=-im*dt/2d0
	z4=-im*dt/4d0
	n2=ni*nj
        nn(1)=ni
        nn(2)=nj
	t=0d0
        corr=(0d0,0d0)
        do 11 k=1,kmax
	 t=t+dt
c------- exp(-im*V*dt/2) ------------------------
	 do 101 l=1,nk
	  do 101 j=1,nj
	   do 101 i=1,ni
101	 psi(i,j,l)=psi(i,j,l)*exp(z2*v(i,j,l)) 
	if(nk.eq.1) goto 301
c------- exp(-im*J^2*dt/2/I) --------------------
	 do 102 j=1,nj
	  do 102 i=1,ni
	   zrot=z4*rot(i,j)
	   do 100 l=1,nk
	    fx(l)=(0d0,0d0)
	    zz=exp(zrot*l*(l-1))
	    do 99 l1=1,nk
99	     fx(l)=fx(l)+zz*tr(l,l1)*psi(i,j,l1)
100	 continue
	   do 103 l=1,nk
	    zx(l)=(0d0,0d0)
	    do 104 l1=1,nk
104          zx(l)=zx(l)+tr(l1,l)*fx(l1)
103	 continue
	   do 105 l=1,nk
105	    psi(i,j,l)=zx(l)
102	continue
c--------- exp(-im*T*dt) ------------------------
301	continue
	do 106 l=1,nk
	 do 107 j=1,nj
	  do 107 i=1,ni
107	   cwf(i,j)=psi(i,j,l)
          call fourn(cwf,nn,2,1)
          call diff(cwf,ak2,dt,ni,nj)
          call fourn(cwf,nn,2,-1)
          do 12 j=1,nj
           do 12 i=1,ni
12        psi(i,j,l)=cwf(i,j)/N2
106 	continue
	if(nk.eq.1) goto 302
c----- exp(-im*J^2*dt/2/I) --------------------
	 do 202 j=1,nj
	  do 202 i=1,ni
	   zrot=z4*rot(i,j)
	   do 200 l=1,nk
	    fx(l)=(0d0,0d0)
	    zz=exp(zrot*l*(l-1))
	    do 199 l1=1,nk
199	     fx(l)=fx(l)+zz*tr(l,l1)*psi(i,j,l1)
200	 continue
	   do 203 l=1,nk
	    zx(l)=(0d0,0d0)
	    do 204 l1=1,nk
204          zx(l)=zx(l)+tr(l1,l)*fx(l1)
203	 continue
	   do 205 l=1,nk
205	    psi(i,j,l)=zx(l)
202	continue
c------- exp(-im*V*dt/2) ------------------------
302	continue
         do 201 l=1,nk
          do 201 j=1,nj
           do 201 i=1,ni
201      psi(i,j,l)=psi(i,j,l)*exp(z2*v(i,j,l))
c-------- correlation function ------------------
        do l=1,nk
        do j=1,nj
        do i=1,ni
        corr(k)=conjg(psib(i,j,l))*psi(i,j,l)*dx*dy+
     &          corr(k)
        enddo
        enddo
        enddo
	call correl(ni,nj,nk,psib,psi,poh,pod,phd,c)
	write(15,1000) t,poh,pod,phd
	write(17,1000) t,c,abs(c)
        write(18,1000) t,corr(k),abs(corr(k))
11	continue
1000    format(20(e14.7,1x))
	return
	end

	subroutine diff(cwf,ak2,ts,ni,nj)
	implicit real*8(a-h,o-z)
	parameter(nij=256) 	
	complex*16 nim,cwf(nij,nij)
	dimension ak2(nij,nij)
	nim=(0.d0,-1.d0)
	do 11 j=1,nj
	 do 11 i=1,ni
11	cwf(i,j)=cwf(i,j)*exp(nim*ts*ak2(i,j))
	return
	end

	subroutine correl(ni,nj,nk,psi0,psi,poh,pod,phd,c)
	implicit real*8(a-h,o-z)
	parameter(ij=256,nth=175,rhd=3.0d0,roh=3.0d0,rod=3.0d0)
	complex*16 psi0(ij,ij,nth),psi(ij,ij,nth),c
	real*8 yr(3)
	common /grid/ xmin,ymin,xmax,ymax,dx,dy,amuh,amud,am(2)
	common/bend/z(nth),w(nth),rot(ij,ij),tr(nth,nth)
ccc compute averages
	c=(0.d0,0.d0)
        poh=0d0
        pod=0d0
        phd=0d0
        pwt=0d0
        pds=0d0
	do 10 j=1,nj
	 y=ymin+(j-1)*dy
	 do 10 i=1,ni
	  x=xmin+(i-1)*dx
	  do 10 k=1,nk
	   c=c+psi(i,j,k)*psi(i,j,k)
	   yr(3)=y
	   yr(1)=dsqrt(x*x+(amuh*y)**2+2d0*z(k)*x*y*amuh)
	   yr(2)=dsqrt(x*x+(amud*y)**2-2d0*z(k)*x*y*amud)
	   wt=abs(psi(i,j,k))**2
           if(yr(3).ge.rhd) then
            if(yr(1).lt.roh.and.yr(2).ge.roh) poh=poh+wt
            if(yr(2).lt.roh.and.yr(1).ge.roh) pod=pod+wt
            if(yr(1).ge.roh.and.yr(2).ge.roh) pds=pds+wt
            if(yr(1).lt.roh.and.yr(2).lt.roh) pwt=pwt+wt
        endif
        if(yr(3).lt.rhd) then
         if(yr(1).lt.roh.or.yr(2).lt.roh) pwt=pwt+wt
         if(yr(1).ge.roh.and.yr(2).ge.roh) phd=phd+wt
        endif
10	continue
	c=c*dx*dy
	phd=phd*dx*dy
	poh=poh*dx*dy
	pod=pod*dx*dy
	pwt=pwt*dx*dy
	pds=pds*dx*dy
	return
	end

        subroutine aver(ni,nj,psi0,psi,c)
        implicit real*8(a-h,o-z)
        parameter(ij=256)
        complex*16 c,psi0(ij,ij),psi(ij,ij)
	common /grid/ xmin,ymin,xmax,ymax,dx,dy,amuh,amud,am(2)
	c=(0.d0,0.d0)
        do 11 j=1,nj
         do 11 i=1,ni
11      c=c+psi(i,j)*conjg(psi0(i,j))
	c=c*dx*dy
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

	subroutine quadrature(np,z,wz,tr)
	implicit real*8 (a-h,o-z)
	parameter(nth=175)
	real*8 z(nth),wz(nth),f(0:nth),p(nth,nth),tr(nth,nth)
	parameter(xmn=-1d0,xmx=1d0)
c-------- Gauss-Legendre quadrature for the bending angle -----
	call gauleg(xmn,xmx,z,wz,np)
	do 10 i=1,np
10	write(*,*) z(i),wz(i)
	write(*,*) 'P MATRIX'
	do 12 k=1,np
	 call legpol(np,z(k),f)
	 do 13 j=1,np
13	  p(j,k)=f(j-1)
12	continue
c	do 14 k=1,np
c14	write(*,*) (p(k,i),i=1,np)
c--------- transformation matrix ----------------
	write(*,*) 'T MATRIX'
	do 15 j=1,np
	 do 15 k=1,np
15	  tr(j,k)=dsqrt(wz(k)*(j-0.5d0))*p(j,k)
c	do 16 k=1,np
c16	write(*,*) (tr(k,i),i=1,np)
	return
	end

      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER*4 n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-15)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592653589793d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END

	subroutine legpol(n,x,f)
	implicit real*8 (a-h,o-z)
	real*8 f(0:n)
	f(0)=1d0
	if(n.eq.1) return
	f(1)=x
	if(n.eq.2) return
	do 10 i=3,n
	 j=i-2
10	 f(i-1)=((2d0*j+1d0)*x*f(j)-j*f(j-1))/(j+1d0)
	return
	end
	
