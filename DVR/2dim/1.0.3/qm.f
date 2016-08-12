	program qm
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccc 2D models cccccccccccccccc
cccccccc use Gaussians instead of the eigenstates ccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit real*8 (a-h,o-z)
	integer*4 ni,nj
	parameter (nmx=64,nx=nmx*nmx)
	real*8 v(nmx,nmx),iso,eig(nx),ham(nx,nx)
	complex*16 im,c,psi(nmx,nmx),psib(nmx,nmx),cr(nx)



	common /prop/alfa,dt,kmax,kout
	common /params/ di,wi
	common /grid/ xmin,ymin,xmax,ymax,dx,dy,ds,am(2)
        common/eig/ac,anc,pi,vmn1,vbar,nh
c---------- constants 
	im=(0.d0,1.d0)
      
      pi=4.d0*atan(1.d0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	open(15,file='proj')
	open(16,file='qwp')
	open(17,file='cor')
	open(18,file='eig')
	open(20,file='wf0')
	open(21,file='wft')
	open(22,file='wft_gnu')
	open(54,file='vec')
	open(55,file='vec_gnu')
	open(56,file='cof')
	open(5,file='INdvr',status='old')
c time propagation parameters
	read(5,*) dt,kmax,kout
c grid spacing
	read(5,*)a3,cut 
c initial wave packet parameters
	read(5,*) alfa,beta,x0,y0,px0,py0 
c number of wavewpackets 
	read(5,*) xmin,ymin
	read(5,*) ame,iso
	read(5,*) ni,nj
	close(5)
c-------- masses -------------------------------------------
        amuh=ame
        amud=ame*iso
        am(1)=amuh
        am(2)=amud
        write(*,*) 'MASSES', am(1),am(2)
        write(*,*) 'GWP:x', alfa,x0,px0
        write(*,*) 'GWP:y', beta,y0,py0
	xmax=-xmin
	ymax=-ymin
        dx=abs(xmax-xmin)/(ni-1)
        dy=abs(ymax-ymin)/(nj-1)
	ds=dx*dy
!--------------------------------------------------------
	call psigaus(ni,nj,psi,psib,v,alfa,beta,x0,y0,px0,py0,a3)
        call ham_dvr(ni,nj,v,psi,cr,ham,eig)
	t=0d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	call plot(20,ni,nj,nk,psi,cut)
cccccccccccccccccc propagate psi ccccccccccccccccccccccccccccccccc
!	 call correl(ni,nj,cr,ham,eig)
220	continue
c	call plot(21,ni,nj,nk,psi,cut)

1000  format(20(e14.7,1x))
	stop
	end
	
	subroutine plot(iwr,ni,nj,nk,psi,cut)
	implicit real*8(a-h,o-z)
c----- write down wavefunction --------------
	parameter(nmx=64)
	complex*16 psi(nmx,nmx)
	common /grid/ xmin,ymin,xmax,ymax,dx,dy,ds,am(2)
	do 10 j=1,nj
	 do 10 i=1,ni
	   if(abs(psi(i,j))**2.gt.cut) then
	    x=xmin+dx*(i-1)
	    y=ymin+dy*(j-1)
	   write(iwr,*) x,y
	  endif
10	continue
	return
	end

	subroutine psigaus(ni,nj,psi,psib,v,alfa,beta,xi,yi,pxi,pyi,a3)
	implicit real*8(a-h,o-z)
cccccccccccc harmonic oscillator wavefunctions ccccccccccccccccc
	parameter (nmx=64)
	complex*16 im,psi(nmx,nmx),psib(nmx,nmx)
	real*8 v(nmx,nmx),gm(2,2),eig(2),work(50)
      
      real*8 :: tran(2,2) ! transformation matrix between coordinates 
      parameter(a1=15d0,a2=15d0,a4=20d0,a5=a2*a2/4d0/a1)
	
      common /grid/ xmin,ymin,xmax,ymax,dx,dy,ds,am(2)
      common/eig/ac,anc,pi,vmn1,vbar,nh
	
      im=(0d0,1d0)
      an=dsqrt(2d0/pi*dsqrt(alfa*beta))
        albe=0d0
!------------------- extrema of the potential ----------------
!---------------- barrier is at (0,0) ------------------------
        vbar=a5
!---------------- donor/acceptor wells (xmn1,ymn1) and (-xmn1,-ymn1)
        xmn1=dsqrt((4d0*a4*a2+a3*a3)/8d0/a4/a1)
        ymn1=-xmn1*a3/a4/2d0
        vmn1=xmn1**2*(a1*xmn1**2-a2)+a3*xmn1*ymn1+a4*ymn1**2+a5
        write(*,*) 'BARRIER: V(0,0)',vbar
        write(*,*) 'WELL:',vmn1
        write(*,*) 'WELL LOCATION:',xmn1,ymn1
        gx=2d0*(6d0*a1*xmn1**2-a2)
        gy=2d0*a4
        gc=a3
        write(*,*) 'DIAG WIDTH PARAMS',sqrt(am(1)*gx)/2d0,sqrt(am(2)*gy)/2d0
!-------------- construct proper Hessian and the well eigenstate -----------
        gx=gx/am(1)
        gy=gy/am(2)
        gc=gc/dsqrt(am(1)*am(2))
!-------------- construct proper Hessian and the well eigenstate -----------
        gm(1,1)=gx/am(1)
        gm(2,2)=gy/am(2)
        gm(1,2)=gc/dsqrt(am(1)*am(2))
        gm(2,1)=gm(1,2)
        call dsyev('V','L',2,gm,2,eig,work,50,info)
        wx=dsqrt(eig(1))/2d0
        wy=dsqrt(eig(2))/2d0
        write(*,*) gm(1,1),gm(1,2)
        write(*,*) gm(2,1),gm(2,2)
        gm(1,1)=gm(1,1)*dsqrt(am(1))
        gm(2,1)=gm(2,1)*dsqrt(am(1))
        gm(1,2)=gm(1,2)*dsqrt(am(2))
        gm(2,2)=gm(2,2)*dsqrt(am(2))
        alx=wx*gm(1,1)**2+wy*gm(2,1)**2
        aly=wx*gm(1,2)**2+wy*gm(2,2)**2
        axy=wx*gm(1,1)*gm(1,2)+wy*gm(2,1)*gm(2,2)
!------------- redefine psi(0) to make it the well eigenstates
        alfa=alx
        beta=aly
        albe=2d0*axy
        an=dsqrt(2d0/pi*dsqrt(alfa*beta-axy**2))
        xi=xmn1
        yi=ymn1
        write(*,*) 'WP width',alx,aly,axy
        write(*,*) 'WP center',xi,yi
!---------------- internal eigenstates info <<<<<<<<<
        nh=8
        ac=dsqrt(am(1)*ak1)/2d0
        anc=dsqrt(dsqrt(2d0*ac/pi))
!----------------------------------------->>>>>>>>>>
	pot=0d0
	anrm=0d0
      eps = 0.08 

      tran(1,1) = 1d0
      tran(1,2) = dsqrt(2d0) 
      tran(2,1) = -dsqrt(2d0)
      tran(2,2) = 1d0 

      tran = tran/dsqrt(3d0) 
        a = 1d0
        b = 4d0
        vc = 4d0


      open(77,file='pot.dat') 

	do j=1,nj
	 y=ymin+dy*(j-1)
	 qy=y-yi
	 do i=1,ni
	  x=xmin+dx*(i-1)
	  qx=x-xi
	  psi(i,j) = an*exp(-alfa*qx**2-beta*qy**2+ 
     &             im*(pyi*qy+pxi*qx)-albe*qx*qy)
          psib(i,j)=conjg(psi(i,j))
!--------------------------------------------------------
!------------- coupled double well + HO -------------------------------
!          vpot=x^2*(a1*x^2-a2)+a3*x*y+a4*y*y+a5
!test          vpot=x*x*15d0+y*y*20d0
!--------------------------------------------------------
!        vpot = x*x*(a1*x*x-a2)+a3*x*y+a4*y*y+a5
!        vpot = (x**2-4.0*dsqrt(2d0)*x*y-y**2)/6d0 +
!     +          (x**4-4.0*dsqrt(2d0)*x**3*y+12.0*x**2*y**2-
!     +           8.0*dsqrt(2d0)*x*y**3 + 4d0*y**4)/18d0 + eps*x*y
!        vpot = (1d0-(x-dsqrt(2d0)*y)**2/3d0)**2/8.0 +
!     +         (dsqrt(2d0)*x-y)**2/6d0 - 1d0/8d0
        

!        q1 = tran(1,1)*x + tran(1,2)*y 
!        q2 = tran(2,1)*x + tran(2,2)*y 
        
!        vpot = (1d0-2d0*x**2)**2/8d0 + y**2/2d0 - 1d0/8d0 + eps*q1*q2  
        v(i,j) = y**2*(a*y**2-b) + vc*(x-y)**2/2d0 + b**2/4d0/a

	  pot=pot+abs(psi(i,j))**2*vpot
	  anrm=anrm+abs(psi(i,j))**2
        write(77,*) x,y,v(i,j)
10	 continue

      enddo 
        write(77,*) ' ' 
      enddo 

      close(77) 

	write(*,*) 'nrm',anrm*ds
	write(*,*) 'pot',pot*ds
	write(*,*) 'kin',(alfa/am(1)+beta/am(2))/2d0
	return
	end

!----------------------------------------------------------
!----------------------------------------------------------
	subroutine correl(ni,nj,cr,ham,eig)
	implicit real*8(a-h,o-z)
	parameter(nmx=64,nx=nmx*nmx)
	complex*16 psi(nx),ot(0:10),c,cr(nx),tc,im,ct(nx)
	real*8 yr(3),ov(0:10),f(0:10),eig(nx),ham(nx,nx)
	common/prop/alfa,dt,kmax,kout
	common/grid/xmin,ymin,xmax,ymax,dx,dy,ds,am(2)
        common/eig/ac,anc,pi,vmn1,vbar,nh
!----------------------------------------------------------
!----------------- autocorrelation function and  probability
        write(*,*) 'COMPUTING COR FUN',dt,kmax,kout
        im=(0d0,1d0)
        nij=ni*nj
        do k=0,kmax
         c=(0d0,0d0)
         prb=0d0
         t=k*dt
         tc=-im*t
         do j=1,nij
          c=c+abs(cr(j))**2*exp(eig(j)*tc)
          ct(j)=cr(j)*exp(tc*eig(j))
         enddo
!------------ compute the wavefunction at t -----------------
         wf_nrm=0d0
         do j=1,nij
          psi(j)=(0d0,0d0)
          do i=1,nij
           psi(j)=psi(j)+ct(i)*ham(j,i)
          enddo
          wf_nrm=wf_nrm+abs(psi(j))**2
          kz=mod(k,kout)
          if(kz.eq.0) then
           ix=mod(j-1,ni)+1
           x=xmin+(ix-1)*dx
           iy=(j-1)/ni+1
           y=ymin+(iy-1)*dy
           t44=abs(psi(j))/dsqrt(ds)
           write(21,1000) x,y,t44
           write(22,1000) x,y,t44
           if(ix.eq.ni) write(22,*) 
          endif
         enddo
!------------ project onto the reactants ---------------------
         do j=1,nij
          ix=mod(j-1,ni)+1
          x=xmin+(ix-1)*dx
          if(x.lt.0d0) prb=prb+abs(psi(j))**2
         enddo
         write(17,1000) t,abs(c),prb,c,wf_nrm
        enddo
        return
!-------------------------------------------------------------
!------------------ old stuff: Hermite polynomials projections
ccc compute averages
        sc=dsqrt(2d0*alfa/pi)
	c=(0.d0,0.d0)
        xav=0d0
        yav=0d0
        x2=0d0
        y2=0d0
        pr=0d0
        do 29 n=0,nh
29       ov(n)=0d0
	do 9 j=1,nj
	 y=ymin+(j-1)*dy
         do 17 n=0,nh
17        ot(n)=(0d0,0d0)
	 do 10 i=1,ni
	  x=xmin+(i-1)*dx
!----------- generate Hermite polynomials <<<<<<<<<<<<<<
          f(0)=1d0
          f(1)=dsqrt(8d0*ac)*x
          do 18 n=2,nh
18         f(n)=2d0*dsqrt(2d0*ac)*x*f(n-1)-2d0*(n-1d0)*f(n-2) 
          fac=anc
          do 19 n=0,nh
           f(n)=f(n)*exp(-ac*x*x)*fac
19         fac=fac/dsqrt((n+1d0)*2d0)
!------------------------------->>>>>>>>>>>>>>>>>>>>>>>>
            if(x.gt.0d0) pr=pr+abs(psi(ij))**2
	    c=c+psi(ij)*psi(ij)
            xav=xav+x*abs(psi(ij))**2
            yav=yav+y*abs(psi(ij))**2
            x2=x2+x*x*abs(psi(ij))**2
            y2=y2+y*y*abs(psi(ij))**2
            do 15 n=0,nh
15           ot(n)=ot(n)+psi(ij)*f(n)
10       continue
         do 16 n=0,nh
16        ov(n)=ov(n)+abs(ot(n))**2
9	continue
	c=c*dx*dy
        xav=xav*dx*dy
        x2=x2*dx*dy
        yav=yav*dx*dy
        y2=y2*dx*dy
        pr=pr*dx*dy
!----------- write projections ---------------
        dv=dx*dx*dy
        write(15,1000) t,(ov(n)*dv,n=0,nh)
        dex=dsqrt(abs(x2-xav**2))
        dey=dsqrt(abs(y2-yav**2))
        evib=0d0
        do 20 n=0,nh
20       evib=evib+ov(n)*dv*(2*n+1)
        evib=evib*ac
        etot=evib+ekin+vxy
	write(16,1000) t,pr,xav,yav,dex,dey
c,ekin,vxy,evib,etot
1000    format(20(e14.7,1x))
	return
	end

        subroutine aver(ni,nj,psi0,psi,c)
        implicit real*8(a-h,o-z)
        parameter(nmx=64)
        complex*16 c,psi0(nmx,nmx),psi(nmx,nmx)
	common /grid/ xmin,ymin,xmax,ymax,dx,dy,ds,am(2)
	c=(0.d0,0.d0)
        do 11 j=1,nj
         do 11 i=1,ni
11      c=c+psi(i,j)*conjg(psi0(i,j))
	c=c*dx*dy
        return
        end

!========================================================================
!========================================================================
        subroutine ham_dvr(ni,nj,v,psi,cr,ham,w)
!-------------set up Hamiltonian matrix using Miller-Colbert DVR --------
        implicit real*8(a-h,o-z)
        parameter(nmx=64,nx=nmx*nmx,lwork=nx*10,eps=1d-6)
        real*8 work(lwork),w(nx),ham(nx,nx),hsv(nx,nx),v(nmx,nmx)
     &          ,ps0(nx),z(nx)
        complex*16 psi(nmx,nmx),im,wfc(nx),cr(nx),zz
        common /grid/ xmin,ymin,xmax,ymax,dx,dy,ds,am(2)
        common/eig/ac,anc,pi,vmn1,vbar,nh
1000    format(81(e14.7,1x))
        im=(0d0,1d0)
        akx=0.5d0/am(1)/dx/dx
        aky=0.5d0/am(2)/dy/dy
        p3=pi*pi/3d0
        n1=ni*nj
c----------------- rewrite potential as a vector -----------
        do j=1,nj
         do i=1,ni
          ij=i+(j-1)*ni
          z(ij)=v(i,j)
         enddo
        enddo
        write(*,*) 'MATRIX SIZE=',ij,ni*nj
c----------------- construct the hamiltonian matrix --------------
        do 10 j2=1,nj
         do 10 i2=1,ni
          j=i2+(j2-1)*ni
          do 10 j1=1,nj
           do 10 i1=1,ni
            i=i1+(j1-1)*ni
            ham(i,j)=0d0
            if(i1.eq.i2.and.j1.eq.j2) ham(i,j)=v(i2,j2)+(akx+aky)*p3
            if(j2.eq.j1.and.i1.ne.i2) then
             ham(i,j)=akx*(-1d0)**(i2-i1)*2d0/(i2-i1)**2
            endif
            if(i2.eq.i1.and.j1.ne.j2) then
             ham(i,j)=aky*(-1d0)**(j2-j1)*2d0/(j2-j1)**2
            endif
10      hsv(i,j)=ham(i,j)
        write(*,*) 'DIAGONALIZATION'
        call dsyev('V','L',n1,ham,nx,w,work,lwork,info)
        write(*,*) 'INFO',info,work(1)
        write(*,*) 'EIGENVALS',w(1),w(2),w(3)
        edel=w(2)-w(1)
        write(*,*) 'SPLITTING=',edel,edel/(w(1)-vmn1)*1d2,'% of ZPE'
        ip=100
        if(ip.gt.n1) ip=n1
        do i=1,ip
         write(18,*) i,w(i)
        enddo
c---------------------- check the eigenvector -----------------
        write(*,*) 'ENTER VECTOR TO TEST'
!        read(*,*)  nt
        nt = 3 
        eig=w(nt)
        dev=0d0
        dds=dsqrt(ds)
        do 30 j=1,n1
         ps0(j)=ham(j,nt)/dds
         if(ps0(j).lt.eps) ps0(j)=eps
         ix=mod(j-1,ni)+1
         iy=(j-1)/ni+1
         x=xmin+(ix-1)*dx
         y=ymin+(iy-1)*dy
         if (x .le. -5.99999d0) write(54,*) ' '    
         write(54,1000) x,y,ham(j,1)**2/ds,ham(j,nt)**2/ds
         
         write(55,1000) x,y,ham(j,1)**2/ds,ham(j,nt)**2/ds
         if(ix.eq.ni) write(55,*) 
         z(j)=0d0
         do 31 i=1,n1
31        z(j)=z(j)+hsv(j,i)*ham(i,nt)/eig
          dev=dev+abs(z(j)-ham(j,nt))**2
30      continue
!--------------------------------------------------------------
!--------------------  TEST the eigenvector --------------
!--------------------------------------------------------------
        write(*,*) 'EIGENVECTOR TEST',dev
        anv=0d0
        do 15 n=1,n1
15         anv=anv+abs(ham(n,nt))**2
        write(*,*) 'VECTOR INDEX and NORM=',nt,anv
!--------------------------------------------------------------
!----------------- EXPRESS psi(0) in terms of eigenvectors ----
!---------------- h is the expansion of 1-heaviside -----------
!--------------------------------------------------------------
        ctot=0d0
        do k=1,n1
         cr(k)=(0d0,0d0)
         do j=1,nj
          do i=1,ni
           ij=i+(j-1)*ni
           cr(k)=cr(k)+psi(i,j)*ham(ij,k)
          enddo
         enddo
         cr(k)=cr(k)*dds
         write(56,*) k,w(k),abs(cr(k))
         ctot=ctot+abs(cr(k))**2
        enddo
        write(*,*) 'SUM OF COEFS IS',ctot
!----------------------------------------------------------------
!----------------- plot the wavefunction ------------------------
        do j=1,nj
         y=ymin+dy*(j-1)
         do i=1,ni
          x=xmin+dx*(i-1)
          ij=i+(j-1)*ni
          zz=(0d0,0d0)
          do k=1,ni*nj
           zz=zz+cr(k)*ham(ij,k)
          enddo
          if(abs(zz).gt.1d-4) write(20,*) x,y,abs(zz)
         enddo
        enddo
        return
        end

