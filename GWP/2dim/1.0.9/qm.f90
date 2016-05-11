      module cdat

      integer*4 :: ipot

      real*8,public,parameter :: PI=4.d0*atan(1.0)

      integer,parameter :: nmax = 16

      complex*16,public,parameter::im=(0d0,1d0)

      real*8 :: al,beta

      integer*4 ::  Ntraj,kmax,kout,idum1,order,nbx,nby
      real*8 :: dt,amx,amy,ax0,ay0,ax,ay,qx0,qy0,px0,py0
      real*8 :: xmin,xmax,ymin,ymax
      save
      end module cdat


      program main
!----------------------------------------------------------------
!     trajectory guided gaussian basis
!----------------------------------------------------------------------
      use cdat
      
      implicit real*8(a-h,o-z)
       
      real*8,allocatable :: ke(:),v(:),x(:),y(:),u(:),px(:),py(:),w(:)
      real*8,allocatable :: du(:),qp(:),dv(:),ddv(:),s(:)      
      
      complex*16,dimension(:,:),allocatable :: h,mat
      complex*16, allocatable :: c(:),c0(:),dc(:),c1(:) 
      
!      complex*16, dimension(nmax,nmax) :: mat

!      complex*16, dimension(nmax) :: c,c0,dc 

     	real*8 :: f(3),ki,ncm

      real*8 :: gasdev

      complex*16 :: z0,z1,z2,zn,z,cor0,cor1,cor2,cor

!      common/wave/q0,p0,a0,sigma
!      common/grid/xmin,xmax,dx,np

      open(100,file='energy.dat')
      open(101,file='phase.dat')
      open(102,file='aver.dat')  
      open(103,file='wft.dat')
      open(104,file='wf0.dat')
      open(106,file='q.dat')
      open(105,file='p.dat')
      open(107,file='c.dat')
      open(108,file='xy0.dat')
      open(109,file='xyt.dat')
      open(110,file='cf0.dat')
      open(111,file='norm.dat')
      open(112,file='force')
      open(113,file='xc')
      open(115,file='cor.dat', status='unknown',action='write')
      open(117,file='cor2.dat')
      open(118,file='cor_asym.dat', status='unknown',action='write')
      open(119,file='cor_sym.dat', status='unknown',action='write')

      open(116,file='qforce.dat', status='unknown',action='write')


      call init_pars()
      
      ep = 1d-2 
      Nb = nbx*nby 
      ntraj = nb 
      af = 0.1 

      write(*,1002) ntraj 
1002  format(/'1D QTM code with Gassian Basis '//,  & 
              'analytic expression for matrix elements'//, & 
              'system: Morse',//        & 
              'number of basis:', i4//)
      
      allocate(y(ntraj),py(ntraj),x(nb),px(nb))

      allocate(ke(Ntraj),v(Ntraj),s(Ntraj), &
               w(Ntraj),du(Ntraj),u(Ntraj),dv(Ntraj),ddv(Ntraj))

      allocate(mat(Ntraj,Ntraj))
      allocate(c(Ntraj),c0(Ntraj),h(Ntraj,Ntraj),dc(Ntraj),c1(ntraj))
      

      pow = 4d0

!     grid size over whole range 
      npx = 100 
      npy = 100 
      xmin = -6d0
      xmax = 6d0
      dx = (xmax-xmin)/dble(npx-1)
      ymin = -6d0 
      ymax = 6d0 
      dy = (ymax-ymin)/(npy-1)

! --- grid to sample initial points    
      g = 0.8d0  
      xl = qx0-g*dsqrt(pow/ax0)
      xr = qx0+g*dsqrt(pow/ax0)
      yl = qy0-g*dsqrt(pow/ay0)
      yr = qy0+g*dsqrt(pow/ay0)

      gx = (xr-xl)/dble(nbx-1)
      gy = (yr-yl)/dble(nby-1)


      px = px0 
      py = py0 

      do i=1,nbx
        do j=1,nby
          k = nby*(i-1)+j
!          x(i) = gasdev(idum1)
!          x(i) = x(i)/dsqrt(2d0*a0)+q0
          x(k) = xl+(i-1)*gx
          y(k) = yl+(j-1)*gy
!          w(k) = exp(-2d0*a0*(x(i)-q0)**2)*dsqrt(2d0*a0/pi)*gx
      enddo 
      enddo 

! --- action function 

      s = 0d0 
      
      do i=1,nb 
        write(108,1000) x(i),y(i) 
      enddo 

! --- define gaussian width as half of bin size 
!      al = 1d0/2d0/sigma**2
      sigx = sqrt(1d0/2d0/ax) 
      sigy = sqrt(1d0/2d0/ay) 


!      xint = 1.8d0*sigma 

!     test distribution
!      an = 0d0
!      a1 = 0d0
!      do i=1,Ntraj
!        an = an+w(i)
!        a1 = a1+w(i)*x(i)
!      enddo
!      write(*,*) an,a1
!      write(*,*) xl,xr,gx

!     print out the initial conditions        
      write(*,1001) sigx,sigy,ax0,ay0,Ntraj,kmax,dt,kout, & 
      amx,amy,px0,py0,qx0,qy0,gx,gy  

1001  format('Initial Conditions'//, &
            'width = ', 2f10.6/, & 
            'a0    = ', 2f10.6/, &
            'Ntraj = ', i6/ ,   &
            'Kmax  = ', i6/ ,   &
            'dt    = ', f10.6/, &
            'kout  = ', i6/ ,   &
            'Mass  = ', 2f10.6/ ,& 
            'p0    = ', 2f10.6/ ,&
            'q0    = ', 2f10.6/, & 
            'dx,dy = ', 2f10.6)

!--------initial coeff---------------------
      write(*,1011) 
1011  format(/'Computing inital coefficients'/)


      call proj(nb,x,y,px,py,s,c0)
      
      c = c0
      do i=1,nb
        write(110,1000) x(i),y(i),abs(c(i)),c(i)
      enddo

      print *,'position = ',(x(i),i=1,6) 
      print *,'init c =',(c(i),i=1,6)  
!------check normalization--------------
      z1 = (0d0,0d0) 

      do i=1,npx
        xi = xmin+dx*(i-1)
        do j=1,npy
          yi = ymin+dy*(j-1)
          call psi(xi,yi,nb,c,x,y,s,z0)
          z1 = z1+conjg(z0)*z0*dx*dy 
!        if(abs(z0).gt.1d-3) write(104,1000) xi,abs(z0)**2
        write(104,*) xi,yi,abs(z0)**2
        enddo
        write(104,*) ' ' 
      enddo 

      write(*,1013)  real(z1) 
1013  format('initial Normalization =', f16.8)
      

!----------------------------------------------

!     print out basis functions
!      do i=1,Np 
!        xi = xmin+dx*dble(i-1)
!        z0 = dsqrt(dsqrt(al/pi))*exp(-al/2d0*(xi-x(1))**2+  &
!             im*p(1)*(xi-x(1)))
!        zn = dsqrt(dsqrt(al/pi))*exp(-al/2d0*(xi-x(2))**2+ &
!             im*p(nb)*(xi-x(2)))
!       
!       if(abs(z0)>1d-3 .or. abs(zn)>1d-3) write(108,1000) xi,abs(z0)**2,abs(zn)**2
!
!!        if(abs(z0).gt.1d-3) write(108,1000) xi,abs(z0)**2
!!        if(abs(zn).gt.1d-3) write(109,1000) xi,abs(zn)**2
!      enddo
!
! --- initial force field 

!      call derivs(x,Ntraj,v,dv,ddv)
!      call qpot(order,nb,ddv,c,am,p,x,s,w,u,du)

!      do i=1,nb
!        write(112,1000) x(i),v(i),dv(i),ddv(i),du(i)
!      enddo

!      call ham(am,w,c,nb,x,p,s,u,du,v,ddv,h)
!     call increm(nb,mat,h,c,dc)

!-----output at t=0-------------------
      t = 0.0
      write(107,1000) t,(abs(c(i)),i=1,Ntraj)
      dt = dt/kout
      dt2 = dt/2d0
      af = 0.0 !friction constant


      call ham0(c,nb,x,y,px,py,s,dc)
      
!      c1 = c 
      c = c + dc*dt2 
!      x = x + p*dt/am  

!-----begin the time step loop----------
      time: do kt = 1,kmax
        
!        print *,'time step', kt

        do 13 kk=1,kout
          
          t = t+dt

!          call norm(nb,x,p,c,mat,an)

          c = c + dc*dt2

!          call norm(nb,x,y,px,py,c,s,an)
!          c = c/dsqrt(an) 

          call mom(nb,x,y,px,py,c,s)

          x = x + px*dt/amx 
          y = y + py*dt/amy 

!-------------update dc/dt--------------------
!          call derivs(nb,x,v,dv,ddv)

        call ham0(c,nb,x,y,px,py,s,dc)

        c = c + dc*dt2


! --------- friction 
          call norm(nb,x,y,px,py,c,s,anm)
          c = c/dsqrt(anm)             
                     
13      enddo
!        call norm(nb,x,y,px,py,c,s,an)
        write(111,1000) t,an 

        call aver(nb,x,y,px,py,s,c,xav,yav)
        write(102,1000) t,xav,yav  

!        if(k == kmax) then 
!        do i=1,nb 
!          write(116,1000) x(i),du(i)
!        enddo 
!        endif 
      
! ----- correlation function 

!        call permute(nb,x,y,px,py,c,s,cor0,cor1,cor2) 
!        write(115,1000) t,cor0,cor1,cor2 

        call corr(nb,x,y,px,py,c,s,cor)
        write(117,1000) t*2d0,cor 


        write(107,1000) t,(abs(c(i)),i=1,nb)

!       print out trajectories in phase space
!        write(101,1000) (x(i),i=1,nb),(p(i),i=1,nb) 
        
        write(105,1000) t,(px(i),i=1,Nb,nby)
        write(106,1000) t,(x(i),i=1,Nb)
!       calculate the total energy
        enk = 0d0
        env = 0d0
        enq = 0d0

!        do i=1,Nb
!          env = env+v(i)*w(i)
!          enk = enk+p(i)**2/2d0/am*w(i)
!          enq = enq+u(i)*w(i)
!        enddo 

!        en = enk+env+enq

!        write(100,1000) t,enk,env,enq,en
!        call flush(100)
        
!-------------average position --------------
!        call aver(Ntraj,x,p,c,xav)
!        write(102,1000) t,xav

      end do time 

!--------print x,c--------------
        do i=1,nb
          write(113,1000) x(i),abs(c(i))  
        enddo
       
      do i=1,nb 
        write(109,1000) x(i),y(i),px(i),py(i)   
      enddo
!---------check matrix singularity--------------
!        do i=1,nb
!          write(*,1018) (mat(i,j),j=1,nb)
!        enddo
!        call eigen(nb,mat)

!1018    format(20(f6.4,1x))

!     print out final wavefunction         
      prob = 0d0

      npx = 100 
      npy = 100 

      do i=1,npx 
        xi = xmin + dx*(i-1)
        do j=1,npy 
          yi = ymin + dy*(j-1)           
        
        call psi(xi,yi,nb,c,x,y,s,z0) 
!        if(abs(z0).gt.1d-3) 
        write(103,*) xi,yi,abs(z0)**2
        prob = prob+abs(z0)**2*dx*dy 
        enddo 
        write(103,*) ' ' 
      enddo

      write(*,*) 'Final prob = ',prob

      open(114,file='temp.dat') 
      write(114,*) t        
      write(114,*) kmax,',',dt,',',kout 
      write(114,*) amx,amy 
      write(114,*) idum1
      write(114,*) gx,gy 
      write(114,*) nb 
      write(114,*) ax,ay,beta 
      
      write(114,*) (x(i),i=1,nb)
      write(114,*) (y(i),i=1,nb)
!      write(114,*) (p(i),i=1,nb)
!      write(114,*) (s(i),i=1,nb)
      write(114,*) (c(i),i=1,nb)
      close(114) 

      close(111)
      close(112)
      close(113)  
      close(115) 

1000 format(2000(e14.7,1x))
      end program main

! ---- autocorrelation using symmetrized wavefunction --------

      subroutine permute(nb,x,y,px,py,c,s,cor0,cor1,cor2) 

      use cdat, only: im,pi,ax,ay,px0,py0,qx0,qy0,ax0,ay0 
      
      implicit real*8(a-h,o-z) 
      real*8,dimension(nb) :: x,y,px,py,s 
      complex*16,dimension(nb) :: c
      complex*16, dimension(nb) ::  z0,z2,z1,z3
      complex*16, dimension(nb,nb) :: mat 
      complex*16 z,cor0,cor1,cor2 

      a = ax
      alfx = 2d0*ax0
      alfy = 2d0*ay0 
      asx = ax+alfx
      asy = ay+alfy 

      avx=ax*alfx/asx
      avy=ax*alfy/asy
      
      an0=dsqrt(2d0*dsqrt(avx/asx))*dsqrt(2d0*dsqrt(avy/asy))

!      an0=dsqrt(dsqrt(a/pi))
!      an2=dsqrt(dsqrt(alfa/pi))

!---- compute c0_k = <psi0|g_k>
      do j=1,nb
        dx = qx0-x(j)
        dy = qy0-y(j)
        dpx = px0+im*ax*dx  
        dpy = py0+im*ay*dy 
!        z = -0.5d0*(dpx**2/asx+dpy**2/asy) - 0.5d0*(ax*dx*dx+ay*dy**2)
        d2 = -0.5d0*avx*dx**2-0.5d0*avy*dy**2
!        d0=-0.5d0/as*dp*dp
!        d1=(a*p0)/as*dq
!        c0(j)=exp(d2+d0+im*d1)*an
        z0(j) = an0*exp(d2)
      enddo

!---- compute c0_k = <permuted psi0|g_k>
      alfx = 2d0*ay0
      alfy = 2d0*ax0 
      asx = ax+alfx
      asy = ay+alfy 

      avx=ax*alfx/asx
      avy=ax*alfy/asy
      
      an1=dsqrt(2d0*dsqrt(avx/asx))*dsqrt(2d0*dsqrt(avy/asy))

      do j=1,nb
        dx = qy0-x(j)
        dy = qx0-y(j)
        dpx = py0+im*ax*dx  
        dpy = px0+im*ay*dy 
!        z = -0.5d0*(dpx**2/asx+dpy**2/asy) -0.5d0*(ax*dx*dx+ay*dy*dy)
      d2 = -0.5d0*avx*dx**2-0.5d0*avy*dy**2
        z1(j) = an1*exp(d2)
      enddo

      cor0 = (0d0,0d0)
      cor1 = (0d0,0d0)
      cor2 = (0d0,0d0)

      do j=1,nb 
        cor0 = cor0 + z0(j)*exp(im*s(j))*c(j)
        cor1 = cor1 + (z0(j)+z1(j))*exp(im*s(j))*c(j)
        cor2 = cor2 + (z0(j)-z1(j))*exp(im*s(j))*c(j) 
      enddo 

      return       
      end subroutine 

! ---- update momentum 

      subroutine mom(nb,x,y,px,py,c,s)

      use cdat, only : ax,ay,pi,beta

      implicit real*8(a-h,o-z) 

      real*8,dimension(nb) ::  x,y,px,py,s 
      complex*16 c(nb),z,z0,dz1,dz2  


! ----- Gaussian convolution 

      an0 = dsqrt(dsqrt(ax*ay/pi**2))
!      alfax = ax 
!      alfay = ay 

      an0 = an0*dsqrt(beta/(beta+ax))*dsqrt(beta/(beta+ay))
      alfax = beta*ax/(beta+ax)
      alfay = beta*ay/(beta+ay)


      do j=1,nb 
        z = (0d0,0d0) 
        dz1 = (0d0,0d0) 
        dz2 = (0d0,0d0)
!        s1 = x(j) 
!        s2 = y(j)
        do i=1,nb 
          dx = x(j)-x(i) 
          dy = y(j)-y(i)
          z0 = c(i)*exp(-0.5*(alfax*dx**2+alfay*dy**2))*an0 
          z = z + z0
          dz1 = dz1 - alfax*dx*z0
          dz2 = dz2 - alfay*dy*z0  
        enddo 
        px(j) = imag(dz1/z)
        py(j) = imag(dz2/z)  
      enddo 

      return 
      end subroutine        
!------------------------------------------------


! --- autocorrelation 
!      double complex function corr(nb,c0,c,x,y,px,py,s)
!
!      use cdat, only: im,pi,ax,ay,px0,py0,qx0,qy0,ax0,ay0  
!
!      implicit real*8 (a-h,o-z) 
!
!      real*8, intent(in),dimension(nb) :: x,y,px,py,s
!      complex*16, intent(in) :: c(nb),c0(nb) 
!
!      complex*16 :: mat(nb,nb),z0,z(nb),dpx,dpy
!!      common/wave/q0,p0,a0,sigma 
!      
!      a = ax
!      alfx = 2d0*ax0
!      alfy = 2d0*ay0 
!      asx = ax+alfx
!      asy = ay+alfy 
!
!      avx=ax*alfx/asx
!      avy=ax*alfy/asy
!      
!      an=dsqrt(2d0*dsqrt(avx/asx))*dsqrt(2d0*dsqrt(avy/asy))
!
!      an0=dsqrt(dsqrt(a/pi))
!      an2=dsqrt(dsqrt(alfa/pi))
!!---- compute c0_k = <psi0|g_k>
!      do j=1,nb
!        dx = qx0-x(j)
!        dy = qy0-y(j)
!        dpx = px0+im*ax*dx  
!        dpy = py0+im*ay*dy 
!        z(j) = -0.5d0*(dpx**2/asx+dpy**2/asy) -0.5d0*(ax*dx*dx+ay*dy**2)
!!        d0=-0.5d0/as*dp*dp
!!        d1=(a*p0)/as*dq
!!        c0(j)=exp(d2+d0+im*d1)*an
!        z(j) = an*exp(z(j))
!      enddo
!
!      corr = (0d0,0d0)
!      do j=1,nb 
!          corr = corr + z(j)*exp(im*s(j))*c(j)
!      enddo 
!      
!      end function 

! --- compute wavefunction value at point x 

      subroutine psi(xi,yi,nb,c,x,y,s,z)
      
      use cdat, only : im,pi,ax,ay
      
      implicit real*8(a-h,o-z) 
      real*8 ,intent(in) :: xi,yi 
      integer*4 :: nb 
      real*8, intent(in), dimension(nb) :: x,y,s
      complex*16,intent(in) :: c(nb) 

      complex*16, intent(out) :: z 
      
      anrm = dsqrt(dsqrt(ax*ay/pi/pi))

      z = (0d0, 0d0) 
      do j=1,nb 
        z = z + c(j)*exp(-ax/2d0*(xi-x(j))**2-ay/2d0*(yi-y(j))**2+im*s(j))*anrm 
      enddo 
      
      return  
      end subroutine

! --- compute non-classical momemtum 

!      double precision function ncm(x,nb,c,q,p,s)
 
!      use cdat, only : im,pi,al 
      
!      implicit real*8(a-h,o-z) 
!      real*8 ,intent(in) :: x 
!      integer*4 :: nb 
!      real*8, intent(in),dimension(nb) :: q,p,s
!      complex*16 :: c(nb),psi,dpsi,z0,z1,gi,gj 
      
!      anrm = dsqrt(dsqrt(al/pi))

!      ncm = 0d0 
!      do i=1,nb 
!        ncm = ncm+abs(c(i))**2*anrm**2*exp(-al*(x-q(i))**2) 
!      enddo 

!      h = 1d-4 
!      z0 = psi(x,nb,c,q,p,s) 
!      xx = x+h 
!      z1 = psi(xx,nb,c,q,p,s)
!      xx = x-h 
!      z2 = psi(xx,nb,c,q,p,s) 
!      ncm = (abs(z2)-abs(z1))/(2d0*h*abs(z0))

!      den = conjg(z0)*z0 
!      den1 = conjg(z1)*z1 
!      ncm = (den1-den0)/(2d0*h*den0) 
            
!      do i=2,nb 
!        do j=1,i-1 
!          dq = q(i)-q(j)
!          dp = p(i)-p(j)
!          xs = (q(i)+q(j))/2d0 
!          ps = (p(i)+p(j))/2d0 
!          ncm = ncm+anrm**2*real(conjg(c(i))*c(j)*(-al*(2d0*x-q(i)-q(j))-im*dp) & 
!                *exp(-al*(x-xs)-im*dp*(x-xs)+im*ps*dq+im*(s(j)-s(i))))
!        enddo 
!      enddo 
!      ncm = ncm/den 
!      do i=1,nb 
!        do j=1,nb 
!          dp = p(i)-p(j)
!          gi = anrm*exp(-al/2d0*(x-q(i))**2+im*p(i)*(x-q(i))+im*s(i))
!          gj = anrm*exp(-al/2d0*(x-q(j))**2+im*p(j)*(x-q(j))+im*s(j))
!          ncm = ncm+conjg(c(i))*c(j)*(-al*(2d0*x-q(i)-q(j))-im*dp)*conjg(gi)*gj
!        enddo 
!      enddo 

!      ncm = ncm/(2d0*den)

! --- analytical way to compute r(x)     
!      z0 = psi(x,nb,c,q,p,s)
!      z1 = dpsi(x,nb,c,q,p,s)
 
!      ncm = real(conjg(z0)*z1/abs(z0)**2)
      
!      end function 

      subroutine dpsi(x,nb,c,q,s,z)
      
      use cdat, only : im,pi,al 
      
      implicit real*8(a-h,o-z) 

      real*8 ,intent(in) :: x 
      integer*4 :: nb 
      real*8, intent(in), dimension(nb) :: q,s
      complex*16,intent(in) :: c(nb) 
      complex*16 :: gj 

      complex*16, intent(out) :: z 
      
      anrm = dsqrt(dsqrt(al/pi))

      z = (0d0, 0d0) 
      do j=1,nb 
        gj = exp(-al/2d0*(x-q(j))**2+im*s(j))*anrm 
        z = z + c(j)*(-al*(x-q(j)))*gj                    
      enddo 
      
      end subroutine
      

!------normalizati-----------------------------
      subroutine norm(nb,x,y,px,py,c,s,anrm)
      
      use cdat, only : pi, al, im  

      implicit real*8(a-h,o-z)

      integer*4, intent(in) :: nb 

      real*8, dimension(nb) :: x,y,px,py,s 

      complex*16, intent(in) :: c(nb) 

      complex*16 :: mat(nb,nb),z

      real*8, intent(out) :: anrm 

      z = (0d0,0d0)

      call multi(nb,x,y,px,py,s,mat)
      

      do j=1,nb
        do k=1,nb
          z = z + conjg(c(j))*mat(j,k)*c(k)
        enddo
      enddo
      
      anrm = real(z)

      return 
      end subroutine

      subroutine corr(nb,x,y,px,py,c,s,z)
      
      use cdat, only : pi, al, im  

      implicit real*8(a-h,o-z)

      integer*4, intent(in) :: nb 

      real*8, dimension(nb) :: x,y,px,py,s 

      complex*16, intent(in) :: c(nb) 

      complex*16 :: mat(nb,nb)

      complex*16, intent(out) :: z 

      z = (0d0,0d0)

      call multi(nb,x,y,px,py,s,mat)
      

      do j=1,nb
        do k=1,nb
          z = z + c(j)*mat(j,k)*c(k)
        enddo
      enddo

      return 
      end subroutine
!     ----------------------------------------------
!     average x 
!     ---------------------------------------------
      subroutine aver(nb,x,y,px,py,s,c,xav,yav)
      use cdat 

      implicit real*8(a-h,o-z)
      integer*4,intent(in) :: nb
      real*8,intent(in),dimension(nb):: x,y,px,py,s
      complex*16, intent(in) :: c(nb)
      real*8,intent(out) :: xav
      complex*16 :: psi,z
!      common/wave/q0,p0,a0,sigma
!      common/grid/xmin,xmax,dx,np
            
      z0 = (0d0,0d0)
      z1 = (0d0,0d0)
      do j=1,nb 
        do k=1,nb
          dx = x(j)-x(k)
          dy = y(j)-y(k) 
          xjk = 0.5d0*(x(j)+x(k))*exp(-ax*dx**2/4d0-ay*dy**2/4d0)
          yjk = 0.5d0*(y(j)+y(k))*exp(-ax*dx**2/4d0-ay*dy**2/4d0) 
          z0 = z0 + conjg(c(j))*xjk*c(k) 
          z1 = z1 + conjg(c(j))*yjk*c(k) 
        enddo   
      enddo

      xav = real(z0)
      yav = real(z1)

      return
      end subroutine
!     -----------------------------------
!     gassian wavepackets overlap matrix, analytic solution
!     -----------------------------------
      subroutine multi(nb,x,y,px,py,s,mat)

      use cdat, only : im,pi,ax,ay 

      implicit real*8(a-h,o-z)
      integer*4, intent(in) :: nb
      real*8   , intent(in) :: x(nb),y(nb),px(nb),py(nb),s(nb)
      complex*16,intent(out) :: mat(nb,nb)
      complex*16 :: z,z0
!      common/wave/q0,p0,a0,sigma 


      anrm = dsqrt(dsqrt(ax*ay/pi**2))
      do j=1,nb
        do k=1,j 
          dqx = x(j)-x(k)
          dqy = y(j)-y(k)
          zx = dqx/2d0
          zy = dqy/2d0
          z0 = exp(-ax*zx**2-ay*zy**2)
          mat(j,k) = z0*exp(im*(s(k)-s(j))) 
        enddo
      enddo

      do j=2,nb 
        do k=1,j-1 
          mat(k,j) = conjg(mat(j,k))
        enddo 
      enddo 
      
      return 
      end subroutine
! --------------------------------------------
      subroutine pot(x,y,v0) 
      implicit real*8(a-h,o-z) 

!      ipot = 0 

!      if (ipot == 0) then 
        ep = 0.2d0

!        v0 = x**2/2d0 + y**2/2d0+ep*x**4/4
!        v11 = 1d0+3d0*ep*x**2 
!        v22 = 1d0 
      v0 = 0.5*(x**2+y**2)+1d0/abs(x+y) 
!      v11 = -4d0/x**3+2d0/(x+y)**3
!      v22 = -4d0/y**3+2d0/(x+y)**3 


!      elseif(ipot == 1) then 
!        eta = 1.3544d0
!        v = 1d0/16d0/eta*x**4-0.5d0*x**2
!        dv(i) = 1d0/4d0/eta*x**3-x
!        ddv = 3d0/4d0/eta*x**2-1d0
!      endif 

      return 
      end subroutine 

      subroutine ham0(c,nb,x,y,px,py,s,dc)
      
      use cdat 
      
      implicit real*8(a-h,o-z)
      
      integer*4, intent(in) :: nb
      complex*16, dimension(nb), intent(in) :: c  
      real*8, intent(in), dimension(nb) :: x,y,px,py,s 
      
      complex*16,intent(out) :: dc(nb)     

      complex*16, dimension(nb,nb) :: h,mat     
      real*8 :: aux(nb,nb) 

      complex*16 :: z,z0,d0,d1,d2,b(nb),d3 
      real*8 b1(nb,2)
     
      call multi(nb,x,y,px,py,s,mat)

      anrm = dsqrt(dsqrt(ax*ay/pi**2))
      
      do j=1,nb
        do k=1,nb
          dx = x(k) - x(j)
          dy = y(k) - y(j)
!          z = dq/2d0
          z0 = exp((-ax*dx**2-ay*dy**2)/4d0)
          d0 = ax/(4d0*amx)*(1d0-ax/2d0*dx**2)+ay/(4d0*amy)*(1d0-ay/2d0*dy**2) ! kinetic 
          xc = (x(j)+x(k)) 
          yc = (y(j)+y(k))
!          call pot(xc,yc,v0,v11,v22)
          
         d1 = 0.5*xc**2-0.5*yc**2  + 4d0 - xc*yc + yc**4/16d0 +  & 
             (3d0/4d0)*yc**2/ay + 1d0/ax -1d0/ay + 3d0/(4d0*ay**2)
!     -----------------------
!          d1 = 0.5d0*(xc**2 + 0.5/ax + yc**2 + 0.5/ay) ! potential 
!          xt = abs(x(j)+x(k))/dsqrt(2d0*ax)
!          yt = abs(y(j)+y(k))/dsqrt(2d0*ay)
!          dxy = (x(j)+x(k))/dsqrt(2d0*ax) - (y(j)+y(k))/dsqrt(2d0*ay) ! need check when ax != ay 
!          dxy = abs(dxy) 
!          d1 = -2d0*(erf(dsqrt(ax)*xt)/xt + erf(dsqrt(ay)*yt)/yt) + erf(dsqrt(ax*ay/(ax+ay))*dxy)/dxy 
!          d3 = erf(dsqrt(ax*ay/(ax+ay))*dxy)/dxy  ! Columb potential
!     --------------------
          d2 = im*px(k)*ax/(2d0*amx)*dx + im*py(k)*ay/(2d0*amy)*dy 
          h(j,k) = z0*(d0+d1+d2) 
        enddo
      enddo
      
      b = matmul(h,c) 

      do i=1,nb 
        b1(i,1) = real(b(i))
        b1(i,2) = imag(b(i))
      enddo 
      
!      call inverse(mat,mativs,nb)


!-------save the overlap matrix-----------
      do i=1,nb
      do j=1,nb
        aux(i,j) = real(mat(i,j))
      enddo 
      enddo 

      call dposv('U', nb, 2, aux, nb, b1, nb, INFO)
      if(INFO .ne. 0) then
        write(*,*) 'dC/dt: Matrix fails, INFO =',INFO
        stop
      endif
!------compute increment dc------------
      do j=1,nb 
        dc(j) = -im*(b1(j,1)+im*b1(j,2))
      enddo 
      
      return
      end subroutine


!     ------------------------------------------

      subroutine increm(nb,mat,h,c,dc)
      use cdat, only:im,pi
      implicit real*8(a-h,o-z)
      complex*16,intent(in) :: mat(nb,nb),h(nb,nb),c(nb)
      complex*16,intent(out) :: dc(nb)
      complex*16 b(nb),aux(nb,nb)
!---------H*c---------------
!      do i=1,Nb
!        b(i) = (0d0,0d0)
!        do j=1,Nb
!          b(i) = b(i)+h(i,j)*c(j)
!        enddo
!      enddo
      b = matmul(h,c)
!      b = -im*b 

!-------save the overlap matrix-----------
      aux = mat
!      do i=1,Nb
!        write(*,1023) (aux(i,j),j=1,nb)
!      enddo
!      write(*,*) (b(i),i=1,Nb)
!1023  format(100(f6.4,1x))
      call zposv('U', Nb, 1, aux, nb, b, nb, INFO)
      if(INFO .ne. 0) then
        write(*,*) 'dC/dt: Matrix fails, INFO =',INFO
        stop
      endif
!------compute increment dc------------
      do j=1,nb
        dc(j) = -im*b(j)
      enddo
!      dc = b 

      return
      end subroutine
!     ----------------------------------------------
!     initial coefficients before gaussian basis
!     ----------------------------------------------
      subroutine proj(nb,x,y,px,py,s,c)
      
      use cdat  

      implicit real*8(a-h,o-z)
      integer*4,intent(in) :: nb
      real*8,intent(in),dimension(nb) :: x,y,px,py,s

      complex*16,intent(out) :: c(nb) 
      
      complex*16 :: c0(nb) 
      complex*16 :: z,z0,mat(nb,nb),psi0,aux(nb,nb)

!      common/wave/q0,p0,a0,sigma 
!      a = ax
      alfx = 2d0*ax0
      alfy = 2d0*ay0 
      asx = ax+alfx
      asy = ay+alfy
      avx = ax*alfx/asx
      avy = ay*alfy/asy 
      an = dsqrt(2d0*dsqrt(avx/asx)*2d0*dsqrt(avy/asy))
      an0=dsqrt(dsqrt(ax/pi)*dsqrt(ay/pi))
      an2=dsqrt(dsqrt(alfx/pi)*dsqrt(alfy/pi))
!------compute c_k=<g_k|psi0>---------
      do j=1,nb
        dx = x(j)-qx0
        dy = y(j)-qy0
        dpx = px(j)- px0
        dpy = py(j)-py0 

        d2=-0.5d0*avx*dx**2-0.5d0*avy*dy**2
        d0=-0.5d0/asx*dpx**2-0.5d0/asy*dpy**2
        d1=(alfx*px(j)+ax*px0)/asx*dx + (alfy*py(j)+ay*py0)/asy*dy  
        c0(j)=exp(d2+d0+im*d1)*an
      enddo


!     intial wavepacket is also a gaussian (p0,q0) with the same width 
!     as basis 
!      do j=1,nb
!        dq = q(j)-q0
!        z = dq/2d0-im/2d0/al*(p(j)-p0)
!        z0 = exp(-al*z*conjg(z)+im/2d0*(p(j)+p0)*dq)
!        c0(j) = z0*dsqrt(dsqrt(al/pi))
!      enddo

!     a0 not equal to al

!      do j=1,nb
!        z = (0d0,0d0)
!        do i=1,Np
!          x = xmin+dx*(i-1)
!          z0 = exp(-al/2d0*(x-q(j))**2+im*p(j)*(x-q(j)))
!          psi0 = exp(-a0/2d0*(x-q0)**2+im*p0*(x-q0))
!          z = z+conjg(z0)*psi0*dx*dsqrt(dsqrt(a0*al)/pi)
!        enddo
!        c0(j) = z
!      enddo
        

!     solve matrix equation, Mc = b    
!      do i=1,nb
!        write(*,1009) (mat(i,j),j=1,Nb)
!      enddo
!1009  format(20(f5.2,1x))
!------save overlap matrix------------------

      call multi(nb,x,y,px,py,s,mat) 

      aux = mat

      call zposv('U', Nb, 1, aux, nb, c0, nb, INFO)

      IF( INFO.ne.0 ) THEN
         WRITE(*,1222) INFO
1222     format('INFO = ', i6/ ,'projection fails.'/)
         STOP
      END IF

      do i=1,nb 
        c(i) = c0(i) 
      enddo 
      
      return
      end subroutine

      
