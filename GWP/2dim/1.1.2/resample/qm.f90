      module cdat

      integer*4 :: ipot

      real*8,public,parameter :: PI=4.d0*atan(1.0)

      integer,parameter :: nx = 40, ny = 40, nmax = nx*ny  

      complex*16,public,parameter::im=(0d0,1d0)

      real*8 :: ax,ay,beta

      integer*4 ::  Ntraj,kmax,kout,idum1,nbx,nby
      real*8 :: dt,amx,amy,ax0,ay0
      real*8 :: xmin,xmax,ymin,ymax
      save
      end module cdat

      program main
!----------------------------------------------------------------
!     trajectory guided gaussian basis
!     add a repulsive force between neibouring trajectories
!     functional force f = d*exp(-c*r)/r
!----------------------------------------------------------------------
      use cdat

      implicit real*8(a-h,o-z)

      real*8, allocatable, dimension(:) :: qx0,px0,qy0,py0,s0 
      complex*16, allocatable, dimension(:) :: c0

      real*8, dimension(nmax) :: xt,yt 

      real*8,allocatable :: ke(:),v(:),x(:),y(:),u(:),px(:),py(:),w(:)
      real*8,allocatable :: du(:),qp(:),dv(:),ddv(:),s(:) 
      
      complex*16,dimension(:,:),allocatable :: h,mat
     	real*8 :: f(3),ki,ncm
      complex*16, allocatable :: c(:),dc(:)
      real*8 :: gasdev 

      complex*16 :: z0,zn,z,z1,cor,z2 

!      common/wave/q0,p0,a0,sigma
!      common/grid/xmin,xmax,dx,np

      open(100,file='energy.dat')
      open(101,file='phase.dat')
      open(102,file='aver')  
      open(103,file='wft.dat')
      open(104,file='wf0.dat')
      open(106,file='q.dat')
      open(105,file='p.dat')
      open(107,file='c.dat')
      open(108,file='basis.dat')
!      open(109,file='basisn.dat')
      open(110,file='cf0.dat')
      open(111,file='norm')
      open(112,file='force')
      open(113,file='xc')
      open(115,file='cor.dat',action='write',access='append')
      
      open(10,file='temp.dat',status='old',action='read')
      read(10,*) t  
      read(10,*) kmax,dt,kout
      read(10,*) amx,amy 
      read(10,*) idum1
      read(10,*) gx,gy  
      read(10,*) nb0 
      read(10,*) ax,ay,beta 

      allocate(qx0(nb0),qy0(nb0),px0(nb0),py0(nb0),s0(nb0),c0(nb0))

      read(10,*) (qx0(i),i=1,nb0)
      read(10,*) (qy0(i),i=1,nb0)
      read(10,*) (c0(i),i=1,nb0)
      close(10) 



      write(*,1002) 
1002  format(/'Description : 2D QTM code with Gassian Basis, & 
              analytic expression for matrix elements'//, & 
              'Author      : Bing Gu'//)

            
            
!     grid size over whole range 
      npx = 100 
      npy = 100 

      xmin = -6d0
      xmax = 6d0
      dx = (xmax-xmin)/dble(npx-1)
      ymin = -6d0 
      ymax = 6d0 
      dy = (ymax-ymin)/(npy-1)

      sigx = dsqrt(1d0/2d0/ax) 
      sigy = dsqrt(1d0/2d0/ay)

      p0 = 0d0 
      s0 = 0d0 



! --- determine number of basis for new simulation 
! --- renormalization 
!       den0 = 0d0 
!       do i=1,np
!         xi = xmin+dx*(i-1)
!         z1 = psi(xi,nb0,c0,q0,p0,s0)
!         den0 = den0+conjg(z1)*z1*dx
!       enddo
!       write(*,1012) den0 
! 1012  format('initial normalization =',f10.6)

!      c0 = c0/dsqrt(den0)

      xl = minval(qx0)
      yl = minval(qy0)

!       xx = xl 
!       z0 = psi(xl,nb0,c0,q0,p0,s0)  
      eps = 5d-3
      k = 0 
      do i=1,nx 
        xi = xl + gx*(i-1)  
        do j=1,ny
          yj = yl + gy*(j-1) 
          call psi(xi,yj,nb0,c0,qx0,qy0,s0,z0)

        if(abs(z0).gt.eps) then 
          k = k+1
          write(*,*) 'z0 = ',xi, yj, abs(z0) 
          xt(k) = xi
          yt(k) = yj 
        endif

        enddo 
      enddo 

      nb = k 
      ntraj = nb 

      allocate(ke(Ntraj),v(Ntraj),x(Ntraj),y(nb),px(Ntraj),py(nb),s(Ntraj), &
               w(Ntraj),du(Ntraj),u(Ntraj),dv(Ntraj),ddv(Ntraj))

      allocate(mat(Ntraj,Ntraj))
      allocate(c(Ntraj),h(Ntraj,Ntraj),dc(Ntraj))

      anrm = dsqrt(dsqrt(al/pi))
      c = (0d0,0d0) 
      s = 0d0 
      w = 0d0
        
      do i=1,nb 
        x(i) = xt(i)
        y(i) = yt(i) 
        px(i) = 0d0 
        py(i) = 0d0 
        do k=1,nb0  
          dqx = x(i)-qx0(k)
          dqy = y(i)-qy0(k)
          z1 = dqx/2d0 
          z2 = dqy/2d0
          z0 = exp(-ax*z1*conjg(z1)-ay*z2*conjg(z2))
          c(i) = c(i) + z0*exp(im*s0(k))*c0(k)
        enddo
      enddo 
      
!      xl = q0-1.2d0*dsqrt(pow/a0)
!      xr = q0+1.2d0*dsqrt(pow/a0)
!      gx = (xr-xl)/dble(Ntraj-1)

!      do i=1,Ntraj
!          x(i) = gasdev(idum1)
!          x(i) = x(i)/dsqrt(2d0*a0)+q0
!          p(i) = p0
!          x(i) = xl+dble(i-1)*gx
!          w(i) = exp(-2d0*a0*(x(i)-q0)**2)*dsqrt(2d0*a0/pi)*gx
!      enddo 

! --- action function 
      
! --- define gaussian width as half of bin size 
!     test distribution
!      an = 0d0
!      a1 = 0d0
!      do i=1,Ntraj
!        an = an+w(i)
!        a1 = a1+w(i)*x(i)
!      enddo
!      write(*,*) an,a1
!      write(*,*) xl,gx

!     print out the initial conditions        
      write(*,1001) sigx,sigy,gx,nb,kmax,dt,kout,am
1001  format('Initial Conditions'//, &
            'variance = ', 2f10.6/, &
            'traj spacing = ', f10.6/, & 
            'number of basis = ', i6/ ,   &
            'Kmax  = ', i6/ ,   &
            'dt    = ', f10.6/, &
            'kout  = ', i6/ ,   &
            'Mass  = ', f10.6/)

!--------initial coeff---------------------
      write(*,1011) 
1011  format(/'Computing inital coefficients'/)



!      do i=1,nb
!        write(*,1024) (mat(i,j),j=1,nb)
!      enddo
!1024  format(100(f6.4,1x))      
      
      call proj(nb,x,y,px,py,s,c)
!------check normalization--------------
      prob = 0d0
      den0 = 0d0 

!------check normalization--------------
      z1 = (0d0,0d0) 

      do i=1,npx
        xi = xmin+dx*(i-1)
        do j=1,npy
          yi = ymin+dy*(j-1)
          call psi(xi,yi,nb,c,x,y,s,z0)
          z1 = z1+conjg(z0)*z0*dx*dy !  if(abs(z0).gt.1d-3) write(104,1000) xi,abs(z0)**2
        write(104,*) xi,yi,abs(z0)**2
        enddo
        write(104,*) ' ' 
      enddo 

      write(*,1013)  real(z1) 
1013  format('initial Normalization =', f16.8)


!-----output at t=0-------------------
      write(107,1000) t,(abs(c(i)),i=1,Ntraj)
      dt = dt/kout
      dt2 = dt/2d0
      
      call ham0(c,nb,x,y,px,py,s,dc)
      
!      c1 = c 
      c = c + dc*dt2 
!      x = x + p*dt/am  

!-----begin the time step loop----------
      time: do kt = 1,kmax
     

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


! ----- correlation function 
        call corr(nb,x,y,px,py,c,s,cor)
        write(115,1000) t*2d0,cor 

        write(107,1000) t,(abs(c(i)),i=1,Ntraj)

!       print out trajectories in phase space
!        do i=1,nb
!          write(101,1000) x(i),p(i)
!        enddo
        
!        write(105,1000) t,(p(i),i=1,Nb)
        write(106,1000) t,(x(i),i=1,Nb)
!       calculate the total energy
!        enk = 0d0
!        env = 0d0
!        enq = 0d0

!      write(*,*), 'u=',u
!        do i=1,Nb
!          env = env+v(i)*w(i)
!          enk = enk+p(i)**2/2d0/am*w(i)
!          enq = enq+u(i)*w(i)
!        enddo 

!        en = enk+env+enq
!        write(100,1000) t,enk,env,enq,en
!        call flush(100)
        

!        write(111,1000) t,an 
!-------------average position --------------
!        call aver(Ntraj,x,p,c,xav)
!        write(102,1000) t,xav

      end do time 

!--------print x,c--------------
        do i=1,nb
          write(113,1000) x(i),c(i)
        enddo
!---------check matrix singularity--------------
!        do i=1,nb
!          write(*,1018) (mat(i,j),j=1,nb)
!        enddo
!        call eigen(nb,mat)

!1018    format(20(f7.4,1x))
!     print out final wavefunction         
      prob = 0d0

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


1000 format(200(e14.7,1x))
      end program main

      
      
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
      
      
      
      
      
!------potential energy ------------------------------------------
!       subroutine pot(x,v,ddv)
!       implicit real*8(a-h,o-z)
! 
!       ipot = 1
! 
!       if (ipot == 0) then
!         ep = 0.0
!         v = x**2/2d0 + ep*x**4/4.0
!         ddv = 1d0 + 3.0*x**2*ep
!       elseif(ipot == 1) then
!         eta = 1.3544d0
!         v = 1d0/16d0/eta*x**4-0.5d0*x**2
! !        dv(i) = 1d0/4d0/eta*x**3-x
!         ddv = 3d0/4d0/eta*x**2-1d0
!       endif
! 
!       return
!       end subroutine

! --- compute wavefunction value at point x


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
  
!-------- derivative of psi at x       
!      double complex function dpsi(x,nb,c,q,p,s) 
!      
!      use cdat, only : im,pi,al  
!      
!      implicit real*8(a-h,o-z) 
!      real*8 ,intent(in) :: x 
!      integer*4 :: nb 
!      real*8, intent(in), dimension(nb) :: q,p,s
!      complex*16,intent(in) :: c(nb) 
!      complex*16 :: gj 
!      
!      anrm = dsqrt(dsqrt(al/pi))
!
!      dpsi = (0d0, 0d0) 
!      do j=1,nb 
!        gj = exp(-al/2d0*(x-q(j))**2+im*s(j))*anrm 
!        dpsi = dpsi+c(j)*(-al*(x-q(j)))*gj                    
!      enddo 
!      
!      end function 

! --- compute non-classical momemtum 

!      double precision function ncm(x,nb,c,q,p,s)
! 
!      use cdat, only : im,pi,al 
!      
!      implicit real*8(a-h,o-z) 
!      real*8 ,intent(in) :: x 
!      integer*4 :: nb 
!      real*8, intent(in),dimension(nb) :: q,p,s
!      complex*16 :: c(nb),psi,z0,z1,gi,gj 
!      
!!      common/wave/q0,p0,a0,sigma
!      anrm = dsqrt(dsqrt(al/pi))
!
!      ncm = 0d0 
!!      do i=1,nb 
!!        ncm = ncm+abs(c(i))**2*anrm**2*exp(-al*(x-q(i))**2) 
!!      enddo 
!
!!      h = 0.05d0 
!      z0 = psi(x,nb,c,q,p,s) 
!!      xx = x+h 
!!      z1 = psi(xx,nb,c,q,p,s)
!
!      den = conjg(z0)*z0 
!!      den1 = conjg(z1)*z1 
!!      ncm = (den1-den0)/(2d0*h*den0) 
!            
!!      do i=2,nb 
!!        do j=1,i-1 
!!          dq = q(i)-q(j)
!!          dp = p(i)-p(j)
!!          xs = (q(i)+q(j))/2d0 
!!          ps = (p(i)+p(j))/2d0 
!!          ncm = ncm+anrm**2*real(conjg(c(i))*c(j)*(-al*(2d0*x-q(i)-q(j))-im*dp) & 
!!                *exp(-al*(x-xs)-im*dp*(x-xs)+im*ps*dq+im*(s(j)-s(i))))
!!        enddo 
!!      enddo 
!!      ncm = ncm/den 
!
!      do i=1,nb 
!        do j=1,nb 
!          dp = p(i)-p(j)
!          gi = anrm*exp(-al/2d0*(x-q(i))**2+im*s(i))
!          gj = anrm*exp(-al/2d0*(x-q(j))**2+im*s(j))
!          ncm = ncm+conjg(c(i))*c(j)*(-al*(2d0*x-q(i)-q(j))-im*dp)*conjg(gi)*gj
!        enddo 
!      enddo 
!
!      ncm = ncm/(2d0*den)
!      end function 
        

      subroutine norm(nb,x,y,px,py,c,s,anrm)
      
      use cdat, only : pi, ax,ay, im  

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
      
      use cdat, only : pi, ax,ay, im  

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
!     average
!     ---------------------------------------------
!      subroutine aver(Ntraj,q,p,s,c,xav)
!      use cdat, only: im,pi,al 
!
!      implicit real*8(a-h,o-z)
!      integer*4,intent(in) :: Ntraj
!      real*8,intent(in)::q(Ntraj),p(Ntraj),s(ntraj)
!      complex*16, intent(in) :: c(Ntraj)
!      real*8,intent(out) :: xav
!      complex*16 :: psi,z
!!      common/wave/q0,p0,a0,sigma
!      common/grid/xmin,xmax,dx,np
!            
!      z = (0d0,0d0)
!      do i=1,Np
!        xi = xmin+dx*(i-1)
!        psi = (0d0,0d0)
!        do j=1,ntraj
!          psi = psi+c(j)*exp(-al/2d0*(xi-q(j))**2+ &
!                im*s(j))*dsqrt(dsqrt(al/pi))
!        enddo
!        z = z+conjg(psi)*psi*dx*xi
!      enddo
!
!      xav = real(z)
!
!      return
!      end subroutine

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
!       alfx = 2d0*ax0
!       alfy = 2d0*ay0 
!       asx = ax+alfx
!       asy = ay+alfy
!       avx = ax*alfx/asx
!       avy = ay*alfy/asy 
!       an = dsqrt(2d0*dsqrt(avx/asx)*2d0*dsqrt(avy/asy))
!       an0=dsqrt(dsqrt(ax/pi)*dsqrt(ay/pi))
!       an2=dsqrt(dsqrt(alfx/pi)*dsqrt(alfy/pi))
! !------compute c_k=<g_k|psi0>---------
!       do j=1,nb
!         dx = x(j)-qx0
!         dy = y(j)-qy0
!         dpx = px(j)- px0
!         dpy = py(j)-py0 
! 
!         d2=-0.5d0*avx*dx**2-0.5d0*avy*dy**2
!         d0=-0.5d0/asx*dpx**2-0.5d0/asy*dpy**2
!         d1=(alfx*px(j)+ax*px0)/asx*dx + (alfy*py(j)+ay*py0)/asy*dy  
!         c0(j)=exp(d2+d0+im*d1)*an
!       enddo


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

      call zposv('U', Nb, 1, aux, nb, c, nb, INFO)

      IF( INFO.ne.0 ) THEN
         WRITE(*,1222) INFO
1222     format('INFO = ', i6/ ,'projection fails.'/)
         STOP
      END IF

      
      return
      end subroutine

!     -------------------------------------------
!     hamiltonian matrix (H-id/dt) with gaussian 
!     basis
!     -------------------------------------------

