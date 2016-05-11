      program main
!----------------------------------------------------------------
!     trajectory guided gaussian basis
!----------------------------------------------------------------------
      use cdat
      
      implicit real*8(a-h,o-z)
       
!      real*8, dimension(nmax) :: ke,v,x,u,p,w,du,qp,dv,ddv,s,q,s0 

      real*8,allocatable :: ke(:),v(:),x(:),u(:),p(:),w(:)
      real*8,allocatable :: du(:),qp(:),dv(:),ddv(:),s(:)      
      
      complex*16,dimension(:,:),allocatable :: h,mat
      complex*16, allocatable :: c(:),c0(:),dc(:),c1(:)
      complex (kind = 8), allocatable ::  cold(:), cnew(:)  
      
!      complex*16, dimension(nmax,nmax) :: mat

!      complex*16, dimension(nmax) :: c,c0,dc 

     	real*8 :: f(3),ki,ncm

      real*8 :: gasdev

      complex*16 :: z0,z1,z2,zn,z,cor 

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
      open(111,file='norm.dat')
      open(112,file='force')
      open(113,file='xc')
      open(114,file='temp.dat',status='unknown',action='write')
      open(115,file='cor.dat', status='unknown',action='write')
      open(116,file='qforce.dat', status='unknown',action='write')
      open(117,file = 'pes.dat') 


      call init_pars()
      
      ep = 1d-2 
      Nb = Ntraj
      Np = 800
      af = 0.1 

      write(*,1002) ntraj 
1002  format(/'1D QTM code with Gassian Basis '//,  & 
              'analytic expression for matrix elements'//, & 
              'system: Morse',//        & 
              'number of basis:', i4//)

      allocate(ke(Ntraj),v(Ntraj),x(Ntraj),p(Ntraj),s(Ntraj), &
               w(Ntraj),du(Ntraj),u(Ntraj),dv(Ntraj),ddv(Ntraj))

      allocate(mat(Ntraj,Ntraj))
      allocate(c(Ntraj),c0(Ntraj),h(Ntraj,Ntraj),dc(Ntraj),c1(ntraj))
      
      allocate(cold(ntraj),cnew(ntraj))

      pow = 6d0

!     grid size over whole range 

      xmin = -10d0
      xmax = 10d0
      dx = (xmax-xmin)/dble(np-1)

! --- plot the potential 
      do i=1,np 
        q = xmin + dx*(i-1) 
        call pot(q,vpot,addv) 
        write(117,1000) q,vpot
      enddo 
      

! --- grid to sample initial points      
      xl = q0-1.d0*dsqrt(pow/a0)
      xr = q0+1.d0*dsqrt(pow/a0)
      gx = (xr-xl)/dble(nb-1)

      do i=1,nb 
!          x(i) = gasdev(idum1)
!          x(i) = x(i)/dsqrt(2d0*a0)+q0
          p(i) = 0d0 
          x(i) = xl+(i-1)*gx
          w(i) = exp(-2d0*a0*(x(i)-q0)**2)*dsqrt(2d0*a0/pi)*gx
      enddo 

! --- action function 

      s = 0d0 

! --- define gaussian width as half of bin size 
!      al = 1d0/2d0/sigma**2
      sigma = sqrt(1d0/2d0/al) 

!      xint = 1.8d0*sigma 

!     test distribution
      an = 0d0
      a1 = 0d0
      do i=1,Ntraj
        an = an+w(i)
        a1 = a1+w(i)*x(i)
      enddo
      write(*,*) an,a1
      write(*,*) xl,xr,gx

!     print out the initial conditions        
      write(*,1001) sigma,a0,Ntraj,kmax,dt,kout,am,p0,q0
1001  format('Initial Conditions'//, &
            'width = ', f10.6/, & 
            'a0    = ', f10.6/, &
            'Ntraj = ', i6/ ,   &
            'Kmax  = ', i6/ ,   &
            'dt    = ', f10.6/, &
            'kout  = ', i6/ ,   &
            'Mass  = ', f10.6/ ,& 
            'p0    = ', f10.6/ ,&
            'q0    = ', f10.6//)

!--------initial coeff---------------------
      write(*,1011) 
1011  format(/'Computing inital coefficients'/)


      call proj(nb,x,p,s,c0)
      
      c = c0
      do i=1,nb
        write(110,1000) x(i),abs(c(i)),c(i)
      enddo
      write(6,1012) x  
1012  format('initial sampling position = '//,20(f5.2,1x)/)  

      write(6,1014) p
1014  format('initial sampling momentum = '//,20(f5.2,1x)/)  

      write(6,1015) c  
1015  format('initial expansion coeff   = '//,20('(',f5.2,',',f5.2,')',1x)/)  

!------check normalization--------------
      z1 = (0d0,0d0) 
      do i=1,np
        xi = xmin+dx*(i-1)
        call psi(xi,nb,c,x,s,z0)
        z1 = z1+conjg(z0)*z0*dx
        if(abs(z0).gt.1d-3) write(104,1000) xi,abs(z0)**2
      enddo
      
      call norm(nb,x,p,c,s,anrm)
      write(*,1013)  anrm,real(z1)   
1013  format('initial normalization =', 2f12.5/)
      


!----------------------------------------------

!     print out basis functions
      do i=1,Np 
        xi = xmin+dx*dble(i-1)
        z0 = dsqrt(dsqrt(al/pi))*exp(-al/2d0*(xi-x(1))**2+  &
             im*p(1)*(xi-x(1)))
        zn = dsqrt(dsqrt(al/pi))*exp(-al/2d0*(xi-x(2))**2+ &
             im*p(nb)*(xi-x(2)))
       
       if(abs(z0)>1d-3 .or. abs(zn)>1d-3) write(108,1000) xi,abs(z0)**2,abs(zn)**2

!        if(abs(z0).gt.1d-3) write(108,1000) xi,abs(z0)**2
!        if(abs(zn).gt.1d-3) write(109,1000) xi,abs(zn)**2
      enddo

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


      call ham0(c,nb,x,p,s,dc)
      
!      c1 = c
      cold = c  
      c = c + dc*dt
!      x = x + p*dt/am  

!-----begin the time step loop----------
      time: do kt = 1,kmax
        
!        print *,'time step', kt

        do 13 kk=1,kout
          
          t = t+dt

!          call norm(nb,x,p,c,mat,an)

!          c = c + dc*dt2

          call mom(nb,x,p,c,s)

          x = x + p*dt/am 

!-------------update dc/dt--------------------
!          call derivs(nb,x,v,dv,ddv)

        call ham0(c,nb,x,p,s,dc)

        cnew = cold + dc*2d0*dt 
        cold = c 
        c = cnew 
        


! --------- friction 
          call norm(nb,x,p,c,s,anm)
          c = c/dsqrt(anm)             
                     
13      enddo

!        if(k == kmax) then 
!        do i=1,nb 
!          write(116,1000) x(i),du(i)
!        enddo 
!        endif 
      
! ----- correlation function 

        call corr(nb,x,p,c,s,cor)

        write(115,1000) 2d0*t,cor,abs(cor) 

        write(107,1000) t,(abs(c(i)),i=1,nb)

!       print out trajectories in phase space
        write(101,1000) (x(i),i=1,nb),(p(i),i=1,nb) 
        
        write(105,1000) t,(p(i),i=1,Nb)
        write(106,1000) t,(x(i),i=1,nb)
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
!---------check matrix singularity--------------
!        do i=1,nb
!          write(*,1018) (mat(i,j),j=1,nb)
!        enddo
!        call eigen(nb,mat)

!1018    format(20(f6.4,1x))

!     print out final wavefunction         
      prob = 0d0
      do i=1,Np
        xi = xmin + dx*(i-1)
        call psi(xi,nb,c,x,s,z0) 
        if(abs(z0).gt.1d-3) write(103,1000) xi,abs(z0)**2
        prob = prob+abs(z0)**2*dx
      enddo

      write(*,*) 'Final prob = ',prob
      write(114,*) t,gx  
      write(114,*) kmax,',',dt,',',kout 
      write(114,*) am 
      write(114,*) idum1
      write(114,*) order 
      write(114,*) nb
      write(114,*) al,beta   
      
      write(114,*) (x(i),i=1,nb)
!      write(114,*) (p(i),i=1,nb)
!      write(114,*) (s(i),i=1,nb)
      write(114,*) (c(i),i=1,nb)

      close(111)
      close(112)
      close(113) 
      close(114) 
      close(115) 

1000 format(2000(e14.7,1x))
      end program main

! ---- update momentum 

      subroutine mom(nb,x,p,c,s) 

      use cdat, only : al,pi,beta

      implicit real*8(a-h,o-z) 

      real*8 x(nb),p(nb),s(nb) 
      complex*16 c(nb),z,z0,dz 


! ----- Gaussian convolution 

!      beta = 16d0 

      an0 = dsqrt(dsqrt(al/pi))
      an0 = an0*dsqrt(beta/(beta+al))
      alfa = beta*al/(beta+al)

      do j=1,nb 
        z = (0d0,0d0) 
        dz = (0d0,0d0) 
        q = x(j) 
        do i=1,nb 
          y = q-x(i) 
          z0 = c(i)*exp(-0.5*alfa*y**2)*an0 
          z = z + z0
          dz = dz - alfa*y*z0 
        enddo 
        p(j) = imag(dz/z) 
      enddo 

      return 
      end subroutine        
!------------------------------------------------


! --- autocorrelation 
!      double complex function corr(nb,c,q,p,s)
!
!      use cdat
!
!      implicit real*8 (a-h,o-z) 
!
!      real*8, intent(in),dimension(nb) :: q,p,s
!      complex*16, intent(in) :: c(nb) 
!
!      complex*16 :: mat(nb,nb),z,z0,c0(nb)
!!      common/wave/q0,p0,a0,sigma 
!      
!      a = al
!      alfa = 2d0*a0
!      as=a+alfa
!      av=a*alfa/as
!      an=dsqrt(2d0*dsqrt(av/as))
!      an0=dsqrt(dsqrt(a/pi))
!      an2=dsqrt(dsqrt(alfa/pi))
!!---- compute c0_k = <psi0|g_k>
!      do j=1,nb
!        dq=q(j)-q0
!!        dp=p(j)-p0  
!        dp = -p0
!        d2=-0.5d0*av*dq*dq
!        d0=-0.5d0/as*dp*dp
!        d1=(a*p0)/as*dq
!        c0(j)=exp(d2+d0+im*d1)*an
!      enddo
!
!      corr = (0d0,0d0)
!      do i=1,nb 
!        corr = corr+conjg(c0(i))*exp(im*s(i))*c(i)
!      enddo 
!      
!      end function 

! --- compute wavefunction value at point x 

      subroutine psi(xx,nb,c,q,s,z)
      
      use cdat, only : im,pi,al
      
      implicit real*8(a-h,o-z) 
      real*8 ,intent(in) :: xx
      integer*4 :: nb 
      real*8, intent(in), dimension(nb) :: q,s
      complex*16,intent(in) :: c(nb) 

      complex*16, intent(out) :: z 
      
      anrm = dsqrt(dsqrt(al/pi))

      z = (0d0, 0d0) 
      do j=1,nb 
        z = z + c(j)*exp(-al/2d0*(xx-q(j))**2+im*s(j))*anrm 
      enddo 
      
      end subroutine

! --- compute non-classical momemtum 

!      double precision function ncm(x,nb,c,q,p,s)
 
!      use cdat, only : im,pi,al 
      
!      implicit real*8(a-h,o-z) 
!      real*8 ,intent(in) :: x 
!      integer*4 :: nb 
!      real*8, intent(in),dimension(nb) :: q,p,s
!      complex*16 :: c(nb),psi,dpsi,z0,z1,gi,gj 
      
!      common/wave/q0,p0,a0,sigma
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
      subroutine norm(nb,x,p,c,s,anrm)
      
      use cdat, only : pi, al, im  

      implicit real*8(a-h,o-z)

      integer*4, intent(in) :: nb 

      real*8, dimension(nb) :: x,p,s 

      complex*16, intent(in) :: c(nb) 

      complex*16 :: mat(nb,nb),z

      real*8, intent(out) :: anrm 

      z = (0d0,0d0)

      call multi(nb,x,p,s,mat)
      

      do j=1,nb
        do k=1,nb
          z = z + conjg(c(j))*mat(j,k)*c(k)
        enddo
      enddo

      anrm = abs(z)

      return 
      end subroutine

      subroutine corr(nb,x,p,c,s,cor)
      
      use cdat, only : pi, al, im  

      implicit real*8(a-h,o-z)

      integer*4, intent(in) :: nb 

      real*8, dimension(nb) :: x,p,s 

      complex*16, intent(in) :: c(nb) 

      complex*16 :: mat(nb,nb),z,cor

      z = (0d0,0d0)

      call multi(nb,x,p,s,mat)

      do j=1,nb
        do k=1,nb
          z = z + c(j)*mat(j,k)*c(k)
        enddo
      enddo

      cor = z 

      return 
      end subroutine

!     ----------------------------------------------
!     average
!     ---------------------------------------------
!      subroutine aver(q,p,s,c,xav)
!      use cdat 
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
!          psi = psi+c(j)*exp(-al/2d0*(xi-q(j))**2)*dsqrt(dsqrt(al/pi))
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
      subroutine multi(nb,q,p,s,mat)

      use cdat, only : im,pi,al 

      implicit real*8(a-h,o-z)
      integer*4, intent(in) :: nb
      real*8   , intent(in) :: q(nb),p(nb),s(nb)
      complex*16,intent(out) :: mat(nb,nb)
      complex*16 :: z,z0
!      common/wave/q0,p0,a0,sigma 


      anrm = dsqrt(dsqrt(al/pi))
      do j=1,nb
        do k=1,j 
          dq = q(j)-q(k)
          z = dq/2d0
          z0 = exp(-al*z*conjg(z))
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
      subroutine pot(x,v,ddv) 
      implicit real*8(a-h,o-z) 

      ipot = 2 

      if (ipot == 0) then 
        ep = 0.0
        v = x**2/2d0 + ep*x**4/4.0
        ddv = 1d0 + 3.0*x**2*ep  

      elseif(ipot == 1) then 

        eta = 1.3544d0
        v = 1d0/16d0/eta*x**4-0.5d0*x**2
        ddv = 3d0/4d0/eta*x**2-1d0
      elseif(ipot == 2) then 
      
        D = 16d0 
!        a = 1.3624d0 
        a = 1d0 
        v = D/cosh(a*x)**2
        ddv = 6d0*D*a**2*sinh(a*x)**2/cosh(a*x)**4-2d0*D*a**2/cosh(a*x)**2
      endif 

      return 
      end subroutine 

      subroutine ham0(c,nb,q,p,s,dc)
      
      use cdat 
      
      implicit real*8(a-h,o-z)
      
      integer*4, intent(in) :: nb
      complex*16, dimension(nb), intent(in) :: c  
      real*8, intent(in), dimension(nb) :: q,p,s 
      
      complex*16,intent(out) :: dc(nb)     

      real*8 :: v,dv,ddv

      complex*16, dimension(nb,nb) :: h,mat     
      real*8 :: aux(nb,nb) 

      complex*16 :: z,z0,d0,d1,d2,b(nb)
      real*8 b1(nb,2)
     
      call multi(nb,q,p,s,mat)

      anrm = dsqrt(dsqrt(al/pi))
      eta = 1.3544d0 
      
      do j=1,nb
        do k=1,nb
          dq = q(k)-q(j)
          z = dq/2d0
          z0 = exp(-al*dq**2/4d0)
          d0 = al/(4d0*am)*(1d0-al/2d0*dq**2)
          y = (q(j)+q(k))/2d0

! --- local harmonic approximation 
          call pot(y,v,ddv)
          d1 = v + ddv/(4d0*al)

! --- analytical potential integration for double-well 
!          d1 = 1d0/(16d0*eta)*(y**4 + 3d0*y**2/al + 3d0/(4d0*al**2)) - 0.5d0*(y**2+1d0/(2d0*al))

          d2 = im*p(k)*al/(2d0*am)*dq 
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


!     -------------------------------------------
!     hamiltonian matrix (H-id/dt) with gaussian 
!     basis
!     -------------------------------------------
      subroutine ham(c0,c,nb,q,p,s)
      
      use cdat, only:pi,im,al,dt 
      
      implicit real*8(a-h,o-z)
      
      integer*4, intent(in) :: nb
      complex*16, dimension(nb) :: c,c0   
      real*8, intent(in), dimension(nb) :: q,p,s 
      
!      complex*16,intent(out) :: dc(nb)     

      real*8 :: v,dv,ddv

      complex*16, dimension(nb,nb) :: h,mat     
      real*8 :: aux(nb,nb) 

      complex*16 :: z,z0,d0,d1,d2,b(nb)
      real*8 b1(nb,2)
     
!      common/wave/q0,p0,a0,sigma
      
      call multi(nb,q,p,s,mat)

      anrm = dsqrt(dsqrt(al/pi))
      
      do j=1,nb
        do k=1,nb
          dq = q(k)-q(j)
          z = dq/2d0
          z0 = exp(-al*dq**2/4d0)
!          d0 = z0*(v(k)-(p(k)**2-al)/2d0/am) ! without action 
          d0 = al/(4d0*am)*(1d0-al/2d0*dq**2)
!          d1 = v(k) + dv(k)*z + 0.5*ddv(k)*(z**2 + 1d0/(2d0*al))
          y = (q(j)+q(k))/2d0 
          call pot(y,v,ddv)
          d1 = v + ddv/(4d0*al)
          d2 = im*p(k)*al/(2d0*am)*dq 
          h(j,k) = z0*(d0+d1+d2) 

        enddo
      enddo

!      print *,'h = ', h(1,2),h(3,4) 
!      print *,'mat =', mat(1,2),mat(3,4)
      
!      b = (0d0,0d0) 
!      do i=1,nb 
!        do j=1,nb 
!          b(i) = b(i) + h(i,j)*c(j)
!        enddo 
!      enddo
      b = -2d0*dt*im*matmul(h,c) + matmul(mat,c0) 

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
      
      c0 = c ! update c0 
 
      do j=1,nb 
        c(j) = -im*(b1(j,1)+im*b1(j,2))
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
      subroutine proj(nb,x,p,s,c)
      
      use cdat  

      implicit real*8(a-h,o-z)
      integer*4,intent(in) :: nb
      real*8,intent(in),dimension(nb) :: x,p,s

      complex*16,intent(out) :: c(nb) 
      
      complex*16 :: c0(nb) 
      complex*16 :: z,z0,mat(nb,nb),psi0,aux(nb,nb)

!      common/wave/q0,p0,a0,sigma 
      a = al
      alfa = 2d0*a0
      as=a+alfa
      av=a*alfa/as
      an=dsqrt(2d0*dsqrt(av/as))
      an0=dsqrt(dsqrt(a/pi))
      an2=dsqrt(dsqrt(alfa/pi))
!------compute c_k=<g_k|psi0>---------
      do j=1,nb
        dq=x(j)-q0
        dp=p(j)-p0
        d2=-0.5d0*av*dq*dq
        d0=-0.5d0/as*dp*dp
        d1=(alfa*p(j)+a*p0)/as*dq
        c0(j)=exp(d2+d0+im*d1)*an
      enddo


!      write(*,*) 'overlap',mat(1,1),mat(1,2),mat(1,3)
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
        

!------save overlap matrix------------------


      call multi(nb,x,p,s,mat) 

      aux = mat

!     solve matrix equation, Mc = b    
!      do i=1,nb
!        write(*,1009) (mat(i,j),j=1,Nb)
!      enddo
!1009  format(20('(',f5.2,',',1x,f5.2,')',2x))
      
      write(6,1010) c0
1010  format('c = ',20('(',f5.2,',',1x,f5.2,')',2x))

      call zposv('U', Nb, 1, aux, nb, c0, nb, INFO)

      IF( INFO.ne.0 ) THEN
         WRITE(*,1222) INFO
1222     format('INFO = ', i6/ ,'Initial projection fails.'/)
         STOP
      END IF

      do i=1,nb 
        c(i) = c0(i) 
      enddo 
      
      return
      end subroutine

      

!      subroutine reproj(nb,x,p,s,c)
!      use cdat
!
!      implicit real*8(a-h,o-z)
!      integer*4,intent(in) :: nb
!      real*8,intent(in),dimension(nmax) :: x,p,s
!
!      complex*16,intent(inout) :: c(nmax) 
!      
!      complex*16 :: c0(nb) 
!
!      complex*16 :: z,z0,mat(Nb,Nb),psi0,aux(nb,nb)
!
!      do i=1,nb 
!        c0(i) = c(i) 
!      enddo 
!
!      call multi(nb,x,p,s,mat) 
!
!      aux = mat
!
!      call zposv('U', Nb, 1, aux, nb, c0, nb, INFO)
!
!      IF( INFO.ne.0 ) THEN
!         WRITE(*,1222) INFO
!1222     format('INFO = ', i6/ ,'projection fails.'/)
!         STOP
!      END IF
!
!      do i=1,nb 
!        c(i) = c0(i) 
!      enddo 
!      
!      return
!      end subroutine


