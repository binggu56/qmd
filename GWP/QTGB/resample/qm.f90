      program main
!----------------------------------------------------------------
!     trajectory guided gaussian basis
!     add a repulsive force between neibouring trajectories
!     functional force f = d*exp(-c*r)/r
!----------------------------------------------------------------------
      use cdat, only:al,pi,im,am,beta 

      implicit real*8(a-h,o-z)
      integer*4, parameter :: nmax=100 

      real*8, allocatable, dimension(:) :: q0,p0,s0 
      complex*16, allocatable, dimension(:) :: c0

      real*8, dimension(nmax) :: xt

      real*8,allocatable :: ke(:),v(:),x(:),u(:),p(:),w(:)
      real*8,allocatable :: du(:),qp(:),dv(:),ddv(:),s(:)
      
      complex*16,dimension(:,:),allocatable :: h,mat
     	real*8 :: f(3),ki,ncm
      complex*16, allocatable, dimension(:) :: c,dc,cold,cnew
      real*8 :: gasdev,norm 
      integer*4 :: order 

      complex*16 :: z0,zn,z,psi,z1,dpsi,cor

!      common/wave/q0,p0,a0,sigma
      common/grid/xmin,xmax,dx,np

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
      read(10,*) t,gx  
      read(10,*) kmax,dt,kout
      read(10,*) am 
      read(10,*) idum1
      read(10,*) order 
      read(10,*) nb0
      read(10,*) al,beta  

      allocate(q0(nb0),p0(nb0),s0(nb0),c0(nb0))

      read(10,*) (q0(i),i=1,nb0)
!      read(10,*) (p0(i),i=1,nb0)
!      read(10,*) (s0(i),i=1,nb0)
      read(10,*) (c0(i),i=1,nb0)
      close(10) 



      write(*,1002) 
1002  format(/'1D QTM code with Gassian Basis '//,  & 
            'analytic expression for matrix elements'//, & 
            'system: Morse',//        & 
            'Author :: Bing Gu'//)

!     grid size
      np = 1000
      xmin = -10d0
      xmax = 10d0
      dx = (xmax-xmin)/dble(np-1)

      sigma = dsqrt(1d0/2d0/al) 

!      gx = 1.8d0*sigma 
      p0 = 0d0 
      s0 = 0d0 



! --- determine number of basis for new simulation 
! --- renormalization 
      den0 = 0d0 
      do i=1,np
        xi = xmin+dx*(i-1)
        z1 = psi(xi,nb0,c0,q0,p0,s0)
        den0 = den0+conjg(z1)*z1*dx
      enddo
      write(*,1012) den0 
1012  format('initial normalization =',f10.6)

      c0 = c0/dsqrt(den0)

      xl  = minval(q0)
      nb = 1 
      xx = xl 
      z0 = psi(xl,nb0,c0,q0,p0,s0)  
      eps = 5d-3

      k = 0 
      do i=1,nmax 
        xi = -6.d0+ gx*(i-1)  
        z0 = psi(xi,nb0,c0,q0,p0,s0)
        if(abs(z0).gt.eps) then 
          k = k+1
          write(*,*) 'z0 = ',xi, abs(z0) 
          xt(k) = xi 
        endif
      enddo 

      nb = k 
      ntraj = nb 

      allocate(ke(Ntraj),v(Ntraj),x(Ntraj),p(Ntraj),s(Ntraj), &
               w(Ntraj),du(Ntraj),u(Ntraj),dv(Ntraj),ddv(Ntraj))

      allocate(mat(Ntraj,Ntraj))
      allocate(c(Ntraj),h(Ntraj,Ntraj),dc(Ntraj))
      allocate(cold(ntraj),cnew(ntraj))

      anrm = dsqrt(dsqrt(al/pi))
      c = (0d0,0d0) 
      s = 0d0 
      w = 0d0
      p = 0d0 
        
      do i=1,nb 
        x(i) = xt(i)
!        z0 = psi(x(i),nb0,c0,q0,p0,s0)
!        z1 = dpsi(x(i),nb0,c0,q0,p0,s0)
        do k=1,nb0  
          dq = x(i)-q0(k)
          z = dq/2d0 
          z0 = exp(-al*z*conjg(z))
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
      write(*,1001) sigma,gx,nb,kmax,dt,kout,am
1001  format('Initial Conditions'//, &
            'variance = ', f10.6/, &
            'traj spacing = ', f10.6/, & 
            'number of basis = ', i6/ ,   &
            'Kmax  = ', i6/ ,   &
            'dt    = ', f10.6/, &
            'kout  = ', i6/ ,   &
            'Mass  = ', f10.6/)

!--------initial coeff---------------------
      write(*,1011) 
1011  format(/'Computing inital coefficients'/)


      call multi(nb,x,p,s,mat)
!      do i=1,nb
!        write(*,1024) (mat(i,j),j=1,nb)
!      enddo
!1024  format(100(f6.4,1x))      
      call proj(nb,x,p,mat,c)

!------check normalization--------------
      prob = 0d0
      den0 = 0d0 

      do i=1,np
        xi = xmin+dx*(i-1)
        z0 = psi(xi,nb,c,x,p,s)
        prob = prob+conjg(z0)*z0*dx
        if(abs(z0).gt.1d-3) write(104,1000) xi,abs(z0)**2
      enddo
      write(*,1013)  prob
1013  format('numerical normalization =',f10.6)

      do i=1,nb
        write(110,1000) x(i),abs(c(i)),c(i)
      enddo
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
     
      call mom(nb,x,p,c,s)
 
      call ham(c,nb,x,p,s,dc)
!      call increm(nb,mat,h,c,dc)

!-----output at t=0-------------------
      write(107,1000) t,(abs(c(i)),i=1,Ntraj)
      dt = dt/kout
      dt2 = dt/2d0


      call ham(c,nb,x,p,s,dc)

      cold = c
      c = c + dc*dt

!-----begin the time step loop----------
      time: do k=1,kmax
        
        do 13 kk=1,kout
          t = t+dt

          call mom(nb,x,p,c,s)
         
! ----  half-step increments of moment, full step increment of positions
          x = x + p*dt/am 
          
          call ham(c,nb,x,p,s,dc)
          cnew = cold + dc*2d0*dt
          cold = c
          c = cnew

!-------------update dc/dt--------------------

        anm  = norm(nb,x,p,c,s)
        c = c/dsqrt(anm) 

13      enddo

! ----- correlation function 

        call corr(nb,x,p,c,s,cor)
        write(115,1000) t*2d0,cor,abs(cor)

        write(107,1000) t,(abs(c(i)),i=1,nb)

!       print out trajectories in phase space
!        do i=1,nb
!          write(101,1000) x(i),p(i)
!        enddo
        
!        write(105,1000) t,(p(i),i=1,Nb)
        write(106,1000) t,(x(i),i=1,nb0)
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
      do i=1,Np
        xi = xmin+dx*(i-1)
        z0 = psi(xi,nb,c,x,p,s) 
        if(abs(z0).gt.1d-3) write(103,1000) xi,abs(z0)**2
        prob = prob+abs(z0)**2*dx
      enddo

      an = norm(nb,x,p,c,s)
      write(*,*) 'final anal & numerical norm = ',an,prob 
      open(114,file='temp.dat')
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
      close(114) 

      write(6,1111) t 
1111  format('ENDING TIME = ',f5.2) 

1000 format(200(e14.7,1x))
      end program main

!------------------------------------------------
!      subroutine pot(x,v,ddv)
!      implicit real*8(a-h,o-z)
!
!      ipot = 1
!
!      if (ipot == 0) then
!        ep = 0.0
!        v = x**2/2d0 + ep*x**4/4.0
!        ddv = 1d0 + 3.0*x**2*ep
!      elseif(ipot == 1) then
!        eta = 1.3544d0
!        v = 1d0/16d0/eta*x**4-0.5d0*x**2
!!        dv(i) = 1d0/4d0/eta*x**3-x
!        ddv = 3d0/4d0/eta*x**2-1d0
!      endif
!
!      return
!      end subroutine

! --- compute wavefunction value at point x 

      double complex function psi(x,nb,c,q,p,s)
      
      use cdat, only : im,pi,al  
      
      implicit real*8(a-h,o-z) 
      real*8 ,intent(in) :: x 
      integer*4 :: nb 
      real*8, intent(in), dimension(nb) :: q,p,s
      complex*16,intent(in) :: c(nb) 
      
      anrm = dsqrt(dsqrt(al/pi))

      psi = (0d0, 0d0) 
      do j=1,nb 
        psi = psi+c(j)*exp(-al/2d0*(x-q(j))**2+im*s(j))*anrm 
      enddo 
      
      end function 


! ---- update momentum 

      subroutine mom(nb,x,p,c,s) 

      use cdat 

      implicit real*8(a-h,o-z) 

      real*8 x(nb),p(nb),s(nb) 
      complex*16 c(nb),z,z0,dz 

! ----- Gaussian convolution 

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

! --- derivative of psi at x       
      double complex function dpsi(x,nb,c,q,p,s) 
      
      use cdat, only : im,pi,al  
      
      implicit real*8(a-h,o-z) 
      real*8 ,intent(in) :: x 
      integer*4 :: nb 
      real*8, intent(in), dimension(nb) :: q,p,s
      complex*16,intent(in) :: c(nb) 
      complex*16 :: zj 
      
      anrm = dsqrt(dsqrt(al/pi))

      dpsi = (0d0, 0d0) 
      do j=1,nb 
        zj = exp(-al/2d0*(x-q(j))**2+im*s(j))*anrm
        dpsi = dpsi+c(j)*(-al*(x-q(j)))*zj                    
      enddo 
      
      end function 

! --- compute non-classical momemtum 

      double precision function ncm(x,nb,c,q,p,s)
 
      use cdat, only : im,pi,al 
      
      implicit real*8(a-h,o-z) 
      real*8 ,intent(in) :: x 
      integer*4 :: nb 
      real*8, intent(in),dimension(nb) :: q,p,s
      complex*16 :: c(nb),psi,z0,z1,gi,gj 
      
!      common/wave/q0,p0,a0,sigma
      anrm = dsqrt(dsqrt(al/pi))

      ncm = 0d0 
!      do i=1,nb 
!        ncm = ncm+abs(c(i))**2*anrm**2*exp(-al*(x-q(i))**2) 
!      enddo 

!      h = 0.05d0 
      z0 = psi(x,nb,c,q,p,s) 
!      xx = x+h 
!      z1 = psi(xx,nb,c,q,p,s)

      den = conjg(z0)*z0 
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

      do i=1,nb 
        do j=1,nb 
          dp = p(i)-p(j)
          gi = anrm*exp(-al/2d0*(x-q(i))**2+im*s(i))
          gj = anrm*exp(-al/2d0*(x-q(j))**2+im*s(j))
          ncm = ncm+conjg(c(i))*c(j)*(-al*(2d0*x-q(i)-q(j))-im*dp)*conjg(gi)*gj
        enddo 
      enddo 

      ncm = ncm/(2d0*den)
      end function 
        

      

!------normalization-----------------------------
      double precision function norm(nb,x,p,c,s)

      implicit real*8(a-h,o-z)

!------compute normalization--------------------

      integer*4, intent(in) :: nb 

      real*8 :: x(nb),p(nb),s(nb)

      complex*16 :: z,mat(nb,nb),c(nb)

      call multi(nb,x,p,s,mat)

      z = (0d0,0d0)
      do j=1,nb
        do i=1,nb
          z = z+conjg(c(j))*mat(j,i)*c(i)
        enddo
      enddo

      norm = real(z) 

      end function 


      subroutine corr(nb,x,p,c,s,cor)

      implicit real*8(a-h,o-z)

!------compute normalization--------------------

      integer*4, intent(in) :: nb 

      real*8 :: x(nb),p(nb),s(nb)

      complex*16 :: z,mat(nb,nb),c(nb),cor 

      call multi(nb,x,p,s,mat)

      z = (0d0,0d0)
      do j=1,nb
        do i=1,nb
          z = z+c(j)*mat(j,i)*c(i)
        enddo
      enddo

      cor = z 
      return 
      end subroutine  
!     ----------------------------------------------
!     -----------------------------------
!     gassian wavepackets overlap matrix, analytic solution
!     -----------------------------------
      subroutine multi(nb,q,p,s,mat)

      use cdat,only : im,pi,al 
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
!          z = dq/2d0-im/2d0/al*(p(j)-p(k)) 
          z = dq/2d0
!          z0 = exp(-al*z*conjg(z)+im/2d0*(p(j)+p(k))*dq)
          z0 = exp(-al*dq**2/4d0)
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

!      common/wave/q0,p0,a0,sigma
      
      call multi(nb,q,p,s,mat)
      
!      call eigen(nb,mat)

      anrm = dsqrt(dsqrt(al/pi))
      
      eta = 1.3544d0 

      do j=1,nb
        do k=1,nb
          dq = q(k)-q(j)
          z = dq/2d0
          z0 = exp(-al*dq**2/4d0)
!          d0 = z0*(v(k)-(p(k)**2-al)/2d0/am) ! without action 
          d0 = al/(4d0*am)*(1d0-al/2d0*dq**2)
!          d1 = v(k) + dv(k)*z + 0.5*ddv(k)*(z**2 + 1d0/(2d0*al))
          y = (q(j)+q(k)) 
!          call pot(y,v,ddv)
!          d1 = v + ddv/(4d0*al)
          d1 = 1d0/(16d0*eta)*(y**4/16d0+3d0*y**2/(4d0*al)+3d0/(4d0*al**2))-0.5d0*(y**2/4d0+1d0/(2d0*al))
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
        write(*,*) 'dc/dt: Matrix fails, INFO =',INFO
        stop
      endif
!------compute increment dc------------
!      do i=1,nb 
!        c(i) = c(i) + dt*dc(i)
!      enddo 
      do j=1,nb 
        dc(j) = -im*(b1(j,1)+im*b1(j,2))
      enddo 
      
      return
      end subroutine




!     -------------------------------------------
!     hamiltonian matrix (H-id/dt) with gaussian 
!     basis
!     -------------------------------------------

!     ----------------------------------------------
!     initial coefficients before gaussian basis
!     ----------------------------------------------
      subroutine proj(nb,x,p,mat,c0)
      use cdat, only:im,pi,al 
      implicit real*8(a-h,o-z)
      
      integer*4,intent(in) :: nb

      real*8,intent(in) :: x(nb),p(Nb)

      complex*16,intent(inout) :: c0(Nb)
      complex*16 :: z,z0,mat(Nb,Nb),psi0
      real*8 :: aux(nb,nb),b(nb,2)
!      common/wave/q0,p0,a0,sigma 
!      a = al
!      alfa = 2d0*a0
!      as=a+alfa
!      av=a*alfa/as
!      an=dsqrt(2d0*dsqrt(av/as))
!      an0=dsqrt(dsqrt(a/pi))
!      an2=dsqrt(dsqrt(alfa/pi))
!!------compute c_k=<g_k|psi0>---------
!      do j=1,nb
!        dq=x(j)-q0
!        dp=p(j)-p0
!        d2=-0.5d0*av*dq*dq
!        d0=-0.5d0/as*dp*dp
!        d1=(alfa*p(j)+a*p0)/as*dq
!        c0(j)=exp(d2+d0+im*d1)*an
!      enddo


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
        

!     solve matrix equation, Mc = b    
!      do i=1,nb
!        write(*,1009) (mat(i,j),j=1,Nb)
!      enddo
!1009  format(20(f5.2,1x))
!------save overlap matrix------------------
      do i=1,nb
      do j=1,nb
      aux(i,j) = real(mat(i,j))
      enddo 
      enddo 
      do i=1,nb 
        b(i,1) = real(c0(i)) 
        b(i,2) = imag(c0(i))
      enddo 

      call dposv('U', Nb, 2, aux, nb, b, nb, INFO)

      IF( INFO.ne.0 ) THEN
         WRITE(*,1222) INFO
1222     format('INFO = ', i6/ ,'projection fails.'/)
         STOP
      END IF
      
      do i=1,nb
      c0(i) = b(i,1)+im*b(i,2)
      enddo 

      return
      end subroutine


