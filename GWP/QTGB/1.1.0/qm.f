      program main
!----------------------------------------------------------------
!     trajectory guided gaussian basis
!     add a repulsive force between neibouring trajectories
!     functional force f = d*exp(-c*r)/r
!----------------------------------------------------------------------
      use cdat
      implicit real*8(a-h,o-z)
      real*8,allocatable :: ke(:),v(:),x(:),u(:),p(:),w(:)
      real*8,allocatable :: du(:),qp(:),dv(:),
     &                      ddv(:),s(:)
      complex*16,allocatable :: h(:,:),mat(:,:),psi(:)
     	real*8 :: f(3),ki
      complex*16,allocatable :: c(:),c0(:),dc(:)
      real :: gasdev
      complex*16 :: z0,zn
      common/wave/al,q0,p0,a0
      common/grid/xmin,xmax,dx,np

      open(100,file='energy.dat')
      open(101,file='phase.dat')
      open(102,file='aver')  
      open(103,file='wf.dat')
      open(104,file='wf0.dat')
      open(106,file='q.dat')
      open(105,file='p.dat')
      open(107,file='c.dat')
      open(108,file='basis1.dat')
      open(109,file='basisn.dat')
      open(110,file='cf0.dat')
      open(111,file='norm')
      open(112,file='force')
      open(113,file='xc')

      open(5,file='IN')
      read(5,*) Ntraj
      read(5,*) kmax,dt,kout
      read(5,*) am
      read(5,*) al,a0
      read(5,*) q0
      read(5,*) p0
      read(5,*) idum1
      close(5)

      Np = 1000

      write(*,1002) Ntraj
1002  format(/'1D QTM code with Gassian Basis '//,
     +        'analytic expression for matrix elements'//,
     +        'system: Morse',//
     +        'number of basis:', i4//,
     +        'version: 1.4'//)
      allocate(ke(Ntraj),v(Ntraj),x(Ntraj),
     & p(Ntraj),s(Ntraj),w(Ntraj),du(Ntraj),
     & u(Ntraj),dv(Ntraj),ddv(Ntraj))

      allocate(mat(Ntraj,Ntraj),psi(np))
      allocate(c(Ntraj),c0(Ntraj),h(Ntraj,Ntraj),dc(Ntraj))
      
      Nb = Ntraj
      pow = 6d0
!     grid size
      xmin = -10d0
      xmax = 10d0

      dx = (xmax-xmin)/dble(np-1)
     
      xl = q0-1.2d0*dsqrt(pow/a0)
      xr = q0+1.2d0*dsqrt(pow/a0)
      gx = (xr-xl)/dble(Ntraj-1)

      do i=1,Ntraj
!          x(i) = gasdev(idum1)
!          x(i) = x(i)/dsqrt(2d0*a0)+q0
          p(i) = p0
          x(i) = xl+dble(i-1)*gx
          w(i) = exp(-2d0*a0*(x(i)-q0)**2)*dsqrt(2d0*a0/pi)*gx
      enddo 

!      w(1) = w(1)/2d0
!      w(nb) = w(nb)/2d0

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
      write(*,1001) al,a0,Ntraj,kmax,dt,kout,am,p0,q0
1001  format('Initial Conditions'// ,
     +       'al    = ', f10.6/,
     +       'a0    = ', f10.6/,
     +       'Ntraj = ', i6/ , 
     +       'Kmax  = ', i6/ ,
     +       'dt    = ', f10.6/ ,
     +       'kout  = ', i6/ , 
     +       'Mass  = ', f10.6/ , 
     +       'p0    = ', f10.6/ ,
     +       'q0    = ', f10.6/)

!--------initial coeff---------------------
      write(*,1011) 
1011  format(/'Computing inital coefficients'/)
      call multi(nb,x,p,mat)
!      do i=1,nb
!        write(*,1024) (mat(i,j),j=1,nb)
!      enddo
!1024  format(20(f6.4,1x))      
      call proj(nb,x,p,mat,c0)
      do i=1,nb
        c(i) = c0(i)
      enddo

      do i=1,nb
        write(110,10000) x(i),abs(c(i)),c(i)
      enddo
!----------------------------------------------

!     print out basis functions
      do i=1,Np 
        xi = xmin+dx*dble(i-1)
        z0 = dsqrt(dsqrt(al/pi))*exp(-al/2d0*(xi-x(1))**2+
     +       im*p(1)*(xi-x(1)))
        zn = dsqrt(dsqrt(al/pi))*exp(-al/2d0*(xi-x(nb))**2+
     +       im*p(nb)*(xi-x(nb)))
        if(abs(z0).gt.1d-3) write(108,10000) xi,abs(z0)**2
        if(abs(zn).gt.1d-3) write(109,10000) xi,abs(zn)**2
      enddo

      call derivs(x,Ntraj,v,dv,ddv)
      call qpot(Nb,ddv,c,am,p,x,w,u,du)
      do i=1,nb
        write(112,10000) x(i),v(i),dv(i),ddv(i),du(i)
      enddo

      call ham(am,w,c,nb,x,p,du,v,ddv,h)
      call increm(nb,mat,h,c,dc)

!------check normalization--------------
      prob = 0d0
      do i=1,Np
        xi = xmin+dx*(i-1)
        psi(i) = (0d0,0d0)
        do j=1,nb
          psi(i) = psi(i)+c0(j)*exp(-al/2d0*(xi-x(j))**2+
     +             im*p(j)*(xi-x(j)))*
     +             dsqrt(dsqrt(al/pi))
        enddo
        prob = prob+abs(psi(i))**2*dx
        if(abs(psi(i)) .gt. 1d-3) write(104,10000) xi,abs(psi(i))**2
      enddo
      write(*,1013)  prob
1013  format('Normalization =', f10.6)
      

!-----output at t=0-------------------
      write(107,10000) t,(abs(c(i)),i=1,Ntraj)
      dt = dt/kout
      dt2=dt/2d0
      t=0d0

!-----begin the time step loop----------
      do 10 k=1,kmax
        
        do 13 kk=1,kout
          t = t+dt

          do i=1,Nb
            c(i) = c(i)+dt*dc(i)
          enddo

!         half-step increments of moment, full step increment of positions
          do i=1,Nb
            p(i)=p(i)+(-dv(i)-du(i))*dt2
            x(i)=x(i)+p(i)*dt/am
          enddo

!          call reduce(nb,x,p,c)
!          print *,'time step',k,nb

          call derivs(x,Ntraj,v,dv,ddv)
          call qpot(Nb,ddv,c,am,p,x,w,u,du)
          
          do i=1,Nb
            p(i)=p(i)+(-dv(i)-du(i))*dt2
          enddo
          
!-------------update dc/dt--------------------
          call multi(nb,x,p,mat)
          call ham(am,w,c,nb,x,p,du,v,ddv,h)
          call increm(nb,mat,h,c,dc)

13      enddo


        write(107,10000) t,(abs(c(i)),i=1,Ntraj)

!       print out trajectories in phase space
        do i=1,nb
          write(101,10000) x(i),p(i)
        enddo
        
        write(105,10000) t,(p(i),i=1,Nb)
        write(106,10000) t,(x(i),i=1,Nb)
!       calculate the total energy
        enk = 0d0
        env = 0d0
        enq = 0d0
	  do i=1,Nb
          env = env+v(i)*w(i)
          enk = enk+p(i)**2/2d0/am*w(i)
          enq = enq+u(i)*w(i)
        enddo
        en = enk+env+enq
        write(100,10000) t,enk,env,enq,en
        call flush(100)
!        call aver(Ntraj,x,p,c,mat,xav)
        call norm(t,nb,x,p,c,mat)
!-------------average position --------------
!        call aver(Ntraj,x,p,c,xav)
!        write(102,10000) t,xav
10    end do
!--------print x,c--------------
        do i=1,nb
          write(113,10000) x(i),c(i)
        enddo
!---------check matrix singularity--------------
        do i=1,nb
          write(*,1018) (mat(i,j),j=1,nb)
        enddo
        call eigen(nb,mat)

1018    format(20(f6.4,1x))
!----------check weights--------------
      ww = 0d0
      do i=1,nb
        ww = ww+w(i)
      enddo
      write(*,*) 'Weigths = ', ww

!     print out final wavefunction         
      prob = 0d0
      do i=1,Np
        xi = xmin+dx*(i-1)
        psi(i) = (0d0,0d0)
        do j=1,nb
          psi(i) = psi(i)+c(j)*exp(-al/2d0*(xi-x(j))**2+
     +          im*p(j)*(xi-x(j)))*
     +          dsqrt(dsqrt(al/pi))
        enddo
        prob = prob+abs(psi(i))**2*dx
      enddo

      write(*,*) 'Final prob = ',prob

      do i=1,np
        xi = xmin+dx*(i-1)
        if(abs(psi(i)) .gt. 1d-3) then 
          write(103,10000) xi,abs(psi(i))**2/prob
        endif
      enddo

10000 format(2000(e14.7,1x))
      end program main

!------------------------------------------------

!------normalization-----------------------------
      subroutine norm(t,nb,x,p,c,mat)
      implicit real*8(a-h,o-z)
!------compute normalization--------------------
      real*8 x(nb),p(nb)
      complex*16 mat(nb,nb),c(nb),z

      z = (0d0,0d0)
      do j=1,nb
        do i=1,nb
          z=z+conjg(c(j))*mat(j,i)*c(i)
        enddo
      enddo

      write(111,1000) t,z

1000  format(20(e12.5,1x))
      return
      end subroutine

!     -----------------------------------------
!     trajectory average
!     -------------------------------------------
!     ----------------------------------------------
!     average
!     ---------------------------------------------
      subroutine aver(Ntraj,q,p,c,xav)
      use cdat, only: im,pi

      implicit real*8(a-h,o-z)
      integer*4,intent(in) :: Ntraj
      real*8,intent(in)::q(Ntraj),p(Ntraj)
      complex*16, intent(in) :: c(Ntraj)
      real*8,intent(out) :: xav
      complex*16 :: psi,z
      common/wave/al,q0,p0,a0
      common/grid/xmin,xmax,dx,np
            
      z = (0d0,0d0)
      do i=1,Np
        xi = xmin+dx*(i-1)
        psi = (0d0,0d0)
        do j=1,ntraj
          psi = psi+c(j)*exp(-al/2d0*(xi-q(j))**2+im*p(j)*(xi-q(j)))*
     +          dsqrt(dsqrt(al/pi))
        enddo
        z = z+conjg(psi)*psi*dx*xi
      enddo

      xav = real(z)

      return
      end subroutine
!     -----------------------------------
!     gassian wavepackets overlap matrix
!     -----------------------------------
      subroutine multi(nb,q,p,mat)
      use cdat
      implicit real*8(a-h,o-z)
      integer*4, intent(in) :: nb
      real*8   , intent(in) :: q(nb),p(nb)
      complex*16,intent(out) :: mat(nb,nb)
      complex*16 :: z,z0
      common/wave/al,q0,p0,a0

      do j=1,nb
        do k=1,nb
!          c = (0d0,0d0)
!          do i=1,Np
!            x = xmin + dx*(i-1)
!            gi = dsqrt(dsqrt(2d0*al/pi))*exp(-al*(x-q(m))**2+
!     &           im*p(m)*(x-q(m)))
!            gj = dsqrt(dsqrt(2d0*al/pi))*exp(-al*(x-q(n))**2+
!     &           im*p(n)*(x-q(n)))
!            c = c + conjg(gi)*gj*dx
!          enddo
          dq = q(j)-q(k)
          z = dq/2d0-im/2d0/al*(p(j)-p(k)) 
          z0 = exp(-al*z*conjg(z)+im/2d0*(p(j)+p(k))*dq)
          mat(j,k) = z0
        enddo
      enddo
      end subroutine
!     -------------------------------------------
!     hamiltonian matrix (H-id/dt) with gaussian 
!     basis
!     -------------------------------------------
      subroutine ham(am,w,c,nb,q,p,du,v,ddv,h)
      use cdat, only: pi,im
      implicit real*8(a-h,o-z)
      integer*4, intent(in) :: nb
      complex*16,intent(in) :: c(Nb)
      real*8, intent(in) :: q(nb),p(nb),w(Nb),am
      real*8 :: du(nb),v(Nb),ddv(Nb),fx(Nb)
      complex*16,intent(out) :: h(nb,nb)
      complex*16 :: z,z0,d0,d1,d2
     
      common/wave/al,q0,p0,a0

!     mat = (0d0,0d0)
!     h = (0d0,0d0)

!      do i=1,nb
!        q(i) = x(i)
!      enddo

!      i = 0
!      j = 0
!      nn = nb
!      do while(i .lt. nb) 
!        i = i+1
!        j = j+1
!        dq = q(i+1)-q(i)
!        if(abs(dq) .gt. 1d-3) then
!          cn(j) = c(i)
!        else
!          cn(j) = c(i)+c(i+1)
!          i = i+1
!          nn = nn-1
!        endif
!      enddo

      do j=1,nb
        do k=1,nb
          dq = q(j)-q(k)
          z = dq/2d0-im/2d0/al*(p(j)-p(k))
          z0 = exp(-al*z*conjg(z)+im/2d0*(p(j)+p(k))*dq)

          d0 = z0*(v(k)-(p(k)**2-al)/2d0/am)
          d1 = -du(k)*z*z0
          d2 = 0.5d0*(ddv(k)-al**2/am)*(z**2+1d0/2d0/al)*z0

!          mat(j,k) = z0          
          h(j,k) = d0+d1+d2 
        enddo
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
      do i=1,Nb
        b(i) = (0d0,0d0)
        do j=1,Nb
          b(i) = b(i)+h(i,j)*c(j)
        enddo
      enddo
!-------save the overlap matrix-----------
      do j=1,nb
        do i=1,nb
          aux(i,j) = mat(i,j)
        enddo
      enddo
!      do i=1,Nb
!        write(*,1023) (aux(i,j),j=1,nb)
!      enddo
!      write(*,*) (b(i),i=1,Nb)
!1023  format(20(f6.4,1x))  
      call zposv('U', Nb, 1, aux, nb, b, nb, INFO)
      
      if(INFO .ne. 0) then
        write(*,*) 'dC/dt: Matrix fails, INFO =',INFO
        stop
      endif
!------compute increment dc------------
      do j=1,nb
        dc(j) = -im*b(j)
      enddo

      return
      end subroutine
!     ----------------------------------------------
!     initial coefficients before gaussian basis
!     ----------------------------------------------
      subroutine proj(nb,x,p,mat,c0)
      use cdat, only:im,pi
      implicit real*8(a-h,o-z)
      integer*4,intent(in) :: nb
      real*8,intent(in) :: x(nb),p(Nb)
      complex*16,intent(out) :: c0(Nb)
      complex*16 :: z,z0,mat(Nb,Nb),psi0,aux(nb,nb)
      common/wave/al,q0,p0,a0
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
        

!     solve matrix equation, Mc = b    
!      do i=1,nb
!        write(*,1009) (mat(i,j),j=1,Nb)
!      enddo
!1009  format(20(f5.2,1x))
!------save overlap matrix------------------
      do j=1,nb
        do i=1,nb
          aux(i,j) = mat(i,j)
        enddo
      enddo

      call zposv('U', Nb, 1, aux, nb, c0, nb, INFO)

      IF( INFO.ne.0 ) THEN
         WRITE(*,1222) INFO
1222     format('INFO = ', i6/ ,         
     &          'Projection fails.'/)
!         WRITE(*,*)'definite; the solution could not be computed.'
         STOP
      END IF
      
      return
      end subroutine

!-----potential------------------------
      subroutine derivs(x,Ntraj,v,dv,ddv)
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Ntraj
      real*8,intent(IN) :: x(Ntraj)
      real*8,intent(OUT) :: v(Ntraj),dv(Ntraj),ddv(Ntraj)

! --- eckart barrier 
      v0 = 16d0
      a = 0.8d0

      do i=1,ntraj
        v(i) = v0/cosh(x(i)/a)**2
        dv(i) = -2d0*v0*tanh(x(i)/a)/a/(cosh(x(i)/a))**2
        ddv(i) = 4d0*v0*tanh(x(i)/a)**2/a**2/cosh(x(i)/a)**2-2d0*v0*
     +           (1d0-tanh(x(i)/a)**2)/a**2/cosh(x(i)/a)**2
      enddo

!------harmonic----------------------
!      ak=1.d0
!      do i=1,Ntraj
!        dv(i) = ak*x(i)
!        v(i)  = ak*x(i)**2/2d0
!        ddv(i) = ak
!      enddo
!------morse potential---------------
!      de = 0.176d0
!      x0 = 1.4d0
!      a = 1.02d0
!      do i=1,Ntraj
!        d = (1d0-exp(-a*(x(i)-x0)))
!        v(i) = de*d**2
!        dv(i) = 2*de*d*a*exp(-a*(x(i)-x0))
!        ddv(i)=2*de*(-d*exp(-a*((x(i)-x0)))*a**2+
!     +         (exp(-a*(x(i)-x0)))**2*a**2)       
!      enddo
!---------doulbe well---------------------------
!      eta = 1.3544d0
!      do i=1,Ntraj
!        v(i) = 1d0/16d0/eta*x(i)**4-0.5d0*x(i)**2
!        dv(i) = 1d0/4d0/eta*x(i)**3-x(i)
!        ddv(i) = 3d0/4d0/eta*x(i)**2-1d0
!      enddo

      return
      end subroutine

