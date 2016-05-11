      subroutine qpot(order,Nb,ddv,c,am,p,x,s,w,u,du)

! --- local linear quantum force, fit ||r - c*f|| under g_k

      use cdat, only: pi,im,al
      implicit real*8(a-h,o-z)

      integer*4,intent(IN) :: nb,order
      integer*4, parameter :: nfit=2 
      complex*16,intent(in) :: c(nb) 
      real*8,intent(in),dimension(nb) :: x,p,w,s
      real*8, intent(OUT) :: du(Nb),u(Nb)
      integer INFO
      real*8 :: u2(nb),du2(nb)
      real*8, dimension(order) :: xq,wq 
      real*8 :: f(2),ddv(nb),mat(2,2),cr(2,1),b(nfit),cf(nfit),ncm 

      common/grid/xmin,xmax,dx,np 

! --- approximate quantum force for each trajectory, s=fXf, b=r*f
      
      iqpot = 1 

      if(iqpot == 0) then 

! --- parameters for Gauss_Hermite 
      alpha = 0d0 
      beta = al

      traj: do k=1,nb 
        a = x(k) 

!  write ( *, '(a)' ) '    Integral ( -oo < x < +oo ) |x-a|^ALPHA ' // &
!    'exp( - b * ( x - a )^2 ) f(x) dx'        
        call hermite(a,alpha,beta,order,xq,wq)    

        mat = 0d0 

      do j=1,nfit 
        do i=1,j 
          do m=1,order 
            f = (/xq(m),1d0/)
!            gk2 = dsqrt(al/pi)*exp(-al*(x(m)-x(k))**2) 
            mat(i,j) = mat(i,j)+f(i)*f(j)*wq(m) 
          enddo 
        enddo 
      enddo 

      do j=2,nfit 
        do i=1,j-1 
          mat(j,i) = mat(i,j) 
        enddo 
      enddo 

        b = 0d0 
!        do i=1,nfit
! ----- the integration over grids can be replaced by Gaussian-Hermite
!       quadrature rules
        do j=1,order
          cc = xq(j)
          arx = ncm(cc,nb,c,x,p,s) ! r value at xx 
!          gk2 = dsqrt(al/pi)*exp(-al*(x(j)-x(k))**2) 
          
          b(1) = b(1)+arx*xq(j)*wq(j)
          b(2) = b(2)+arx*wq(j)
!          if(k == nb) print *,'r=',cc,arx 
        enddo 

! ----- compute matrix c(t)
!       print *,'b=',b(1),b(2),arx 

        call DPOSV('U',2,1,mat,2,b,2,INFO)
        if(INFO/=0) then
          write(*,*) 'info =',info, & 
                     'QPOT matrix fails for',k,'th trajectory'
          stop
        endif



!      cf(1) = (b(1)*mat(2,2)-b(2)*mat(1,2))/(mat(2,2)-mat(1,2))
!      cf(2) = (b(2)-b(1)*mat(1,2))/(mat(2,2)-mat(1,2))

        ar = b(1)*x(k)+b(2)
        dr = b(1)
        u(k) = -(ar**2+dr)/2d0/am 
        du(k) = -ar*dr/am 

      enddo traj 

      elseif(iqpot == 1) then     ! quantum force from trajectories    

      x1 = 0d0
      x2 = 0d0
      an = 0d0 
      do i=1,nb
        an=an+w(i)
        x1=x1+w(i)*x(i)
        x2=x2+w(i)*x(i)**2
      enddo

      b1=-0.5d0*an*an/(an*x2-x1*x1)
      b2 = 0.5d0*x1*an/(an*x2-x1*x1)
      do i=1,Nb
        r=b1*x(i)+b2
        u(i) = -(r*r+b1)/2d0/am
        du(i) = -b1*r/am
      enddo

      endif 

! ---- define c matrix 
!      do i=1,ntraj 
!        f = (/x(i),1d0/)
!        c(i) = r
!
!! quantum force
!        s=0d0
!        do i=1,Ntraj
!          f=(/x(i),1d0/)
!        do m=1,2
!          do n=1,2
!            s(m,n)=w(i)*f(m)*f(n)+s(m,n)
!          end do
!        end do
!        end do
!

!      call rep(nb,w,x,p,du2)
      
!      do i=1,nb
!        du(i) = du(i)+du2(i)
!      enddo
!
      return
      end subroutine

!------------------------------------------------------
!------extra repulsion added like quantum potential
!----------------------------------------------------
      subroutine rep(nb,w,x,p,du2)
      use cdat
      implicit real*8(a-h,o-z)
      real*8 x(nb),p(nb),dxp(nb),r0(nb,nb),w(nb),du2(nb)
      
      d = 0.16d0 

      r = x(2)-x(1)
      du2(1) = repul(r)

      r = x(nb)-x(nb-1)
      du2(nb) = -repul(r)

      do i=2,nb-1
        r1 = x(i+1)-x(i)
        r2 = x(i)-x(i-1)
        du2(i) = -repul(r2)+repul(r1) 
      enddo
      
      return 
      end subroutine 

! --- repulsive force function over distance 

      double precision function repul(r)

      implicit real*8(a-h,o-z)

      real*8, intent(in) :: r 
      common/wave/al,q0,p0,a0,sigma
      
      r0 = sigma
      d = 1.6d-4
      beta = 2d0
      x = abs(r)

      if(r > r0) then 
        repul = 0d0 
      else
!        v = 1d0 - 1d0/(1d0+exp(-2d0*beta*(r-r0)))
!        v = v*2d0**4
!       repul = d/(exp(beta*x)-1.d0)
        repul = d/x

      endif 

      return 
      end function 

!      c = 1d6
!      d = 1d3
!
!      r = x(2)-x(1)
!      du2(1) = d*exp(-c*r**2)
!
!      r = x(nb)-x(nb-1)
!      du2(nb) = -d*exp(-c*r**2)
!
!      do i=2,nb-1
!          r1 = x(i+1)-x(i)
!          r2 = x(i)-x(i-1)
!          du2(i) = -d*exp(-c*r2**2)+d*exp(-c*r1**2)
!      enddo

!      a = 0.0004d0
!      d = 16d1
!
!      r = x(2)-x(1)
!      du2(1) = d/cosh(r/a)**2
!
!      r = x(nb)-x(nb-1)
!      du2(nb) = -d/cosh(r/a)**2
!
!      do i=2,nb-1
!          r1 = x(i+1)-x(i)
!          r2 = x(i)-x(i-1)
!          du2(i) = -d/cosh(r2/a)**2+d/cosh(r1/a)**2
!      enddo

