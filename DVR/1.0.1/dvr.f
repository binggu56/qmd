        program main
!========================================================================
!------------------ Miller-Colbert DVR (grid DVR) in one dimension ------
!------------------ equidistant grid, arbitrary number of points --------
!------------------ gives eigenvalues and eigenvectors of H -------------
!------------------ symmetric double well is defined --------------------
!========================================================================
	implicit real*8(a-h,o-z)
	parameter(nx=4096,aw=14d0,bw=20d0,vbot=bw*bw/4d0/aw)
      REAL*8 x(nx),v(nx)
!========================================================================
      OPEN(15,FILE='INdvr',STATUS='OLD')
      READ(15,*) n1,xmin,xmax
      read(15,*) am
      close(15)

      write(*,*) 'GRID',n1,xmin,xmax
      write(*,*) 'MASS=',am
      pi=4d0*atan(1d0)

	open(55, file='wf_dat')
	open(54, file='pot'   )
	open(62, file='eigs'  )
      open(101,file='r'     )
      open(11,file='pes.dat')
!=======================================================================
!cccccccccccccccccccccccccc define the grid cccccccccccccccccccccccccccccc
      xn=xmax-xmin
      dx=xn/(n1-1)
      write(*,*) 'GRID SPACING=',dx
	write(*,*) 'x range',xmin,xmax
cccccccccccccccccccccccccc define the potential ccccccccccccccccccccccccc
	De = 0.176d0
      epi = 0.5d0

c      a0 = -2.22384d-2
      a0 = 0d0
      a1 = 6.49557d-6
      a2 = -4.50827d-6
      a3 = 1.56397d-6
      a4 = 3.61582d-8
c      a12 = 6.5047d-6
      aww = 0.25d0
      bww = 2.d0
      alfa = 0.05d0
      ipot = 0 

      if(ipot == 0) then 
        do 11 i=1,N1
         x(i)=xmin+dx*(i-1)
!	   if(x(i) .lt. 0.653594d0) then
!	     v(i)=0.1831104*x(i)**2-0.186772608d0*x(i)**3
!	   else
!	     v(i) = 0.026074d0
!	   endif
!	  v(i) = De*(1d0-exp(-1.02d0*x(i)))**2
C	  v(i) = x(i)**2/2d0-0.05d0*x(i)**3
c         read(11,*) a,v(i)
         v(i) = (1d0-2d0*x(i)**2)**2/8d0 
c        v(i) = a0+a1*x(i)**2+a2*x(i)**4+a3*x(i)**6+a4*x(i)**
        write(54,*) x(i),v(i)
11      continue
        close(54)
      endif 

!=======================================================================
!------------------- construct and diagonalize H -----------------------
        call ham_dvr(n1,pi,am,x,v)

        close(62)
        close(55)
        STOP
        END
!========================================================================
!========================================================================
        subroutine ham_dvr(n1,pi,am,x,v)
!------ construct the DVR Hamiltonian and solve the eigenvalue problem --
!========================================================================
        implicit real*8(a-h,o-z)
        parameter(nx=4096,lwork=nx*10)
        real*8 x(*),v(*),work(lwork),w(nx),ham(nx,nx),hsv(nx,nx),z(nx)
        real*8 r(n1),dr(n1),ddr(n1)
        a=x(1)
        b=x(n1)
        ba=pi*a/(b-a)
        dx=x(2)-x(1)
        an=0.5d0/am/dx/dx
        anf=pi/(b-a)**2/am
c----------------- construct the hamiltonian matrix --------------
        do 10 j=1,n1
         ham(j,j)=an*pi*pi/3d0+v(j)
         do 9 i=j+1,n1
          ham(i,j)=an*(-1d0)**(i-j)*2d0/(i-j)**2
          hsv(i,j)=ham(i,j)
          hsv(j,i)=ham(i,j)
9       continue
10      hsv(j,j)=ham(j,j)
        write(*,*) ham(1,1),ham(1,2),ham(2,1),ham(2,2)
        call dsyev('V','L',n1,ham,nx,w,work,lwork,info)
        write(*,*) 'INFO',info,work(1)
        eig=w(1)
        write(*,*) 'EIGENVALS',w(1),w(2),w(3)
!---------------------- check the eigenvector --------------------
!------ ham(j,k)/dsqrt(dx) is the eigenfunction in coord space ---
        dev=0d0
        write(*,*) 'NUMBER OF EIGENFUNCTIONS TO PRINT'
        read(*,*) ns
        do 30 j=1,n1
         write(62,*) j,w(j)
         write(55,1000) x(j),(ham(j,k)/dsqrt(dx),k=1,ns)

         z(j)=0d0
         do 31 i=1,n1
31        z(j)=z(j)+hsv(j,i)*ham(i,1)/eig
          dev=dev+abs(z(j)-ham(j,1))**2
30      continue

c     print out quantum potential 
        do j=1,n1-1
          r(j) = (abs(ham(j+1,1))-abs(ham(j,1)))
     +                  /dx/abs(ham(j,1))
        enddo

        do j=1,n1-2
          dr(j) = (r(j+1)-r(j))/dx
        enddo

        do j=1,n1-3
          ddr(j) = (dr(j+1)-dr(j))/dx
          write(101,1000) x(j),r(j),dr(j),ddr(j)
        enddo
        

        write(*,*) 'EIGENVECTOR TEST',dev
c--------------------  construct the eigenvector ------------
         x2 = 0d0
         x33 = 0d0
         x22 = 0d0
         do i=1,n1
           x2 = x2+x(i)**2*abs(ham(i,1))**2
           x33 = x33+x(i)**2*abs(ham(i,3))**2
           x22 = x22+x(i)**2*abs(ham(i,2))**2
         enddo
         write(*,*) 'Expectation of x^2 of ground state',x2
         write(*,*) 'Expectation of x^2 of 1st excited state',x22
         write(*,*) 'Expectation of x^2 of 2nd excited state',x33

        anv=0d0
        do 15 n=1,n1
15         anv=anv+abs(ham(n,1))**2
        write(*,*) 'VECTOR NORM=',anv
1000	  format(81(e14.7,1x))
        return
        end
!========================================================================
