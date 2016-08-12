      subroutine fitp(am,x,px,rx,w,Ntraj,dpx,ddp,dr,ddr,u,ap,ar,err)
      
      implicit real*8(a-h,o-z)

      integer*4,intent(IN) :: Ntraj
      integer*4,parameter :: nb=3
      real*8 :: px(Ntraj), w(Ntraj),x(Ntraj),rx(ntraj)
      real*8 :: dpx(Ntraj),dr(ntraj),ddr(ntraj)
      real*8 :: f(nb),s(nb,nb),c(nb,2),cr(nb,1),mat(nb,nb),u(ntraj), ap(ntraj),ar(ntraj),ddp(ntraj)
      real*8 :: err(2) 

      c = 0d0
      s = 0d0
      cr = 0d0

   do j = 1, nb 
      do i=1,Ntraj

!        f=(/x(i),x(i)**2,x(i)**3,1d0/)
        c(j,1) = c(j,1) + px(i) * w(i) * x(i)**(j-1) 
        c(j,2) = c(j,2) + rx(i) * w(i) * x(i)**(j-1) 
!        c(2,1) = c(2,1)+px(i)*f(2)*w(i)
!        c(3,1) = c(3,1)+f(3)*px(i)*w(i)
!        c(4,1) = c(4,1)+f(4)*px(i)*w(i)
      enddo
	enddo 

!      do i=1,Ntraj
!		do j=1,nb 
!        f=(/x(i),x(i)**2,x(i)**3,1d0/)
!        cr(j,1) = cr(j,1)+rx(i) * w(i) * x(i)**(j-1) 
!        cr(2,1) = cr(2,1)+rx(i)*f(2)*w(i)
!        cr(3,1) = cr(3,1)+f(3)*rx(i)*w(i)
!        cr(4,1) = cr(4,1)+f(4)*rx(i)*w(i)
!      enddo!
!		enddo 
    
      do m=1,nb
        do n=1,nb
          do i=1,Ntraj
            !f=(/x(i),x(i)**2,x(i)**3,1d0/)
            s(m,n) = w(i)* x(i)**(m+n-2) + s(m,n)            
          end do
        end do
      end do
      

! calculate matrix c(t)
        mat = s
        call DPOSV('U',nb,2,mat,nb,c,nb,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, 'S = ', s
        print *, 'c = ', c 
        print *, "polynomial fitting of p : matrix fails" 
        stop
        end if

!        call DPOSV('U',nb,1,s,nb,cr,nb,INFO)
!        if(INFO/=0) then
!        print *, "info=",info
!        print *, "polynomial fitting of r : matrix fails"
!        stop
!        end if

!     check the error of polyfit

! the momentum operator r=cf
	ap = 0d0 
	ar = 0d0 
	dpx = 0d0 
	dr = 0d0 
	ddp = 0d0 
	ddr = 0d0 
	
      do i=1,Ntraj
		do j=1,nb 
        ap(i) = ap(i) + c(j,1)*x(i)**(j-1) 
        ar(i) = ar(i) + c(j,2)*x(i)**(j-1) !+cr(2,1)*x(i)**2+cr(3,1)*x(i)**3+cr(4,1)
        enddo 
        
        do j=2,nb 
        dpx(i) = dpx(i) + c(j,1) * (j-1) * x(i)**(j-2) !+2d0*c(2,1)*x(i)+c(3,1)*3d0*x(i)**2
        dr(i) = dr(i) + c(j,2) * (j-1) * x(i)**(j-2) ! + 2d0*x(i)*cr(2,1) + cr(3,1)*3d0*x(i)**2
        enddo 
        
        do j=3,nb 
        ddp(i) = ddp(i) + c(j,1) * (j-1)*(j-2) * x(i)**(j-3) !2d0*c(2,1) + 6d0*x(i)*c(3,1)
        ddr(i) = ddr(i) + c(j,2) * (j-1)*(j-2) * x(i)**(j-3) !2d0*cr(2,1) + 6d0*x(i)*cr(3,1)
		enddo 
	enddo 
	
	do i=1,ntraj 
        u(i) = -(ar(i)**2 + dr(i))/(2d0*am)
    enddo
      
      err = 0d0  
      do i=1,ntraj 
        err(1) = err(1) + (px(i)-ap(i))**2*w(i) 
        err(2) = err(2) + (rx(i)-ar(i))**2*w(i) 
      enddo
!c     &             py(i)*(c(1,2)+2d0*c(3,2)*x(i)+c(5,2)*y(i))/am2
!c          dpy(i) = px(i)*(c(2,1)+2d0*c(4,1)*y(i)+c(5,1)*x(i))/am1+
!c     &             py(i)*(c(2,2)+2d0*c(4,2)*y(i)+c(5,2)*x(i))/am2
!! calculate quantum potential
!c          qp(i) = (c(1,1)+2d0*c(2,1)*x(i)+c(3,1)*3d0*x(i)**2)/2d0/am1
!c     &            (c(2,2)+2d0*c(4,2)*y(i)+c(5,2)*x(i))/2d0/am2
!c          qfx(i) = -(2d0*c(2,1)+6d0*x(i)*c(3,1))/2d0/am1
!c     qfy(i) = -(c(5,1)/2d0/am1+c(4,2)/am2)

      return
      end subroutine

!     --------------------------------
!     check the error of poly-fit of p
!     --------------------------------
      subroutine errp(nt,p,x,w,c,error)
      implicit real*8(a-h,o-z)

      integer*4, intent(in)  :: nt
      real*8,    intent(in)  :: w(nt),p(nt),x(nt),c(4,1)
      real*8,    intent(out) :: error 
      
      error = 0d0

      do i=1,nt
        pfit = c(1,1)*x(i)+c(2,1)*x(i)**2+c(3,1)*x(i)**3+c(4,1)
        error = error + (p(i)-pfit)**2*w(i)
      enddo
      

      return
      end subroutine


