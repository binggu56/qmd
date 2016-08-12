      subroutine fitp(x,px,w,Ntraj,dpx,err)
      implicit real*8(a-h,o-z)

      integer*4,intent(IN) :: Ntraj
      real*8, intent(IN)   :: px(Ntraj), w(Ntraj),x(Ntraj)
      real*8, intent(OUT)  :: dpx(Ntraj),err
      real*8 :: f(4),s(4,4),c(4,1)

      c = 0d0
      s = 0d0
      nb = 4

      do i=1,Ntraj
        f=(/x(i),x(i)**2,x(i)**3,1d0/)
        c(1,1) = c(1,1)+px(i)*f(1)*w(i)
        c(2,1) = c(2,1)+px(i)*f(2)*w(i)
        c(3,1) = c(3,1)+f(3)*px(i)*w(i)
        c(4,1) = c(4,1)+f(4)*px(i)*w(i)
      enddo



      do m=1,nb
        do n=1,nb
          do i=1,Ntraj
            f=(/x(i),x(i)**2,x(i)**3,1d0/)
            s(m,n)=w(i)*f(m)*f(n)+s(m,n)
          end do
        end do
      end do

! calculate matrix c(t)
        call DPOSV('U',nb,1,s,nb,c,nb,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, "polynomial fitting : matrix fails"
        stop
        end if
c     check the error of poly-fit
      call errp(ntraj,px,x,c,err)  

! the momentum operator r=cf
        do i=1,Ntraj
          dpx(i) = c(1,1)+2d0*c(2,1)*x(i)+c(3,1)*3d0*x(i)**2
c     &             py(i)*(c(1,2)+2d0*c(3,2)*x(i)+c(5,2)*y(i))/am2
c          dpy(i) = px(i)*(c(2,1)+2d0*c(4,1)*y(i)+c(5,1)*x(i))/am1+
c     &             py(i)*(c(2,2)+2d0*c(4,2)*y(i)+c(5,2)*x(i))/am2
! calculate quantum potential
c          qp(i) = (c(1,1)+2d0*c(2,1)*x(i)+c(3,1)*3d0*x(i)**2)/2d0/am1
c     &            (c(2,2)+2d0*c(4,2)*y(i)+c(5,2)*x(i))/2d0/am2
c          qfx(i) = -(2d0*c(2,1)+6d0*x(i)*c(3,1))/2d0/am1
c     qfy(i) = -(c(5,1)/2d0/am1+c(4,2)/am2)
        enddo


      return
      end subroutine

c     --------------------------------
c     check the error of poly-fit of p
c     --------------------------------
      subroutine errp(nt,p,x,c,err)
      implicit real*8(a-h,o-z)

      integer*4, intent(in)  :: nt
      real*8,    intent(in)  :: p(nt),x(nt),c(4,1)
      real*8,    intent(out) :: err
      
      common/smaple/xmin,xmax,dx
     
      err = 0d0
      do i=1,nt
        pfit = c(1,1)*x(i)+c(2,1)*x(i)**2+c(3,1)*x(i)**3+c(4,1)
        err = err+(p(i)-pfit)**2*dx
      enddo

      return
      end subroutine


