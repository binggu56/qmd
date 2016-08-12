	subroutine xex(ni,nj,psi,xav,yav,xx,yy,ex)
	implicit none
	integer*4 i,j,ni,nj
	real*8 :: xmin,ymin,xmax,ymax,dx,dy,x,y,wx,wy 
	complex*16 :: psi(ni,nj)
	real*8,intent(out) :: xav,yav,xx,yy,ex
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
	common /para1/ wx,wy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	xav=0d0
	yav=0d0
        xx=0d0
        yy=0d0
	ex=0d0
c	evx=0d0
c	evy=0d0

c 	calculate the average value of x,y 
	do j=1,nj
	y=ymin+(j-1)*dy
	do i=1,ni
	x=xmin+(i-1)*dx
	xav=xav+abs(psi(i,j))**2*dx*dy*x
        xx=xx+abs(psi(i,j))**2*dx*dy*x**2
	yav=yav+abs(psi(i,j))**2*dx*dy*y
        yy=yy+abs(psi(i,j))**2*dx*dy*y**2
c	evx=evx+abs((psi(i,j))**2*dx*dy*wx**2*x**2/2d0
c	evy=evy+abs((psi(i,j))**2*dx*dy*wy**2*y**2/2d0
	ex=ex+abs(psi(i,j))**2*dx*dy*(x**2*y-y**3/3d0)
c        if ( x.lt.-5d0) then
c        prob=abs(psi(i,j))**2*dx*dy+prob
c        endif
	enddo
	enddo

	end subroutine
	
	
