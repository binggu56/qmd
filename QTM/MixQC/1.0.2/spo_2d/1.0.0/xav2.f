	subroutine xex2(ni,nj,psi,evx,evy)
	implicit none
	integer*4 i,j,ni,nj
	real*8 :: xmin,ymin,xmax,ymax,dx,dy,x,y,wx,wy 
c	integer*4,parameter :: mn=256
	complex*16 :: psi(ni,nj)
	real*8,intent(out) ::evx,evy
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
	common /para1/ wx,wy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	xav=0d0
c	yav=0d0
c        xx=0d0
c        yy=0d0
c        prob=0d0
	evx=0d0
	evy=0d0

c 	calculate the average value of x,y 
	do j=1,nj
	y=ymin+(j-1)*dy
	do i=1,ni
	x=xmin+(i-1)*dx
c	xav=xav+abs(psi(i,j))**2*dx*dy*x
c        xx=xx+abs(psi(i,j))**2*dx*dy*x**2
c	yav=yav+abs(psi(i,j))**2*dx*dy*y
c        yy=yy+abs(psi(i,j))**2*dx*dy*y**2
	evx=evx+abs(psi(i,j))**2*dx*dy*(wx**2)*(x**2)/2d0
	evy=evy+abs(psi(i,j))**2*dx*dy*(wy**2)*(y**2)/2d0
c        if ( x.lt.-5d0) then
c        prob=abs(psi(i,j))**2*dx*dy+prob
c        endif
	enddo
	enddo

	end subroutine
	
	
