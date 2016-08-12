	subroutine xsqure(ni,nj,psi,xx,yy)
	implicit none
	integer*4 i,j,ni,nj
	real*8 :: xmin,ymin,xmax,ymax,dx,dy,x,y 
	integer*4,parameter :: mn=256
	complex*16 :: psi(mn,mn)
	real*8,intent(out) :: xx,yy
	common /grid/ xmin,ymin,xmax,ymax,dx,dy
	xx=0d0
	yy=0d0

c 	calculate the average value of x,y 
	do j=1,nj
	y=ymin+(j-1)*dy
	do i=1,ni
	x=xmin+(i-1)*dx
	xx=xx+abs(psi(i,j))**2*dx*dy*x**2
	yy=yy+abs(psi(i,j))**2*dx*dy*y**2
	enddo
	enddo

	end subroutine
	
	
