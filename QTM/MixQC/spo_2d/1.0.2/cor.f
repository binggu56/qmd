	subroutine cor(ni,nj,psi0,psi,auto,au_x,au_y)
	implicit none
	integer*4 i,j,ni,nj
	real*8 :: xmin,ymin,xmax,ymax,dx,dy,x,y,wx,wy 
	complex*16,intent(IN) :: psi0(ni,nj),psi(ni,nj)
	complex*16 :: im
	real*8,intent(out) :: auto
	real*8 :: qx0,qy0,ax,ay,au_x,au_y,px0,py0
      real*8 cor_d

	common /grid/ xmin,ymin,xmax,ymax,dx,dy
	common /para1/ wx,wy
      common /ini/ qx0,qy0,px0,py0
      common /wav/ ax,ay
      common /correlation/ cor_d
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	auto=0d0
	au_x=0d0
	au_y=0d0
	im=(0d0,1d0)
      cor_d=0d0
c------correlation function-------------------------------
        do j=1,nj
                y=ymin+(j-1)*dy
        do i=1,ni
                x=xmin+(i-1)*dx
      auto=auto+conjg(psi(i,j))*psi0(i,j)*dx*dy
	au_x=au_x+conjg(psi(i,j))*exp(-ax*x**2+im*px0*(x-qx0))*dx*dy
      au_y=au_y+conjg(psi(i,j))*exp(-ay*y**2+im*py0*(x-qy0))*dx*dy
      cor_d=cor_d+abs(psi0(i,j))**2*abs(psi(i,j))**2*dx*dy
      enddo
	enddo
c       yy=yy+abs(psi(i,j))**2*dx*dy*y**2
c	ex=ex+abs(psi(i,j))**2*dx*dy*(x**2*y-y**3/3d0)

	end subroutine
	
	
