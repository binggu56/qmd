!     polynomial fitting for {r}
!     local linear fitting with domains 

! --- J. Chem. Phys. 137, 074115, 2012 

      subroutine fit(am,x,p,w,ndim,ntraj,u,du)

      use cdat
      use domain 

      implicit real*8(a-h, o-z)

!      integer (kind = 4), parameter :: nd = 4, nq = 1 
!      integer (kind = 4), parameter :: nb = ndim+1
      real*8, dimension(ntraj)      :: w
      real*8, dimension(ndim)     :: am
      real*8, dimension(ndim,ntraj) :: x,p,du

      real (kind = 8) du_dom(ndim,ntraj) 

      real (kind = 8), dimension(ntraj) ::dom,ddom,dom1,ddom1,dom2,ddom2 

      real (kind = 8) r(nq,ntraj) 
!      real*8, dimension(ndim,ndim,ntr) :: ddp,ddr,dp,dr

      real*8 f(ndim+1),s(ndim+1,ndim+1),cr(ndim+1,nq),mat(ndim+1,ndim+1)
      integer info

      real (kind = 8), dimension(nd+1) :: rt  ! trajectories ratio in each domain 
 
!---------matrix S1 & c for linear fitting------------------
      nb = ndim+1

      cut = 0.02d0

      lc = 0d0  

      rt = 0d0 
      do l=1,nd 
        do i=1,ntraj 
          rt(l) = rt(l) + step(ydom(l),ydom(l+1),x(id,i))
        enddo 
      enddo

      rt = rt/dble(ntraj)  
      rt(nd+1) = 1d0 

      open(110, file = 'ratio.dat', access = 'append') 

      write(110,1000) (rt(l),l=1,nd)
!1001  format('RATIO FOR EACH DOMAIN',20(f9.4,1x)/)

!      err = 0d0
!      call basis(Ndim,Ntr,x,f)
      du = 0d0 
      u = 0d0 
     
      dom1 = 0d0 
      ddom1 = 0d0 
 
      l = 1 
      do while(l < nd+1) ! loop over domains    

        if( rt(l) > cut .and. rt(l+1) > cut) then 
         
        lc = lc+1 


! --- count trajectories in l-th domain 
!        nn = 0 
!        do i=1,ntraj 
!          if(x(id,i) > ydom(l) .and. x(id,i) < ydom(l+1)) 
!             nn = nn+1 
!          endif 
!        enddo 
        
!      if(nn > ncut) then 
          call domfunc(l,ndim,ntraj,x,dom,ddom) 

          dom = dom + dom1 
          ddom = ddom + ddom1 

        cr = 0d0 
      
      d1 = 0d0 
      do i=1,ntraj 
        d1 = d1 + w(i)*dom(i) 
      enddo 
!      do i = 1,ntraj
        
!        do j=1,ndim 
!          f(j) = x(j,i) 
!        enddo 
!        f(ndim+1) = 1d0 
        
!        do j = 1,ndim+1
!        do k = 1,nq 
!            cp(j,k) = cp(j,k)+f(j,i)*p(k,i)*w(i) 
!            cr(j,k) = cr(j,k)+f(j)*r(k,i)*w(i)
!        enddo
!        enddo
!      enddo 
      do k=1,nq
          cr(k,k) = -0.5d0*d1  
      enddo 

! --- overlap matrix  <f*f> computed with domain function 
      mat = 0d0 

      do i = 1,ntraj        
        do j=1,ndim 
          f(j) = x(j,i) 
        enddo 
        f(ndim+1) = 1d0 

        do k1=1,ndim+1 
        do k2 = k1,ndim+1 
            mat(k1,k2) = mat(k1,k2)+f(k1)*f(k2)*w(i)*dom(i) 
        enddo
        enddo
      enddo

      do j = 2,nb 
        do k = 1,j-1 
          mat(j,k) = mat(k,j)
        enddo 
      enddo  
          
!---------------solve matrix M*c = f*p-------------
!      mat = s1
!      call dposv('U',nb,ndim,mat,nb,cp,nb,INFO)
!      if(info/=0) then
!        write(*,*) 'linear fitting of p failed.'
!        stop
!      endif
!-------------same for r--------------------------
      s = mat 

!      print *,s(1,1),s(1,3),s(2,2),cr(1,1),cr(2,1)
      call dposv('U',nb,nq,mat,nb,cr,nb,INFO)
      if(info/=0) then
        write(*,*) 'linear fitting of r failed for',l,'th domain.'
        stop
      endif

!-------------output for {ap,ar}--------------------
!      dp = 0d0
!      ap = 0d0
!      ar = 0d0
!      dr = 0d0
!      ddp = 0.d0
!      ddr = 0.d0
!
!      do i=1,ntr
!        do j=1,ndim
!          do k=1,nb
!            ap(j,i) = ap(j,i)+cp(k,j)*f(k,i)
!            ar(j,i) = ar(j,i)+cr(k,j)*f(k,i)
!          enddo
!        enddo
!      enddo
!!------first-order derivative of p,r, second order is 0 for linear fitting---------------
!      do i=1,ntr
!        do j=1,ndim
!          do k=1,ndim
!            dp(k,j,i) = cp(k,j)
!            dr(k,j,i) = cr(k,j)
!          enddo
!        enddo
!      enddo

      r = 0d0 

      do i=1,ntraj 
        do j=1,ndim
          f(j) = x(j,i)
        enddo
        f(ndim+1) = 1d0
     
        do j=1,nq  
          do k=1,ndim+1  
            r(j,i) = r(j,i) + f(k)*cr(k,j) 
          enddo 
        enddo 

      enddo 
! --- quantum potential 
      udom = 0d0 
      do m=1,nq 

      d = 0d0 
      do k=1,ndim+1 
        do j=1,ndim+1 
          d = d + s(k,j)*cr(k,m)*cr(j,m) 
        enddo 
      enddo 


      udom = udom -half/am(m)*(d + cr(m,m)*d1)  ! quantum potential average on one domain 
      enddo 

! --- quantum force (quantum DoF and classical DoF separately) 

      du_dom = 0d0 
      do i=1,ntraj 
        do j=1,ndim 
          do k=1,nq 
            du_dom(j,i) = du_dom(j,i)-r(k,i)*cr(j,k)/am(k)  
          enddo 
          du_dom(j,i) = du_dom(j,i)*dom(i) 
        enddo 
      enddo 

      do i=1,ntraj
        d0 = 0d0 
        do k=1,nq 
          d0 = d0 + (-(r(k,i)**2 + cr(k,k))/2d0/am(k))  
        enddo 
       
        du_dom(id,i) = du_dom(id,i)+ d0*ddom(i) 
       
      enddo 
              
      
      du = du + du_dom 
      u = u + udom 

      l = l+1 
      dom1 = 0d0 
      ddom1 = 0d0 
      
      else 
 
        call domfunc(l,ndim,ntraj,x,dom2,ddom2) ! combine the domain l and domain l+1 
!        call domfunc(l+1,ndim,ntraj,x,dom,ddom)

        dom1 = dom1 + dom2 
        ddom1 = ddom1 + ddom2  

        rt(l+1) = rt(l) + rt(l+1) 
        l = l+1  

      endif 

      end do  ! domain loop 

!      print *,'qpot',u,du(1,3),du(2,1)

!      do i=1,ntr
!        do j=1,ndim
!          u(i) = u(i)-(ar(j,i)**2+dr(j,j,i))/2d0/am(j)
!        enddo
!        
!        do j=1,ndim
!
!          do k=1,ndim
!            du(j,i) = du(j,i)-(ar(k,i)*dr(j,k,i))/am(k)
!            fr(j,i) = fr(j,i)-dp(j,k,i)*ar(k,i)/am(k)
!          enddo
!      
!          du(j,i) = du(j,i)-ddr(j,j,i)/2d0/am(j)
!          fr(j,i) = fr(j,i)-ddp(j,j,i)/2d0/am(j)
!        enddo
!
!      enddo

      write(6,*) 'computed number of domains',lc 
1000  format(20(e14.7,1x)) 
      return
      end subroutine
!----------------------------------------
!     linear basis 
!---------------------------------------
!      subroutine basis(Ndim,Ntr,x,f)
!      implicit real*8(a-h,o-z)
!      integer*4  :: Ndim,Ntr
!      real*8      :: x(Ndim,Ntr)
!      real*8,intent(OUT) :: f(Ndim+1,Ntr)
!
!      f = 0d0
!
!!---basis vector f = ((x(i),i=1,Ndim),1) for each trajectory-------------------
!      do i=1,Ntr
!        do j=1,Ndim
!          f(j,i)=x(j,i)
!        enddo
!          f(Ndim+1,i)=1d0
!       enddo
!     
!      return
!      end subroutine
!--------------------------------------------------------------
! --- domain function depending on the domain index l 

      subroutine domfunc(l,ndim,ntraj,x,dom,ddom) 
      use cdat
      use domain  

      implicit double precision (a-h,o-z) 

      real (kind = 8) x(ndim,ntraj)
      real (kind = 8) y(ntraj),dom(ntraj),ddom(ntraj) 
      
      ak = 5d0 

      do i=1,ntraj 
        y(i) = x(id,i)      
      enddo 
      
      do i=1,ntraj 
        dom(i) = 0.5d0*(tanh(ak*(y(i)-ydom(l)))-tanh(ak*(y(i)-ydom(l+1))))
        ddom(i) = ak*half*(-tanh(ak*(y(i)-ydom(l)))**2 + tanh(ak*(y(i)-ydom(l+1)))**2)
      enddo     

      return 
      end subroutine

