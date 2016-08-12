c ----------------------------------------------------------------------

c     this subroutine computes the local energy and potential energy
c     of a configuration.

c ----------------------------------------------------------------------

      subroutine local(q, vloc, dv)

      implicit double precision (a-h, o-z)

      include 'sizes.h'

      include 'qsats.h'
      
      real*8, intent(out) :: vloc

      common /bincom/ bin, binvrs, r2min

c --- alpha is the exponential parameter in psi:

c     psi = N * exp(-alpha*(r-r0)**2) * Jastrow

c --- bb is the exponential parameter in Jastrow:

c     ln Jastrow(ij) = -0.5 * (bb/rij)**5

      dimension q(NATOM3), dlng(NATOM3), d2lng(NATOM3),dv(NATOM3)

!      do i=1, NATOM3
!         dlng(i)=0.0d0
!         d2lng(i)=0.0d0
!      end do
!
!      do i=1, NATOMS
!
!         xx=q(3*i-2)
!         yy=q(3*i-1)
!         zz=q(3*i)
!
!         dlng(3*i-2)=dlng(3*i-2)-2.0d0*aa*xx
!         dlng(3*i-1)=dlng(3*i-1)-2.0d0*aa*yy
!         dlng(3*i)  =dlng(3*i)  -2.0d0*aa*zz
!
!         d2lng(3*i-2)=d2lng(3*i-2)-2.0d0*aa
!         d2lng(3*i-1)=d2lng(3*i-1)-2.0d0*aa
!         d2lng(3*i)  =d2lng(3*i)  -2.0d0*aa
!
!      end do

! --- loop over all interacting pairs.
! --- potential for one configuration

      vloc=0.0d0

      do n=1, nvpair2

         i=ivpair2(1, n)
         j=ivpair2(2, n)

         dx=-((q(3*j-2))+vpvec2(1, n)+(-q(3*i-2)))
         dy=-((q(3*j-1))+vpvec2(2, n)+(-q(3*i-1)))
         dz=-((q(3*j))  +vpvec2(3, n)+(-q(3*i))  )

         r2=dx*dx+dy*dy+dz*dz
         r = dsqrt(r2)

         ibin=int((r2-r2min)*binvrs)+1

         if (ibin.gt.0) then
            dr=(r2-r2min)-bin*dble(ibin-1)
            vloc=vloc+v(1, ibin)+v(2, ibin)*dr
         else
            vloc=vloc+v(1, 1)+(r2-r2min)*v(2,1)
         end if
         
      end do
        
! ----- classical force for each DoF (NATOM3)
      dv = 0d0

      do n=1, nvpair2

        i = ivpair2(1,n)
        j = ivpair2(2,n)
         
        dx=-((q(3*j-2))+vpvec2(1, n)+(-q(3*i-2)))
        dy=-((q(3*j-1))+vpvec2(2, n)+(-q(3*i-1)))
        dz=-((q(3*j))  +vpvec2(3, n)+(-q(3*i))  )

        r2=dx*dx+dy*dy+dz*dz
        r = dsqrt(r2)

        ibin=int((r2-r2min)*binvrs)+1
          

        if (ibin.gt.1) then
           dr = (r2-r2min)-bin*dble(ibin-1)
           z = v(2,ibin)+(v(2,ibin+1)-v(2,ibin-1))/2d0*dr*binvrs
      
           dv(3*i-2) = dv(3*i-2)+z*dx*2d0
           dv(3*i-1) = dv(3*i-1)+z*dy*2d0
           dv(3*i)   = dv(3*i)  +z*dz*2d0

           dv(3*j-2) = dv(3*j-2)-z*dx*2d0
           dv(3*j-1) = dv(3*j-1)-z*dy*2d0
           dv(3*j)   = dv(3*j)  -z*dz*2d0

        else
!           dr = (r2-r2min)-bin*dble(ibin-1)
           z = v(2,1)

!           z = 0d0
           dv(3*i-2) = dv(3*i-2)+z*dx*2d0
           dv(3*i-1) = dv(3*i-1)+z*dy*2d0
           dv(3*i)   = dv(3*i)  +z*dz*2d0

           dv(3*j-2) = dv(3*j-2)-z*dx*2d0
           dv(3*j-1) = dv(3*j-1)-z*dy*2d0
           dv(3*j)   = dv(3*j)  -z*dz*2d0
!           write(*,*) "Out of range of r2."
         end if

       enddo         
!         end if

C         br2=bb*bb/r2
C
C         br5=br2*br2*sqrt(br2)
C
C         br52=br5/r2
C
C         dlng(3*i-2)=dlng(3*i-2)+2.5d0*br52*dx
C         dlng(3*i-1)=dlng(3*i-1)+2.5d0*br52*dy
C         dlng(3*i)  =dlng(3*i)  +2.5d0*br52*dz
C
C         d2lng(3*i-2)=d2lng(3*i-2)+2.5d0*br52*
C     *                             (1.0d0-7.0d0*dx**2/r2)
C         d2lng(3*i-1)=d2lng(3*i-1)+2.5d0*br52*
C     *                             (1.0d0-7.0d0*dy**2/r2)
C         d2lng(3*i)  =d2lng(3*i)  +2.5d0*br52*
C     *                             (1.0d0-7.0d0*dz**2/r2)


c --- now sum up the kinetic energy components.

!      do i=1, NATOM3
!         tloc=tloc+d2lng(i)+dlng(i)**2
!      end do

c --- account for mass factor and for double-counting of pairs.

!      tloc=-0.5d0*tloc/amass
!      vloc=0.5d0*vloc
      
!      print *,'vloc =',vloc
      return
      end subroutine

