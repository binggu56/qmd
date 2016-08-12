
c     QSATS version 1.0 (3 March 2011)

c     file name: eloc.f

c ----------------------------------------------------------------------

c     this computes the total energy and the expectation value of the
c     potential energy from the snapshots recorded by QSATS.

c ----------------------------------------------------------------------

      program eloc

      implicit double precision (a-h, o-z)

      include 'sizes.h'

      include 'qsats.h'
      
      common /bincom/ bin, binvrs, r2min


      character*4 :: atom(NATOMS)

c --- this common block is used to enable interpolation in the potential
c     energy lookup table in the subroutine local below.



      dimension q(NATOM3), vtavg(NREPS), vtavg2(NREPS),
     +          etavg(NREPS), etavg2(NREPS),dv(NATOM3)
      dimension idum(NATOM3)

      parameter (half=0.5d0)
      parameter (one=1.0d0)


      open(100, file='energy.dat')
      open(101, file='xoutput')
      open(102, file='traj4')
      open(103, file='pot.dat')
      open(104, file='force.dat')
      open(105, file='pes')
      open(106, file='temp.dat')

      Ndim = NATOM3

c --- initialization.

!      call tstamp

      write (6, 6001) NREPS, NATOMS, NATOM3, NATOM6, NATOM7,
     +                NVBINS, RATIO, NIP, NPAIRS
6001  format ('compile-time parameters:'//,
     +        'NREPS  = ', i6/,
     +        'NATOMS = ', i6/,
     +        'NATOM3 = ', i6/,
     +        'NATOM6 = ', i6/,
     +        'NATOM7 = ', i6/,
     +        'NVBINS = ', i6/,
     +        'RATIO  = ', f6.4/,
     +        'NIP    = ', i6/,
     +        'NPAIRS = ', i6/)

!      call input
      den = 4.61421d-3

      call vinit(r2min, bin)
!-----------------------------------------------------
     
      binvrs=one/bin

c --- read crystal lattice points.
      ltfile = 'ltfile'
      write (6, 6200) ltfile
6200  format ('READING crystal lattice from ', a16/)

      open (8, file='lattice-file-180', status='old')

      read (8, *) nlpts

      if (nlpts.ne.NATOMS) then
         write (6, *) 'ERROR: number of atoms in lattice file = ', nlpts
         write (6, *) 'number of atoms in source code = ', NATOMS
         stop
      end if

c --- read the edge lengths of the supercell.

      read (8, *) xlen, ylen, zlen

c --- compute a distance scaling factor.

      den0=dble(NATOMS)/(xlen*ylen*zlen)

c --- scale is a distance scaling factor, computed from the atomic
c     number density specified by the user.

      scale=exp(dlog(den/den0)/3.0d0)

      write (6, 6300) scale
6300  format ('supercell scaling factor computed from density = ',
     +        f12.8/)

      xlen=xlen/scale
      ylen=ylen/scale
      zlen=zlen/scale

      write (6, 6310) xlen, ylen, zlen
6310  format ('supercell edge lengths [bohr]         = ', 3f10.5/)

      dxmax=half*xlen
      dymax=half*ylen
      dzmax=half*zlen

      do i=1, NATOMS

         read (8, *) xtal(i, 1), xtal(i, 2), xtal(i, 3)

         xtal(i, 1)=xtal(i, 1)/scale
         xtal(i, 2)=xtal(i, 2)/scale
         xtal(i, 3)=xtal(i, 3)/scale

      end do

      close (8)

      write (6, 6320) xtal(NATOMS, 1), xtal(NATOMS, 2),
     +                xtal(NATOMS, 3)
6320  format ('final lattice point [bohr]            = ', 3f10.5/)

c --- this variable helps us remember the nearest-neighbor distance.

      rnnmin=-1.0d0

      do j=2, NATOMS

         dx=xtal(j, 1)-xtal(1, 1)
         dy=xtal(j, 2)-xtal(1, 2)
         dz=xtal(j, 3)-xtal(1, 3)

c ------ this sequence of if-then-else statements enforces the
c        minimum image convention.

         if (dx.gt.dxmax) then
            dx=dx-xlen
         else if (dx.lt.-dxmax) then
            dx=dx+xlen
         end if

         if (dy.gt.dymax) then
            dy=dy-ylen
         else if (dy.lt.-dymax) then
            dy=dy+ylen
         end if

         if (dz.gt.dzmax) then
            dz=dz-zlen
         else if (dz.lt.-dzmax) then
            dz=dz+zlen
         end if

         r=sqrt(dx*dx+dy*dy+dz*dz)

         if (r.lt.rnnmin.or.rnnmin.le.0.0d0) rnnmin=r

      end do

      write (6, 6330) rnnmin
6330  format ('nearest neighbor (NN) distance [bohr] = ', f10.5/)

      write (6, 6340) xlen/rnnmin, ylen/rnnmin, zlen/rnnmin
6340  format ('supercell edge lengths [NN distances] = ', 3f10.5/)

c --- compute interacting pairs.

      do i=1, NATOMS
         npair(i)=0
      end do

      nvpair=0

      do i=1, NATOMS
      do j=1, NATOMS

         if (j.ne.i) then

            dx=xtal(j, 1)-xtal(i, 1)
            dy=xtal(j, 2)-xtal(i, 2)
            dz=xtal(j, 3)-xtal(i, 3)

c --------- this sequence of if-then-else statements enforces the
c           minimum image convention.

            if (dx.gt.dxmax) then
               dx=dx-xlen
            else if (dx.lt.-dxmax) then
               dx=dx+xlen
            end if

            if (dy.gt.dymax) then
               dy=dy-ylen
            else if (dy.lt.-dymax) then
               dy=dy+ylen
            end if

            if (dz.gt.dzmax) then
               dz=dz-zlen
            else if (dz.lt.-dzmax) then
               dz=dz+zlen
            end if

            r2=dx*dx+dy*dy+dz*dz

            r=sqrt(r2)

c --------- interacting pairs are those for which r is less than a
c           certain cutoff amount. 

            if (r/rnnmin.lt.RATIO) then

               nvpair=nvpair+1

               ivpair(1, nvpair)=i
               ivpair(2, nvpair)=j

               vpvec(1, nvpair)=dx
               vpvec(2, nvpair)=dy
               vpvec(3, nvpair)=dz

               npair(i)=npair(i)+1

               ipairs(npair(i), i)=nvpair

            end if

         end if

      end do
      end do

      write (6, 6400) npair(1), nvpair
6400  format ('atom 1 interacts with ', i3, ' other atoms'//,
     +        'total number of interacting pairs = ', i6/)

c --- initialization.

c      loop=0
!      do k=1, NREPS
!         vtavg(k)=0.0d0
!         etavg(k)=0.0d0
!         vtavg2(k)=0.0d0
!         etavg2(k)=0.0d0
!      end do

!      open (10, file=spfile, form='unformatted')

c --- this loops reads the snapshots saved by QSATS.

!300   loop=loop+1

!     do k=1, NREPS, 11

!        read (10, end=600) (path(i, k), i=1, NATOM3)

!------reduce the number of pairs-------------
c      call reduce('REDUCE')

c      print *,'nvpair2=',nvpair2
c      print *,ivpair2(1,1),ivpair2(2,1)

c ------ compute the local energy and the potential energy.
      

      dx = 0.01d0

      q = 0d0
      call local(q,v0)

      q(1) = dx
      call local(q,v1)

      q(1) = 0d0
      q(4) = dx
      call local(q,v2)

      q = 0d0
      q(1) = dx
      q(4) = dx
      call local(q,v3)
     
      coup = (v3+v0-v1-v2)/dx**2
      write(*,1111) coup,v1-v0,v3-v0
1111  format('linear coupling coeff ',e14.6/,
     +       'deviation of one step ',e14.6/,
     +       'deviation of two steps',e14.6/)



      q = 0d0

      do j=1,201

        q(1) = -5d0+0.05d0*dble(j-1)

        call local(q,vloc)

        write(105,1000) q(1),vloc
      enddo


1000  format(20(e14.7,1x))
      return
      end program


c ----------------------------------------------------------------------

c     quit is a subroutine used to terminate execution if there is
c     an error.

c     it is needed here because the subroutine that reads the parameters
c     (subroutine input) may call it.

c ----------------------------------------------------------------------

      subroutine quit

      write (6, *) 'termination via subroutine quit'

      stop
      end subroutine
