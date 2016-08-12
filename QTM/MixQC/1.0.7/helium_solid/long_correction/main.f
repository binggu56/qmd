      program main
      implicit real*8(a-h,o-z)
      
      include 'size.h'
      include 'qsats.h'
      integer*4, parameter :: NATOMS = 180
      real*8,    parameter :: rcut = 13.8d0

      real*8 atom(3,NATOMS)


      
