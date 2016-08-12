c     QSATS version 1.0 (3 March 2011)

c     file name: sizes.h

c --- number of replicas or beads in the VPI polymer chain.

      parameter (NREPS=430)

c --- number of atoms in the system.

      parameter (NATOMS=180)

c --- various multiples of NATOMS.

      parameter (NATOM3=NATOMS*3)
      parameter (NATOM6=NATOMS*6)
      parameter (NATOM7=NATOMS*7)

c --- number of points on the interatomic potential energy curve, for
c     linear interpolation of the potential energy function.

      parameter (NVBINS=20000)

c --- "radius" of the interacting-pair region, in multiples of the
c     nearest-neighbor distance.

      parameter (RATIO=2.05)

c --- number of interacting pairs for each atom.

      parameter (NIP=56)

c --- total number of interacting pairs in the simulation box.

      parameter (NPAIRS=NATOMS*NIP)
