c -----------------------------------------------
c     long-range correction for the potential
c     sum up all long-range interactions
c     units in Kelvin
c -----------------------------------------------

      program main
      implicit real*8(a-h,o-z)
      
      integer*4, parameter :: NL=18
      real*8,    parameter :: rcut = 13.8d0
c --- define unit vector in xy plane
      real*8 :: rx(2),ry(2)
      real*8 :: atom(3,50000)
      integer*4 :: il(NL)
      integer   :: ipair(1,50000)
c --- energy units converstion : hartree to Kelvin
      h2k = 3.157746d5

      open(10,file='xyz')
      open(11,file='pot')

      scale = 0.10487673d0

      do i=1,NL/2
        il(i) = -NL/2+(i-1)
        il(i+NL/2) = -il(i)
      enddo

      atom(1,1) = 0d0
      atom(2,1) = 0d0
      atom(3,1) = 0d0
      
      rx = (/0.0d0,0.7071068d0/)
      ry = (/0.6123722d0,0.3535534d0/) 
      
      nx = 16
      ny = 16

      k = 1
!      write(10,1000) atom(1,1),atom(2,1),atom(3,1)

      do i=-nx,nx,1
        do j=-ny,ny,1
          if (i**2 + j**2 /= 0) then
            k = k+1
            atom(1,k) = atom(1,1)+i*rx(1)+j*ry(1)
            atom(2,k) = atom(2,1)+i*rx(2)+j*ry(2)
            atom(3,k) = atom(3,1)
          endif
        enddo
      enddo

      do i=k+1,2*k,1
        atom(1,i) = atom(1,i-k)+0.4082485d0
        atom(2,i) = atom(2,i-k)+0.0d0
        atom(3,i) = atom(3,i-k)+0.5773507
      enddo

!
!      do i=2*k+1,3*k
!        atom(1,i) = atom(1,i-2*k)
!        atom(2,i) = atom(2,i-2*k)
!        atom(3,i) = atom(3,i-2*k)+1.1547004d0
!      enddo
!
!      do i=3*k+1,4*k
!        atom(1,i) = atom(1,i-2*k)
!        atom(2,i) = atom(2,i-2*k)
!        atom(3,i) = atom(3,i-2*k)+1.1547004d0
!      enddo
!
!      j=1
!      do i=4*k+1,5*k
!        atom(1,i) = atom(1,j)
!        atom(2,i) = atom(2,j)
!        atom(3,i) = atom(3,j)-1.1547004d0
!        j = j+1
!      enddo
!      
!      j = k
!      do i=5*k+1,6*k
!        j=j+1
!        atom(1,i) = atom(1,j)
!        atom(2,i) = atom(2,j)
!        atom(3,i) = atom(3,j)-1.1547004d0
!      enddo  


      do i=1,NL
        do j=1,k
          atom(1,2*i*k+j) = atom(1,j)
          atom(2,2*i*k+j) = atom(2,j)
          atom(3,2*i*k+j) = atom(3,j)+il(i)*1.1547004d0
          atom(1,(2*i+1)*k+j) = atom(1,j+k)
          atom(2,(2*i+1)*k+j) = atom(2,j+k)
          atom(3,(2*i+1)*k+j) = atom(3,j+k)+il(i)*1.1547004d0
        enddo
      enddo
      

      NATOMS = (NL+1)*2*k
     
      write(*,6002) k,NATOMS
6002  format('Num of atoms in one plane ',i6/ ,
     +       'Num of atoms in total ', i6/)

      atom = atom/scale

c --- print out the geometry file [a.u.]       
      write(10,6001) NATOMS
6001  format(i6,/)
      do i=1,natoms
        write(10,1000) atom(1,i),atom(2,i),atom(3,i)
      enddo      
1000  format('He ', 3(f12.6,1x))
c -----long-range pair interactions------------------

      NPAIR = 0
      do i=2,natoms
        r2 = (atom(1,i)-atom(1,1))**2+
     +       (atom(2,i)-atom(2,1))**2+
     +       (atom(3,i)-atom(3,1))**2
        r = sqrt(r2)
        if(r > rcut) then
!          vl = vl+hfdbhe(r2)
          NPAIR = NPAIR+1
          IPAIR(1,NPAIR)=i
        endif
      enddo
      
      write(*,6003) NPAIR
6003  format('Num of long-range interactions ',i6) 
c --- add displacements----------------
      do ii=1,60
        
        atom(1,1) = -3d0+(ii-1)*1d-1
        
        do jj=1,60
          
          atom(2,1) = -3d0+(jj-1)*1d-1
          
          do kk=1,60
            
            atom(3,1) = -3d0+(kk-1)*1d-1

          vlr = 0d0
          do i=1,NPAIR
            j = IPAIR(1,i)
            r2 = (atom(1,j)-atom(1,1))**2+
     +           (atom(2,j)-atom(2,1))**2+
     +           (atom(3,j)-atom(3,1))**2
            vlr = vlr +  hfdbhe(r2)
          enddo

c --- double counting------------------

          vlr = vlr/2d0    
        
          write(11,1001) atom(1,1),atom(2,1),atom(3,1),vlr*h2k
          
          enddo

        enddo
      enddo



6000  format('Long term potential correction [K] ',f10.6)

      

1001  format(20(e14.7,1x))
      stop
      end program 


c----------------------------------------------------------------------

c     this function computes the He-He potential energy as a function
c     of the squared interatomic distance.

c     the potential energy function is described in detail by R.A. Aziz
c     et al., Molecular Physics, volume 61, p. 1487 (1987).

c ----------------------------------------------------------------------

      double precision function hfdbhe(r2)

      implicit real*8 (a-h, o-z)

c --- parameters for the HFD(B) He-He potential energy function

      parameter (astar=1.8443101d5)
      parameter (alstar=10.43329537d0)
      parameter (bestar=-2.27965105d0)
      parameter (d=1.4826d0)
      parameter (c6=1.36745214d0)
      parameter (c8=0.42123807d0)
      parameter (c10=0.17473318d0)

      parameter (rm=5.59926d0)
      parameter (eps=10.948d0)

c --- hartree to kelvin conversion factor

      parameter (hart=315774.65d0)

      r=sqrt(r2)

      x=r/rm

      vstar=astar*exp(-alstar*x+bestar*x**2)

      vd=c6/x**6+c8/x**8+c10/x**10

      if (x.lt.d) vd=vd*exp(-(d/x-1.0d0)**2)

      hfdbhe=(vstar-vd)*eps/hart

      return
      end

     
