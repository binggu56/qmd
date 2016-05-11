!        subroutine fft_damp(data,Nt,dt,out)
      program main  
      implicit real*8(a-h,o-z)
      
      complex*16 :: im,z  
      complex*16,dimension(:),allocatable :: cor
      real*8, dimension(:),allocatable :: t
      complex*16 out(1000) 

!        common/four/Emin,Emax,de
!        Emax = 8d-2
        im = (0d0,1d0)
        Pi = 4.0d0*atan(1d0)
!        alfa = 1d0/2d0**6
!      alfa = 0d0 
      np = 200
      emin = 0d0 
      emax = 6d0
      de = (emax-emin)/(np-1)
!      alfa = 0d0
!        De = 2d0/(Nt*dt)
!        Emin = 0.015d0

      print *,'how many time steps?'
      read(*,*) nt 
!      print *,'what is time interval?'
!      read(*,*) dt 

!      nt = 4000 

      allocate(cor(nt),t(nt)) 

      open(10,file='cor.dat') 
      do i=1,nt 
        read(10,*) t(i),a,b,c
        cor(i) = a+im*b 
      enddo
      close(10)  


      dt = t(2)-t(1) 
      alfa = -log(1d-4)/t(nt)**2
      print *,'damping constant',alfa 

      out = (0d0,0d0)

        i = 1
        do while (Emin+(i-1)*De<Emax.and.i<=Nt)
        en = Emin + (i-1)*de
        do j=1,Nt
          t = 0d0 +j*dt
          out(i) = 2d0*real(cor(j)*exp(-im*en*t(j)))*exp(-alfa*t(j)**2)*dt + out(i)
        enddo
        i=i+1
        end do

        out = out/(2.d0*pi)

      open(100,file='fft.dat') 
      i = 1
      do while (Emin+(i-1)*De<Emax.and.i<=Nt)
        en = Emin + (i-1)*de
        write(100,1000) en,out(i) 
        i=i+1
      end do
            
1000  format(20f14.6,1x)

        stop
        end program 

        subroutine rfft_damp(data,Nt,dt,out)

        use cdat, only : pi, im 
        
        implicit real*8(a-h,o-z)
        
        real*8 :: data(Nt),out(Nt)

        common/four/Emin,Emax,de
!        Emax = 8d-2        
        alfa = 1d0/2d0**6
!      alfa = 0d0
!        De = 2d0/(Nt*dt)
!        Emin = 0.015d0
        out = 0d0 

        i = 1
        do while (Emin+(i-1)*De<Emax.and.i<=Nt)
        en = Emin + (i-1)*de
        do j=1,Nt
          t = 0d0 +j*dt
          out(i) = 2d0*real(data(j)*exp(im*en*t))*exp(-alfa*t**2)*dt + out(i)
        enddo
        i=i+1
        end do

        out = out/(2.d0*pi)

        return
        end subroutine


        subroutine fft(data,Nt,dt,out)
        implicit real*8(a-h,o-z)
        complex*16 :: im,data(Nt),out(Nt)
        common/four/Emin,Emax,de
!        Emax = 8d-2
        im = (0d0,1d0)
        Pi = 4.0d0*atan(1d0)
      alfa = 0d0

!        De = 2d0/(Nt*dt)
!        Emin = 0.015d0
        out = (0d0,0d0)
      
        i = 1
        do while (Emin+(i-1)*de<Emax.and.i<=Nt)
        en = Emin + (i-1)*de
          do j=1,Nt
            t = 0d0 +j*dt
            out(i) = 2d0*real(data(j)*exp(im*en*t))*dt + out(i)
          enddo
          i=i+1
        end do
        out = out/(2.d0*pi) 
        return
        end subroutine


! --- cosine FFT 
      subroutine fft_cos(data,Nt,dt,out)
      
      use cdat, only : pi, im 

      implicit real*8(a-h,o-z)

      complex*16 :: data(Nt),out(Nt)

      common/four/Emin,Emax,de
!      Emax = 8d-2
      alfa = 0d0
!        De = 2d0/(Nt*dt)
!        Emin = 0.015d0
        out = (0d0,0d0)
      
        i = 1
        do while (Emin+(i-1)*de<Emax.and.i<=Nt)
        en = Emin + (i-1)*de
          do j=1,Nt
            t = 0d0 +j*dt
            out(i) = 2d0*data(j)*cos(en*t)*dt + out(i)
          enddo
          i=i+1
        end do
        out = out/(2.d0*pi) 
        return
        end subroutine
      
      subroutine fft_sin_damp(data,Nt,dt,out)
      
      use cdat, only : pi, im 

      implicit real*8(a-h,o-z)

      real*8 :: data(Nt),out(Nt)

      common/four/Emin,Emax,de
      
      alfa = 1d0/2d0**6 
      
      out = 0d0 
      
        i = 1
        do while (Emin+(i-1)*de<Emax.and.i<=Nt)
        en = Emin + (i-1)*de
          do j=1,Nt
            t = 0d0 +j*dt
            out(i) = 2d0*data(j)*sin(en*t)*exp(-alfa*t**2)*dt + out(i)
          enddo
          i=i+1
        end do
        out = out/(2.d0*pi) 
        return
        end subroutine
