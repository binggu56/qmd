        subroutine derivs(x,Ntraj,fx,v,ddv)
        implicit real*8(a-h,o-z)
        parameter(ak=1.2d0)
        integer*4,intent(IN) :: Ntraj
        real*8,intent(IN) :: x(Ntraj)
        real*8,intent(OUT) :: fx(Ntraj),v(Ntraj),ddv(Ntraj)
        
        do i=1,Ntraj
          fx(i) = -ak*x(i)
          v(i)  = ak*x(i)**2/2d0
          ddv(i) = ak
        enddo

        return
        end subroutine

