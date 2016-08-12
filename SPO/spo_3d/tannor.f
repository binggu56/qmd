! Required common block for using this subroutine:
! common/mass/ma,mb,mc
! common/iw/iwc,iws [file numbers for output cf(Nt) and Sab]
	subroutine tw(corr,alpha,beta,px0,pxp,x0,xp,Nt,dt,sab)
	implicit real*8(a-h,o-z)
	real*8 ma,mb,mc,ka,kb
	real*8 m_at,m_bt,m_av,m_bv,w_at,w_bt
	real*8 at,av,bt,bv,pi,v_a,v_b
	complex*16 eta_a,eta_b,sab,im
	complex*16 corr(Nt),cf(Nt)
	common/mass/ma,mb,mc
!	common/width/alfa,beta
!	common/par_p/px0,py0,pxp,pyp
	common/four/Emin,Emax,de
	common/iw/iwc,iws
!	open(26,file='corf')
!	open(27,file='smat')
	im=(0d0,1d0)
        fra = 4401.2d0
        frb = 4401.2d0
	pi = 4.d0*atan(1.d0)
	sol= 3.d10
        au_s = 2.418884d-17

        m_at = ma*(mb+mc)/(ma+mb+mc)
        m_bt = (ma+mb)*mc/(ma+mb+mc)
        m_av = mb*mc/(mb+mc)
        m_bv = ma*mb/(ma+mb)
        print *,'Mass for translational motion in channel a =',m_at
        print *,'Mass for vibrational motion in channel a =',m_av
        print *,'Mass for translational motion in channel b =',m_bt
        print *,'Mass for vibrational motion in channel b =',m_bv

        at = dsqrt(dsqrt(2d0*alpha/pi))
        av = dsqrt(dsqrt(2d0*beta/pi))
        bv = dsqrt(dsqrt(2d0*beta/pi))
        bt = dsqrt(dsqrt(2d0*alpha/pi))
        print *,'Normalization constants'
        print *,'Norm_r1',av
        print *,'Norm_R1',at
        print *,'Norm_r2',bv
        print *,'Norm_R2',bt


        call fourier(corr,Nt,dt,cf)
        do i=1,Nt
        en = Emin + (i-1)*de
        write(iwc,1000) en,cf(i)
        enddo
! asymptotic potential for reactant and product, v_r & v_p      
        v_a = 0.d0
        v_b = 0.d0
	w_at=alpha
	w_bt=alpha
! transfer frequency to ground state energy
        Ea = Pi*sol*fra*au_s
        Eb = Pi*sol*frb*au_s
        print *,'Ground State Energy'
        print *,'E0(H-H) =',Ea
        print *,'E0(O-H) =',Eb
! calculate s-matrix
        i = 1
        do while (Emin+(i-1)*de < Emax)
        en = Emin + (i-1)*de
        E_at = En-Ea-V_a
        E_bt = En-Eb-V_b
! judge if k is lower than 0
!        if( E_at < 0d0 .or. E_bt < 0d0) then
!        print *,'ERROR: Energy minimum is too small.'
!        endif
        ka = dsqrt(2d0*m_at*E_at)
        kb = dsqrt(2d0*m_bt*E_bt)
        cona = at*dsqrt(m_at/(2d0*ka*w_at))
        conb = bt*dsqrt(m_bt/(2d0*kb*w_bt))
        eta_a=exp(-(ka+px0)**2/(4d0*w_at)+im*x0*ka)*cona
        eta_b=exp(-(kb-pxp)**2/(4d0*w_bt)-im*xp*kb)*conb
        sab = cf(i)/(conjg(eta_b)*eta_a)
        write(iws,1000) en,eta_a,eta_b,sab,abs(sab)**2
        i=i+1
        end do

1000	format(20(e14.7,1x))
	return
        end subroutine

