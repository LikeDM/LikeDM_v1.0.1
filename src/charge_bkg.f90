!############################################################################
!#
!# Note: Module of DM e+ and antiproton fluxes computation  
!#
!# Author: Yue-Lin Sming Tsai, Qiang Yuan
!# Email: smingtsai@gmail.com, yuanq@pmo.ac.cn
!# Date: 2013-10-21
!# Date: 2014-08-13
!############################################################################


module charge_bkg
  use ReadTable
  use MathLib
  implicit none
  

contains
!=======================================================================

!###########################################################
!# Example to modulate the e+, e-, and antiproton fluxes
!###########################################################

  subroutine AfterModul(ptype,datapts,ein,flux,flux_mod)
     implicit none
     CHARACTER(LEN=*), intent(in) :: ptype 
     integer, intent(in) :: datapts   ! number of data point stored in the file
     real*8, intent(in) :: ein(datapts),flux(datapts)
     real*8, intent(out) :: flux_mod(datapts)
     !--- dummy varible ---
     real*8 smod       ! stored in module ReadTable 

     !--- calling the solar modulation ----
     if (trim(ptype).eq.'e+'.or.trim(ptype).eq.'e-') then
       smod=SMOD__ep
       call SOLAR_MOD(0,1,ein,flux,flux_mod,datapts,smod)

     elseif (trim(ptype).eq.'pbar') then 
       smod=SMOD__ap
       call SOLAR_MOD(1,-1,ein,flux,flux_mod,datapts,smod)

     endif

     return  
  end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++ solar modulation with phi (GV) of CRs in force field approximation
!++ Reference: Gleeson & Axford, 1968, ApJ, 154, 1011
!++ E and E_mod are kinetic energy per nucleon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine SOLAR_MOD(A,Z,E,F,F_mod,n,phi0)
    implicit none
    integer, intent(in) :: A,Z,n
    real*8, intent(in) ::  E(n),F(n)
    real*8, intent(in) :: phi0
    real*8, intent(out) :: F_mod(n)
    real*8 phi
    real*8 ELIS
    real*8 E_tmp(n),F_tmp(n)
    !--- dummy varibles ---
    real*8 m
    integer i
    real*8 ss

    if (A.ne.0) then
         m=0.9382d0
         phi=phi0*abs(Z)/A
    else
         m=0.511d-3
         phi=phi0
    endif

    if (phi>=0) then  !++ modulation
         do i=1,n,1
            E_tmp(i)=E(i)
            F_tmp(i)=F(i)
         enddo
         do i=1,n,1
            ELIS=E(i)+phi
            call linInt2p(E_tmp,F_tmp,n,.true.,ELIS,ss,.false.)
            F_mod(i)=(E(i)*(E(i)+2*m))/(ELIS*(ELIS+2*m))*ss
         enddo
    else  !++ demodulation
         do i=1,n,1
            E_tmp(i)=E(i)-phi
            F_tmp(i)=(((E(i)+m-phi)**2-m*m)/((E(i)+m)**2-m*m)*F(i))
         enddo
         do i=1,n,1
            if (E(i)>E_tmp(1)) then
               call linInt2p(E_tmp,F_tmp,n,.true.,E(i),ss,.false.)
               F_mod(i)=ss
            endif
         enddo
    endif

    return

  end subroutine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine phenom_model_epem(alpha,beta,En,cem,cep,cs,&
                                     gem1,gep1,gs,      &
                                     gem2,gep2,         &
                                     ebrem,ebrep,ebrs,  &                                  
                                     phiem,phiep)
    implicit none
    real*8, intent(in) :: alpha(3),beta(3)
    real*8, intent(in) :: En                 ! input energy 
    real*8, intent(in) :: cem,cep,cs         ! coefficient 
    real*8, intent(in) :: gem1,gep1,gs       ! index gamma_1
    real*8, intent(in) :: gem2,gep2          ! index gamma_2
    real*8, intent(in) :: ebrem,ebrep,ebrs   ! Energy E_{br}
    real*8, intent(out) :: phiem,phiep       ! final return result for ep,em,ap
    real*8 :: phiep_,phiem_, phis 

    phiem_=alpha(1)*cem*En**(beta(1)-gem1)/(1d0+(En/ebrem)**gem2)
    phiep_=alpha(2)*cep*En**(beta(2)-gep1)/(1d0+(En/ebrep)**gep2)
    phis=alpha(3)*cs*En**(beta(3)-gs) * dexp(-En/ebrs)

    phiep=phiep_+phis
    phiem=phiem_+0.6d0*phiep_+phis

    return

  end subroutine

  subroutine phenom_model_ap(alpha,beta,En,cap,gap1,gap2,ebrap,phiap)
    implicit none
    real*8, intent(in) :: alpha,beta
    real*8, intent(in) :: En         ! input energy 
    real*8, intent(in) :: cap        ! coefficient 
    real*8, intent(in) :: gap1,gap2  ! index gamma_1
    real*8, intent(in) :: ebrap      ! Energy E_{br}
    real*8, intent(out) :: phiap     ! final return result for ep,em,ap

    phiap=alpha*cap*En**(beta+gap1)/(1d0+En/ebrap)**gap2

    return

  end subroutine



!=======================================================================
end module charge_bkg

