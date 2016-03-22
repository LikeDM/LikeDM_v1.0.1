!############################################################################
!#
!# Note: Module of DM e+ and antiproton fluxes computation  
!#
!# Author: Yue-Lin Sming Tsai, Qiang Yuan
!# Email: smingtsai@gmail.com, yuanq@pmo.ac.cn
!# Date: 2013-10-21
!# Date: 2015-06-29
!############################################################################


module charge_lepton
  use ReadTable
  use MathLib
  use charge_bkg
  use charge_data
  implicit none
  


contains




!###########################################################
!# compute DM signal for e+ and pbar from Green function. 
!###########################################################
subroutine DM_signal_epem(mx,sv,npts,Ein,fluxout)
  implicit none 
  integer, intent(in) ::  npts
  real*8, intent(in) :: Ein(npts)
  real*8, intent(out) ::  fluxout(npts)
  !--- DM property ---
  real*8, intent(in) :: mx           
  real*8, intent(in) :: sv !sv: annihilation
  real*8 dnde

  !--- Flux function ---
  real*8 DMF(GALPROP__ne)
  !--- dummy varible ---
  integer i,j,k
  integer, parameter :: nb=20
  real*8 DMpart(nb),AstroPart(nb),DMAstro,X(nb),Y(nb)
  integer :: Fstate       !Fstate=1 is e+ final state. Fstate=2 is anti-proton final state.
  real*8 resol            ! energy resolution defined by user. 

  Fstate=1

  !--- Green function output from GALPROP ---- 
  !    Epost(3,GALPROP__ne),Epre(3,GALPROP__nm), GreenFun(3,GALPROP__ne,GALPROP__nm)

  do i=1,GALPROP__ne,1
    DMF(i)=0.0
    do j=1,GALPROP__nm,1

      resol=1.1788 ! square root of kernel bin width
   
      !--- compute dnde before propagation ---
      call Gen_linspace(GALPROP__Epre(Fstate,j)/resol,GALPROP__Epre(Fstate,j)*resol,nb,.false.,X)

      !--- average the contribution inside the kernel bin ---
      do k=1,nb,1

        call Gen_dnde_from_Table('positron',mx,X(k),dnde) ! 2 for e+
        !--- compute DM particle annihilation part contribution ---
        if (sets__DMdecay) then
          DMpart(k)=dnde*(sv/1d-26)/(mx/2d3)
        else
          DMpart(k)=dnde*(sv/1d-26)/(mx/1d3)**2
        endif

        !--- compute annihilation green function for dnde after propagation
        AstroPart(k)=0d0
        if (X(k)<GALPROP__Epre(Fstate,j)) then
          if (j==1) then
            call Between_two_pts(GALPROP__Epre(Fstate,j),GALPROP__GreenFun(Fstate,i,j), &
                                 GALPROP__Epre(Fstate,j+1),GALPROP__GreenFun(Fstate,i,j+1),&
                                 .true.,X(k),AstroPart(k))
          else
            call Between_two_pts(GALPROP__Epre(Fstate,j),GALPROP__GreenFun(Fstate,i,j), &
                                 GALPROP__Epre(Fstate,j-1),GALPROP__GreenFun(Fstate,i,j-1),&
                                 .true.,X(k),AstroPart(k))
          endif ! j==1
        else
          if (j==GALPROP__nm) then
            call Between_two_pts(GALPROP__Epre(Fstate,j),GALPROP__GreenFun(Fstate,i,j), &
                                 GALPROP__Epre(Fstate,j-1),GALPROP__GreenFun(Fstate,i,j-1),&
                                 .true.,X(k),AstroPart(k))
          else
            call Between_two_pts(GALPROP__Epre(Fstate,j),GALPROP__GreenFun(Fstate,i,j), &
                                 GALPROP__Epre(Fstate,j+1),GALPROP__GreenFun(Fstate,i,j+1),&
                                 .true.,X(k),AstroPart(k))
          endif ! j==GALPROP__nm
        endif
        Y(k)=DMpart(k)*AstroPart(k)
      enddo
      call TableArea(X,Y,nb,DMAstro)

      !--- integrate all the parameter space ---
      DMF(i)=DMF(i)+DMAstro

    enddo
  enddo

  !--- intepolation ---
  do i=1,npts,1
    call linInt2p(GALPROP__Epost(2,:),DMF,GALPROP__ne,.true.,Ein(i),fluxout(i),.false.) 
  enddo


  return 
end subroutine




  !############################################################
  !# Generate DM induced spectrum of e+/e- after prop. with respect 
  !# to experimental energy points (before solar modulation)
  !############################################################

  subroutine Gen_epem_full_spectrum(mx,sv)
    implicit none
    real*8, intent(in) :: mx,sv
    real*8 Ein(100),flux(100)
    real*8 Eup,Elo

    Elo=0.1d0
    Eup=2d3

    call Gen_linspace(Elo,Eup,100,.true.,Ein)
    call DM_signal_epem(mx,sv,100,Ein,flux)

    PreMod__DMpred(1,1,:)=Ein(:)
    PreMod__DMpred(1,2,:)=flux(:)


    return
  end subroutine 







  !###########################################################
  !# Giving alpha and beta to resphape diffuse background for 
  !# bkg_elec, bkg_posi, psr_posi
  !# then return theor which depends on expeiments 
  !# 1. Reshape(bkg_elec)+Reshape(bkg_posi)+Reshape(psr_posi)+ DM_e+, 
  !# 2. Reshape(bkg_posi)+Reshape(psr_posi)+ DM_e+,
  !###########################################################

  subroutine adjustShape_epem(alpha,beta,ep_mod,em_mod)
    implicit none

    real*8, intent(in) :: alpha(3),beta(3) ! at most 3 combinations  
    real*8, intent(out) ::  ep_mod(100),em_mod(100)

    call ep_em_atEarth(100,alpha,beta,ep_mod,em_mod)

    return  
  end subroutine


  !-----------------------------------------------------------------
  !
  !  To compute the total fluxes of e+ and e- at the earth  
  !
  !-----------------------------------------------------------------
  subroutine ep_em_atEarth(pts,alpha,beta,ep_mod,em_mod)
    implicit none
    integer, intent(in) :: pts
    real*8, intent(in) :: alpha(3),beta(3) ! at most 3 combinations 
    real*8, intent(out) ::  ep_mod(100),em_mod(100) 
    real*8 ep(100),em(100)
    real*8 phiem,phiep
    integer i

    ep(:)=PreMod__DMpred(1,2,:)              ! DM_posi      
    em(:)=PreMod__DMpred(1,2,:)              ! DM_elec

    call AfterModul('e+',pts,PreMod__DMpred(1,1,:),ep,ep_mod)
    call AfterModul('e-',pts,PreMod__DMpred(1,1,:),em,em_mod)

    do i=1,pts,1
       call phenom_model_epem(alpha,beta,PreMod__DMpred(1,1,i),     &
            21.6701d0,1.4991d0,0.6526d0,0.9344d0,0.9024d0,2.3390d0, &
            2.3734d0,2.3647d0,3.6390d0,2.8434d0,6.5289d2,           &
            phiem,phiep)
       ep_mod(i)=ep_mod(i)+phiep
       em_mod(i)=em_mod(i)+phiem
    enddo

    return  
  end subroutine


  !-----------------------------------------------------------------
  ! Input: reshape parameters alpha(3), beta(3), e+/e- energy E(ne)
  ! Output: bkg_em(ne), bkg_ep(ne), dm_ep(ne) (after solar mod.)
  !-----------------------------------------------------------------
  subroutine output_ep_em_atEarth(mx,sv,alpha,beta,ne,E,bkg_em,bkg_ep,dm_ep)
    implicit none
    integer, intent(in) :: ne
    real*8, intent(in) :: alpha(3),beta(3) ! at most 3 alpha/beta  
    real*8, intent(in) :: mx,sv,E(ne)
    real*8, intent(out) :: bkg_em(ne),bkg_ep(ne),dm_ep(ne)
    real*8 ep(ne)
    integer i

    call DM_signal_epem(mx,sv,ne,E,ep)
    call AfterModul('e+',ne,E,ep,dm_ep)

    do i=1,ne,1
       call phenom_model_epem(alpha,beta,E(i),                      &
            21.6701d0,1.4991d0,0.6526d0,0.9344d0,0.9024d0,2.3390d0, &
            2.3734d0,2.3647d0,3.6390d0,2.8434d0,6.5289d2,           &
            bkg_em(i),bkg_ep(i))
    enddo

    return  
  end subroutine

!=======================================================================
! likelihood subroutines for 4 e+/e- experiments from AMS02.
! Here we consider reshaping fluxes by using a factor a*E^b. 
! a:[0.1,10], b[-0.5,0.5] 
! Since it is required to connect Minuit, we have to just group to two 
! datasets, e+e- and antiproton  
!=======================================================================


  subroutine epemChisq(dataset,alpha,beta,chisq)
    CHARACTER(LEN=*), intent(in) :: dataset 
    real*8, intent(in) :: alpha(3),beta(3)
    real*8, intent(out) :: chisq
    real*8 chisq_
    real*8 ep_mod(100) ,em_mod(100)  !theoretical prediction


    !--- Re-shape the spectrum with alpha, beta, and solar module ---
    call adjustShape_epem(alpha,beta,ep_mod,em_mod)

    chisq=0.0d0

    if (trim(dataset).eq.'AMS02efr') then

      !# AMS02 positron fraction (PRL, 2013, 110, 141102)
      call AMS02_posifrac(PreMod__DMpred(1,1,:),em_mod,ep_mod,chisq_)
      chisq=chisq+chisq_

    elseif (trim(dataset).eq.'AMS02e+') then 

      !# AMS02 positron spectra (PRL, 2014, 113, 121102)
      call AMS02_positron(PreMod__DMpred(1,1,:),em_mod,ep_mod,chisq_)
      chisq=chisq+chisq_

    elseif (trim(dataset).eq.'AMS02e-') then 

      !# AMS02 electron spectra (PRL, 2014, 113, 121102)
      call AMS02_electron(PreMod__DMpred(1,1,:),em_mod,ep_mod,chisq_)
      chisq=chisq+chisq_


    elseif (trim(dataset).eq.'AMS02e+e-') then 

      !# AMS02 total electron/positron spectra (PRL, 2014, 113, 221102)
      call AMS02_total_elec(PreMod__DMpred(1,1,:),em_mod,ep_mod,chisq_)
      chisq=chisq+chisq_

    else
      write(*,*)'We do not have this dataset:', dataset
      stop
    endif



    !# very slow!  
    if (sets__seebug>=6) then
      write(*,*)''
      write(*,*)'--------- ',dataset,' result >>>>>'
      write(*,'(A,3G12.4)')'alpha[1:3]=',alpha(:)      
      write(*,'(A,3G12.4)')'beta[1:3]=',beta(:)
      write(*,'(A,1E12.4)')' chisq/-2ln(likelihood)= ',chisq
      write(*,*)''
    endif


    return  
  end subroutine



!=======================================================================
end module charge_lepton

