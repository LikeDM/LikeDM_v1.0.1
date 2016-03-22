!############################################################################
!#
!# Note: Module of antiproton fluxes computation  
!#
!# Author: Yue-Lin Sming Tsai, Qiang Yuan
!# Email: smingtsai@gmail.com, yuanq@pmo.ac.cn
!# Date: 2013-10-21
!# Date: 2015-06-30
!############################################################################


module charge_antip
  use ReadTable
  use MathLib
  use charge_bkg
  use charge_data
  implicit none
  


contains
!=======================================================================



!###########################################################
!# compute DM signal for pbar from Green function
!###########################################################
subroutine DM_signal_ap(mx,sv,npts,Ein,fluxout)
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
  real*8 resol            ! energy resolution defined by user. 
  integer :: Fstate       !Fstate=1 is e+ final state. Fstate=2 is anti-proton final state.

  Fstate=2

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

        call Gen_dnde_from_Table('antiproton',mx,X(k),dnde)

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
  !# Generate DM induced spectrum of antiproton after prop.
  !# (before solar modulation)
  !############################################################

  subroutine Gen_ap_full_spectrum(mx,sv)
    implicit none
    real*8, intent(in) :: mx,sv
    real*8 Ein(100),flux(100)
    real*8 Eup,Elo


    Elo=0.1d0
    Eup=2d3

    call Gen_linspace(Elo,Eup,100,.true.,Ein)
    call DM_signal_ap(mx,sv,100,Ein,flux)

    PreMod__DMpred(2,1,:)=Ein(:)
    PreMod__DMpred(2,2,:)=flux(:)



    return
  end subroutine 



  !###########################################################
  !# Giving alpha and beta to resphape diffuse background for 
  !# bkg_pbar then return theor which depends on expeiments 
  !# 3. Reshape(bkg_pbar)+DM_ap
  !###########################################################

  subroutine adjustShape_ap(alpha,beta,theor)
    implicit none
    real*8, intent(in) :: alpha,beta
    real*8, intent(out) :: theor(100)       !theoretical prediction based on Ecv
    integer i
    !--- dummy varible ---
    real*8 ap(100)
    real*8 ap_mod(100)
    real*8 phiap

    do i=1,100,1
      ap(i)=PreMod__DMpred(2,2,i)            ! DM_ap
    enddo

    call AfterModul('pbar',100,PreMod__DMpred(2,1,:),ap,ap_mod)

    !Flux(m^-2 s^-1 sr^-1 GeV^-1)
    do i=1,100,1
      call phenom_model_ap(alpha,beta,PreMod__DMpred(2,1,i),0.0995d0,1.844d0,5.077d0,2.849d0,phiap)
      theor(i)=ap_mod(i)+phiap
    enddo 

    return  
  end subroutine



  !-----------------------------------------------------------------
  ! Input: reshape parameters alpha, beta, pbar energy E(ne)
  ! Output: bkg_ap(ne), dm_ap(ne) (after solar mod.)
  !-----------------------------------------------------------------
  subroutine output_ap_atEarth(mx,sv,alpha,beta,ne,E,bkg_ap,dm_ap)
    implicit none
    integer, intent(in) :: ne
    real*8, intent(in) :: alpha,beta,mx,sv,E(ne)
    real*8, intent(out) :: bkg_ap(ne),dm_ap(ne)
    real*8 ap(ne)
    integer i

    call DM_signal_ap(mx,sv,ne,E,ap)
    call AfterModul('pbar',ne,E,ap,dm_ap)

    do i=1,ne,1
      call phenom_model_ap(alpha,beta,E(i),0.0995d0,1.844d0,5.077d0,2.849d0,bkg_ap(i))
    enddo

    return
  end subroutine

  !=======================================================================
  ! likelihood subroutines for PAMELA pbar experiments from halo. 
  ! Here we consider reshaping fluxes by using a factor a*E^b. 
  ! a:[0.1,10], b[-0.5,0.5] 
  !=======================================================================


  subroutine apChisq(alpha,beta,chisq)
    real*8, intent(in) :: alpha,beta
    real*8, intent(out) :: chisq
    real*8 theor(100),pflux(100)  !theoretical prediction based on Ecv

    pflux(:)=0d0
    !--- Re-shape the spectrum with alpha, beta, and solar module ---
    call adjustShape_ap(alpha,beta,theor)

    call PAMELA_antip(PreMod__DMpred(2,1,:),pflux,theor,chisq)

      

    !# very slow!
    if (sets__seebug>=6) then
      write(*,*)''
      write(*,*)'--------- antiproton result >>>>>'
      write(*,'(A,2G12.4)')'alpha,beta=',alpha,beta
      write(*,'(A,1E12.4)')' chisq/-2ln(likelihood)= ',chisq
      write(*,*)''
    endif


    return  
  end subroutine




!=======================================================================
end module charge_antip

