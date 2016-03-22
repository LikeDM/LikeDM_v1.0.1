!############################################################################
!#
!# Note: Module of DM gamma fluxes computation  
!#
!# Author: Xiaoyuan Huang, Yue-Lin Sming Tsai, Qiang Yuan
!# Email: huangxiaoyuan@gmail.com, smingtsai@gmail.com, yuanq@pmo.ac.cn
!# Date: 2013-10-21
!# Date: 2014-08-13
!############################################################################


module dSphs_gamma
  use ReadTable
  use MathLib
  implicit none
  

contains
!=======================================================================




  !###########################################################
  ! likelihood subroutines for Gamma-ray from dSphs
  !
  ! input: mass mx and corss section sigmav
  ! output: DeltaChisq
  !###########################################################
  subroutine dSphs_calclike(mx,sigmav,totchisq)
    implicit none
    real*8, intent(in) :: mx,sigmav      !# DM mass (GeV), sigmav (ann: cm^3 s^-1 but decay: s)
    real*8, intent(out) :: totchisq
    real*8,parameter :: pi=dacos(-1.0d0)
    real*8 dnde,EsqF_,en_cv,tmp
    real*8 logJ(dsphs__ndsphs),logJerr(dsphs__ndsphs),logJtry(100)
    real*8 chisq_i,chisq_J,chisq_k,totchisq_k
    real*8 prefac
    integer i,j,k
    real*8 best_fit_logJ(dsphs__ndsphs)


    logJ(:)=dsphs__logJ(:)
    logJerr(:)=dsphs__logJerr(:)
    best_fit_logJ(:)=-5d3


    !... To setup the flux for giving mx,sigmav ...
    totchisq=0d0
    do i=1,dsphs__ndsphs,1

      if (dsphs__use(i)) then 
          call Gen_linspace(logJ(i)-3d0*logJerr(i),logJ(i)+3d0*logJerr(i),100,.False.,logJtry)

          chisq_i=1d80

          do j=1,100,1
            !--- total chisq for all the energy bins, giving logJtry  
            totchisq_k=0d0
            do k=1,dsphs__nEbin,1 
              en_cv=dsphs__en(k)
              call Gen_dnde_from_Table('gamma_ray',mx,en_cv,dnde) 

              if (sets__DMdecay) then
                prefac=sigmav/mx
              else
                prefac=sigmav/mx/mx/2d0
              endif

              EsqF_=10.0**(logJtry(j))*prefac*dnde*(en_cv**2.0d0)/(4.0d0*pi)   

 
              call linInt2p(dsphs__EsqF,dsphs__chisqMap(k,1:dsphs__nFbin,i),dsphs__nFbin,.true., &
                            EsqF_,chisq_k,.false.,.false.)
              totchisq_k=totchisq_k+chisq_k
            enddo

            chisq_J=((logJtry(j)-logJ(i))/logJerr(i))**2.0d0
            tmp=dlog(10.0d0)*dsqrt(2.0d0*pi)*logJerr(i)*(10.0d0**logJ(i)) !# normalization


            totchisq_k=totchisq_k+chisq_J-2.0d0*dlog(1.0d0/tmp)


            if (totchisq_k<chisq_i) then
              chisq_i=totchisq_k
              if (sets__seebug==5) best_fit_logJ(i)=logJtry(j)
            endif

          enddo

          totchisq=totchisq+chisq_i
      endif

    enddo


    if (dsphs__TS_bkg==0d0) call dSphs_bkg_only(dsphs__TS_bkg)
    totchisq=totchisq-dsphs__TS_bkg


    
    !*************************************************************
    !*
    !*       Debug flag for dSphs Gamma Ray  
    !*
    !*************************************************************
    !best_fit_logJ(:)=best_fit_logJ(:)
    write(*,*)'************************************************************'
    write(*,*)'--------- dSphs result: delta [chisq/-2ln(likelihood)] >>>>>'
    write(*,*) ' Fermi_dSphs:',totchisq
    !write(*,'(A,1E12.4)')'  TS(s+b)-TS(b) = ',totchisq
    !write(*,'(A,1E12.4)')'  TS(b) = ',dsphs__TS_bkg
    write(*,*)'--------- dSphs result: delta [chisq/-2ln(likelihood)] <<<<<'
    write(*,*)'************************************************************'


   if (sets__seebug==5) then 
     write(*,*)'************************************************************'
     write(*,*)'seebug=',sets__seebug,':  individual dSph spectrum >>>>>'
     do i=1,dsphs__ndsphs,1

       if (dsphs__use(i)) then 

         write(*,*)'---------------------------------------------------'
         write(*,*)'                dSph Name: ',trim(dsphs__lab(i))
         write(*,*)' The best fit log10[J factor] is:', best_fit_logJ(i)
         write(*,*)'               Energy(GeV) vs E^2dPhidE '
         write(*,*)'---------------------------------------------------'


         do k=1,dsphs__nEbin,1 
              en_cv=dsphs__en(k)
              call Gen_dnde_from_Table('gamma_ray',mx,en_cv,dnde) 
              write(*,'(2G12.4)')en_cv,10d0**best_fit_logJ(i)*prefac*dnde*(en_cv**2.0d0)/(4.0d0*pi)   
         enddo
       endif

     enddo 
     write(*,*)''

   endif



    return
  end subroutine dSphs_calclike



  subroutine dSphs_bkg_only(totchisq)
    implicit none
    real*8, intent(out) :: totchisq
    real*8,parameter :: pi=dacos(-1.0d0)
    real*8 tmp
    real*8 logJ(dsphs__ndsphs),logJerr(dsphs__ndsphs)
    real*8 chisq_J
    integer i,j,k

    logJ(:)=dsphs__logJ(:)
    logJerr(:)=dsphs__logJerr(:)

    !... To setup the flux for giving mx,sigmav ...
    totchisq=0d0
    chisq_J=0d0
    do i=1,dsphs__ndsphs,1

      if (dsphs__use(i)) then 

            do k=1,dsphs__nEbin,1 
              totchisq=totchisq + dsphs__chisqMap(k,1,i)
              !write(*,*)trim(dsphs__lab(i)),k,dsphs__chisqMap(k,1,i),dsphs__chisqMap(k,dsphs__nFbin,i), &
              !           minval(dsphs__chisqMap(k,1:dsphs__nFbin,i))

            enddo

            chisq_J=0d0
            tmp=dlog(10.0d0)*dsqrt(2.0d0*pi)*logJerr(i)*(10.0d0**logJ(i)) !# normalization
            totchisq=totchisq+chisq_J-2.0d0*dlog(1.0d0/tmp)

      endif

    enddo

    return
  end subroutine dSphs_bkg_only






!=======================================================================
end module dSphs_gamma

