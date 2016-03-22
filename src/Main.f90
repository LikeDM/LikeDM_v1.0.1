!###############################################################
!#
!# Example to calculate likelihood (chisq)
!#                                                          
!# Author: Xiaoyuan Huang, Yue-Lin Sming Tsai, Qiang Yuan
!# Email: huangxiaoyuan@gmail.com, smingtsai@gmail.com, yuanq@pmo.ac.cn
!# Date: 2013-10-21    
!# Date: 2014-08-13             
!###############################################################


program LikeDM_fortran
  use MathLib
  use ReadTable
  use charge_lepton
  use charge_antip
  use dSphs_gamma
  implicit none

  integer i

  CHARACTER(LEN=300) buf
  CHARACTER(LEN=100) input_ini
  real*8 mx,sv
  logical use_pppc4
  integer WhatHalo,WhatGALPROP
  real*8 dSphs,MNBF(8)
  real*8 chisq(9),totchisq,tmp(100)
  real*8 aep(3),bep(3) ! at most 3 combinations  
  real*8 aap,bap
  integer nEbin


  !--- read ini-file and dnde from argv ---    
  CALL GETARG(1, buf )
  read(buf,'(A)') input_ini

  call read_setting(input_ini,mx,sv,use_pppc4,WhatHalo,WhatGALPROP)

  !#--- Gamma ray likelihood for dSphs ---
  if (sets__use_dSphs) then 
    call load_dsphs_LikeMap
    call dsphs_calclike(mx,sv,chisq(1))
  endif


  if (sets__use_ep .or. sets__use_ap) then 
    call loadGreen(WhatHalo,WhatGALPROP,sets__DMdecay)
    call Gen_epem_full_spectrum(mx,sv)
    call Gen_ap_full_spectrum(mx,sv)

    if (sets__use_ep) then
    call minuit_interface('positron',MNBF(1:6))
    aep(1:3)=MNBF(1:3)
    bep(1:3)=MNBF(4:6)

    call epemChisq('AMS02efr',aep,bep,chisq(3))
    call epemChisq('AMS02e+',aep,bep,chisq(4))
    call epemChisq('AMS02e-',aep,bep,chisq(5))
    call epemChisq('AMS02e+e-',aep,bep,chisq(6))


    endif

    if (sets__use_ap) then
    call minuit_interface('pbar',MNBF(1:6))  
    aap=MNBF(1)
    bap=MNBF(2)
    call apChisq(aap,bap,chisq(7))
    !write(*,*)'PAMELA_pbar:',chisq(7)
    endif
  endif 


  call print_debug_info(mx,sv,aep,bep,aap,bap,chisq)


end program LikeDM_fortran


!#################################################################
  subroutine read_setting(fname,mx,sv,use_pppc4,WhatHalo,WhatGALPROP)
!#################################################################
    use IniFile
    use ReadTable
    use PYTHIA_PPPC4
    implicit none
    CHARACTER(LEN=*), intent(in) :: fname 
    real*8, intent(out) :: mx,sv
    logical, intent(out) :: use_pppc4
    integer, intent(out) :: WhatHalo,WhatGALPROP
    CHARACTER(LEN=80) :: dnde_file
    CHARACTER(LEN=80) :: buf
    logical bad
    integer i


    call Ini_Open(fname, 1, bad, .false.)       !open input file
    if (bad) then 
      write(*,*)'Do not find input file', fname
      STOP
    endif
    sets__DMdecay=Ini_Read_Logical('decayDM',.false.)
    mx=Ini_Read_Double('DMmass',0.0d0)

    if (sets__DMdecay) then 
      sv=1d0/Ini_Read_Double('decay_time',0.0d0)
    else
      sv=Ini_Read_Double('sigmav',0.0d0)
    endif

    sets__seebug=Ini_Read_Int('seebug', 0)
    sets__halo=Ini_Read_Int('WhatHalo', 0)
    sets__use_dSphs=Ini_Read_Logical('use_dSphs',.false.)
    sets__use_ep=Ini_Read_Logical('use_ep',.false.)
    sets__use_ap=Ini_Read_Logical('use_ap',.false.)

    call game_initializing(Ini_Read_String('dsphs_map'), & 
                           Ini_Read_String('output_name') )

    SMOD__ep=Ini_Read_Double('epmod',0.5d0)
    SMOD__ap=Ini_Read_Double('apmod',0.5d0)

    WhatHalo=sets__halo
    WhatGALPROP=Ini_Read_Int('WhatGALPROP', 0)


  !--- set reshape parameters alpha, beta --- 
    if (sets__use_ep .or. sets__use_ap) then
      do i=1,4,1
  !--- cv step min max --- 
        GALPROP__alpha_prior(i)='1.0 0.02 0.1 10'
        GALPROP__beta_prior(i)='0.0 0.01 -0.5 0.5'
      enddo
    endif

    use_pppc4=Ini_Read_Logical('use_pppc4',.false.)
    if (use_pppc4) then
      do i=1,28,1
        write(buf,*) i
        PYTHIA__BR(i)=Ini_Read_Double('BR_'//trim(adjustl(buf)),0.0d0)
      enddo 
      call ReadPPPC4(PYTHIA__BR)
      if (sets__DMdecay) then 
        call dNdE_From_PPPC4(mx,mx)
      else
        call dNdE_From_PPPC4(mx,2d0*mx)
      endif

      call Load_dnde('copy from PPPC4')
    else
      call GETARG(2, buf )
      read(buf,'(A)') dnde_file
      call Load_dnde(trim(dnde_file))
    endif 


    call Ini_Close !close input file





    return
  end subroutine 



  subroutine minuit_interface(pX,BF)
    use ReadTable
    implicit none
      CHARACTER(LEN=*) pX
      integer nprs
      integer i,IERFLG
      external fcn
      integer NPRM(6)
      real*8 VSTRT(6),STP(6) 
      real*8 ARGLIS(12)
      real*8 UPLim(6),LOLim(6)
      real*8 BF(6)
      integer j,k
      real*8 ans
      CHARACTER*10 PNAM(6)
      DATA PNAM /'alpha(1)', 'alpha(2)', 'alpha(3)','beta(1)','beta(2)','beta(3)'/
      common /experiments/ nprs

      
      open (unit=64, file=trim(sets__foutput)//trim(pX)//'.mn.out', status='replace')

      CALL MNINIT(5,64,7)



      if (trim(pX).eq.'positron') then
        nprs=6
        do i=1,nprs/2,1
          !#cv, step, lower, upper    
          read(GALPROP__alpha_prior(i),*) VSTRT(i),STP(i),LOLim(i),UPLim(i)  
          read(GALPROP__beta_prior(i),*) VSTRT(i+nprs/2),STP(i+nprs/2),LOLim(i+nprs/2),UPLim(i+nprs/2) 
        enddo




      elseif (trim(pX).eq.'pbar') then
        read(GALPROP__alpha_prior(4),*) VSTRT(1),STP(1),LOLim(1),UPLim(1)  
        read(GALPROP__beta_prior(4),*) VSTRT(2),STP(2),LOLim(2),UPLim(2) 
        nprs=2
        PNAM(1)='alpha(4)'
        PNAM(2)='beta(4)'
      endif


  
      do i=1,nprs,1
        NPRM(i)=i
        CALL MNPARM(NPRM(i),PNAM(i),VSTRT(i),STP(i),LOLim(i),UPLim(i),IERFLG)
        IF (IERFLG .NE. 0) THEN
          WRITE (6,*) ' UNABLE TO DEFINE PARAMETER NO.',I
          STOP
        endif
      enddo

      CALL MNSETI('LikeDM minuit_interface')

      ARGLIS(1)=1 ! minimization strategy=1
      CALL MNEXCM(fcn,'SET STR',ARGLIS,1,IERFLG)
      ARGLIS(1)=500 ! maximum iteration
      ARGLIS(2)=0.01 ! tolerence
      CALL MNEXCM(fcn,'MIGRAD',ARGLIS,2,IERFLG)


      do i=1,nprs,1
        CALL MNPOUT(i,PNAM(i),ans,j,LOLim(i),UPLim(i),k)
        BF(i)=ans
      enddo
      close(5)
      close(64)
      close(7)
    return
  end subroutine




  subroutine fcn( npar, gin, totchisq, X, iflag, poly)
    use ReadTable
    use charge_lepton
    use charge_antip
    implicit none
    integer npar, iflag
    integer nprs
    real*8 gin(*), X(*)
    real*8 chi2, totchisq, poly
    real*8 alpha(3),beta(3)
    external poly
    common /experiments/ nprs
    
    if (nprs==6) then
      alpha(1) = X(1)
      alpha(2) = X(2)      
      alpha(3) = X(3)
      beta(1) = X(4)
      beta(2) = X(5)
      beta(3) = X(6)

      call epemChisq('AMS02efr',alpha,beta,chi2)
      totchisq= chi2
      call epemChisq('AMS02e+',alpha,beta,chi2)
      totchisq= totchisq+chi2
      call epemChisq('AMS02e-',alpha,beta,chi2)
      totchisq= totchisq+chi2
      call epemChisq('AMS02e+e-',alpha,beta,chi2)
      totchisq= totchisq+chi2

    elseif (nprs==2) then
      alpha(1) = X(1) 
      beta(1) = X(2)
      call apChisq(alpha(1),beta(1),chi2)
      totchisq= chi2
    endif
   
 
    return 
  end subroutine

!=======================================================================

