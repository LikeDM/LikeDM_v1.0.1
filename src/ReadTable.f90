!###############################################################
!#
!# Note: Module of ReadTable
!#
!# Author: Yue-Lin Sming Tsai, Qiang Yuan
!# Email: smingtsai@gmail.com, yuanq@pmo.ac.cn
!# Date: 2013-10-21
!# Date: 2014-08-13
!###############################################################
  
module ReadTable
  use MathLib
  use PYTHIA_PPPC4
  use charge_data
  implicit none
  
  !# runing setting  
   integer :: sets__seebug
   integer :: sets__halo
   logical :: sets__use_dSphs
   logical :: sets__use_ep
   logical :: sets__use_ap
   logical :: sets__use_DD
   character*80 :: sets__foutput
   character*80 :: sets__dSphs_MAP
   !--- flag for decay or annihilation --- 
   logical :: sets__DMdecay


  !# group pythia table 
   real*8 :: PYTHIA__BR(28)
   integer :: PYTHIA__nrows
   !--- Energy, gamma, e+, pbar ---
   real*8 :: PYTHIA__dnde(4,300)
   CHARACTER(LEN=100) PYTHIA__dnde_file



  !--- number of energy bins BEFORE propagation ---- 
  integer,parameter :: GALPROP__nm=50
  !--- number of energy bins AFTER propagation ---- 
  integer,parameter :: GALPROP__ne=99
  
  !--- Green function output --- 
  !--- 1: ep , 2: ap , 3: bkg e  ---
  real*8 GALPROP__Epost(3,GALPROP__ne),GALPROP__Epre(3,GALPROP__nm)
  real*8 GALPROP__GreenFun(3,GALPROP__ne,GALPROP__nm)

  integer :: GALPROP__prop_id
  character (LEN = 80) :: GALPROP__alpha_prior(4),GALPROP__beta_prior(4)

  real*8 :: SMOD__ep ! positron modulation potential
  real*8 :: SMOD__ap ! antiproton modulation potential
  
  !--- background for e+, e- and ap (not used!!!) ---
  !    Ek(GeV), bkg_elec(GeV^-1 m^-2 s^-1 sr^-1),bkg_posi, psr_posi, bkg_anti
  integer nbkp
  parameter(nbkp=99)
  real*8 charged_bgnd(5,nbkp)
  
  !--- The parameters of dSphs Likelihood map --- 
  integer,parameter :: dsphs__nEbin=24
  integer,parameter :: dsphs__nFbin=40
  integer,parameter :: dsphs__ndsphs=15  
  real*8 :: dsphs__en(dsphs__nEbin),dsphs__EsqF(dsphs__nFbin)
  real*8 :: dsphs__chisqMap(dsphs__nEbin,dsphs__nFbin,dsphs__ndsphs)  !# New chisq map including J unc. 


  integer :: dsphs__id_Bootes1
  integer :: dsphs__id_Canes_Venatici2 
  integer :: dsphs__id_Carina 
  integer :: dsphs__id_Coma
  integer :: dsphs__id_Draco
  integer :: dsphs__id_Fornax
  integer :: dsphs__id_Hercules
  integer :: dsphs__id_LeoII
  integer :: dsphs__id_LeoIV
  integer :: dsphs__id_Sculptor
  integer :: dsphs__id_Segue1
  integer :: dsphs__id_Sextans
  integer :: dsphs__id_UrsaMajor2
  integer :: dsphs__id_UrsaMinor
  integer :: dsphs__id_WillmanI

  logical :: dsphs__use(dsphs__ndsphs)
  real*8 :: dsphs__logJ(dsphs__ndsphs)
  real*8 :: dsphs__logJerr(dsphs__ndsphs)

  real*8 :: dsphs__TS_bkg

  CHARACTER*20 :: dsphs__lab(dsphs__ndsphs)



  
       
contains
!=======================================================================


  subroutine game_initializing(dSphs_MAP,foutput)
    use PYTHIA_PPPC4
    implicit none
    character*80, intent(in) :: foutput
    character*80, intent(in) :: dSphs_MAP 

    write(*,*)''
    write(*,*)' LikeDM (version 1.0)'
    write(*,*)''
    sets__foutput= trim(foutput)
    sets__dSphs_MAP= trim(dSphs_MAP)

    PreMod__DMpred(:,:,:)=0.0d0

    return
  end subroutine 



  !###########################################################
  ! Reading dnde from Table.
  ! Please keep in mind that this subroutine must be called always 
  ! right after change mx and dnde.  
  !###########################################################
  subroutine Load_dnde(fname)
    implicit none
    integer nrows,i
    CHARACTER(LEN=*), intent(in) ::  fname



    if (trim(fname).eq.'copy from PPPC4') then
      PYTHIA__nrows=179
      PYTHIA__dnde(1,1:PYTHIA__nrows)=PPPC4_dnde(1,1:PYTHIA__nrows)
      PYTHIA__dnde(2,1:PYTHIA__nrows)=PPPC4_dnde(2,1:PYTHIA__nrows)
      PYTHIA__dnde(3,1:PYTHIA__nrows)=PPPC4_dnde(3,1:PYTHIA__nrows)
      PYTHIA__dnde(4,1:PYTHIA__nrows)=PPPC4_dnde(4,1:PYTHIA__nrows)

    else

      i=1
      open(unit=58,file=trim(fname),status='old')
      do while(.true.)
         !--- Energy, gamma, e+, pbar ---
         read(unit=58,fmt=*,end=230) PYTHIA__dnde(1:4,i)
         i = i + 1
      enddo
230  close(58)
      nrows=i-1
      PYTHIA__nrows=nrows

    endif

    return
  end subroutine  
  
  

   
  !###########################################################
  ! Reading the chisq map for dSphs 10 sources
  !###########################################################
  subroutine load_dsphs_LikeMap()
    implicit none 
    integer i,j
    character (LEN = 10) :: buf
    real*8,parameter :: Msun=1.116d57    !# GeV
    real*8,parameter :: kpc2cm=3.086d21  !# cm

    !# nEbin=11, nFbin=31
    call Gen_linspace(0.57739d0,432.982156d0,dsphs__nEbin,.true.,dsphs__en)
    call Gen_linspace(1d-11,1d-7,dsphs__nFbin,.true.,dsphs__EsqF)


    !### reading the data ###
    open(unit=19,file=sets__dSphs_MAP,status='old')
    !read(19,*) buf
    do i=1,dsphs__nEbin,1
       do j=dsphs__nFbin,1,-1 ! inverse j order to make EsqF from small to big
         read(unit=19,fmt=*) dsphs__chisqMap(i,j,1:dsphs__ndsphs) 
       enddo !j
    enddo !i
    close(19)
  
    dsphs__chisqMap(:,:,:)=dsphs__chisqMap(:,:,:)*2d0
  
    !--- J-factor Table ---

    if (sets__DMdecay) then 
      ! Taken from Table II in paper 1504.02048 , unit Msun*kpc-2
      !dsphs__logJ=(/ 3.7d0,2.3d0,3.2d0,3.6d0,3.8d0,3.6d0,3.0d0,2.6d0,1.5d0,3.6d0,1.9d0,4.0d0,4.6d0,4.0d0,2.8d0 /)
      !dsphs__logJerr=(/ 0.4d0,0.5d0,0.1d0,0.5d0,0.3d0,0.1d0,0.4d0,0.2d0,0.6d0,0.1d0,1.1d0,0.2d0,0.5d0,0.1d0,0.5d0 /)

      dsphs__logJ=(/ 5.3d0,4.9d0,4.7d0,6.0d0,5.7d0,4.6d0,4.7d0,4.4d0,3.5d0,4.5d0,3.9d0,5.1d0,6.3d0,5.0d0,5.0d0 /)
      dsphs__logJerr=(/ 0.9d0,1.3d0,0.6d0,1.3d0,0.8d0,0.3d0,0.9d0,0.9d0,1.4d0,0.4d0,2.2d0,0.5d0,1.0d0,0.5d0,1.7d0 /)


      ! working in progress http://arxiv.org/pdf/1510.00389
      !                     http://arxiv.org/pdf/1408.0002v2

      do i=1,dsphs__ndsphs,1
        dsphs__logJ(i)= dsphs__logJ(i)+dlog10(Msun/kpc2cm/kpc2cm)
        dsphs__logJerr(i)= dsphs__logJerr(i)  
      enddo 

    else
      ! Taken from Table II in paper 1504.02048 , unit Msun*kpc-2
      dsphs__logJ=(/ 11.9d0,11.0d0,11.2d0,12.5d0,12.2d0,11.1d0,10.7d0,11.2d0,9.0d0,11.9d0,9.8d0,11.2d0,13.4d0,12.4d0,12.3d0 /)
      dsphs__logJerr=(/ 0.6d0,0.6d0,0.1d0,0.6d0,0.3d0,0.1d0,0.7d0,0.2d0,1.2d0,0.1d0,2.0d0,0.3d0,0.7d0,0.1d0,0.5d0 /)
      

      do i=1,dsphs__ndsphs,1
        dsphs__logJ(i)=dsphs__logJ(i)+dlog10(Msun**2/kpc2cm**5)
        dsphs__logJerr(i)= dsphs__logJerr(i)
      enddo 
    endif


    !~~~ to decide which dSph is to use and location in the table ~~~!
    dsphs__id_Bootes1=1
    dsphs__id_Canes_Venatici2=2
    dsphs__id_Carina=3 
    dsphs__id_Coma=4
    dsphs__id_Draco=5
    dsphs__id_Fornax=6
    dsphs__id_Hercules=7
    dsphs__id_LeoII=8
    dsphs__id_LeoIV=9
    dsphs__id_Sculptor=10
    dsphs__id_Segue1=11
    dsphs__id_Sextans=12
    dsphs__id_UrsaMajor2=13
    dsphs__id_UrsaMinor=14
    dsphs__id_WillmanI=15

    !dsphs__use(dsphs__id_Bootes1)=.true.
    !dsphs__use(dsphs__id_Carina)=.true.
    dsphs__use(1:15)=.true.



    dsphs__lab(dsphs__id_Bootes1)='Bootes I'
    dsphs__lab(dsphs__id_Canes_Venatici2)='Canes Venatici II'
    dsphs__lab(dsphs__id_Carina)='Carina'
    dsphs__lab(dsphs__id_Coma)='Coma'
    dsphs__lab(dsphs__id_Draco)='Draco'
    dsphs__lab(dsphs__id_Fornax)='Fornax'
    dsphs__lab(dsphs__id_Hercules)='Hercules'
    dsphs__lab(dsphs__id_LeoII)='Leo II'
    dsphs__lab(dsphs__id_LeoIV)='Leo IV'
    dsphs__lab(dsphs__id_Sculptor)='Sculptor'
    dsphs__lab(dsphs__id_Segue1)='Segue I'
    dsphs__lab(dsphs__id_Sextans)='Sextans'
    dsphs__lab(dsphs__id_UrsaMajor2)='Ursa Major II'
    dsphs__lab(dsphs__id_UrsaMinor)='Ursa Minor'
    dsphs__lab(dsphs__id_WillmanI)='Willman I'

    dsphs__TS_bkg=0d0


    return 
  end subroutine
  
  
  !###########################################################
  !# Reading the Green function for e+ and pbar. 
  !###########################################################
  subroutine loadGreen(WhatHalo,WhatGALPROP,DMdecay)
    implicit none 
    integer, intent(in) ::  WhatHalo, WhatGALPROP
    logical, intent(in) ::  DMdecay ! DM decay or annihilation  
    integer Fstate
    character (LEN = 80) :: table_name
  
    !--- dummy varible ---
    character (LEN = 80) :: string_tmp
    integer i
  
    GALPROP__prop_id=WhatGALPROP
  
    do Fstate=1,2,1
  
      if (Fstate==1) then 
        table_name='./dat/GreenFun/posi_'
      elseif (Fstate==2) then 
        table_name='./dat/GreenFun/pbar_'
      else
        write(*,*)'Fstate=1 is DM anni. positron.'
        write(*,*)'Fstate=2 is DM anni. antiproton.'
        write(*,*)'The rest is not supported yet.'
        stop
      endif
  
    
  
      if (WhatHalo==1) then 
        table_name=trim(table_name)//'nfw'
      elseif (WhatHalo==2) then 
        table_name=trim(table_name)//'ein'
      elseif (WhatHalo==3) then 
        table_name=trim(table_name)//'iso'
      else
        write(*,*)'WhatHalo=1 is NFW profile.'
        write(*,*)'WhatHalo=2 is Einasto profile.'
        write(*,*)'WhatHalo=3 is isothermal profile.'
        write(*,*)'The rest is not supported yet.'
        stop
      endif    

      if (WhatGALPROP<=6) then
        write(string_tmp,*) WhatGALPROP
        table_name=trim(table_name)//'_prop'//trim(adjustl(string_tmp))
      else
        write(*,*)'The WhatGALPROP>6 is not supported yet.'
        stop
      endif 


      !--- Decay or annihilation --- 

      if (DMdecay) then
        table_name=trim(table_name)//'_d.txt'
      else
        table_name=trim(table_name)//'.txt'
      endif 



  
      !--- Reading the table from file= table_name ---
      if (sets__seebug==1) write(*,*)' Reading Table from ',table_name
      open(unit=19,file=table_name,status='old')
      do i=1,GALPROP__ne,1
        !++ Epost(i) in [GeV], g(i,j) in [GeV^-1 m^-2 s^-1 sr^-1]
        read(19,*) GALPROP__Epost(Fstate,i),GALPROP__GreenFun(Fstate,i,1:GALPROP__nm)
      enddo
      close(19)
  
      do i=1,GALPROP__nm,1
        GALPROP__Epre(Fstate,i)=1d-2*10**(7.0/(GALPROP__nm-1)*(i-1)) ! GeV
      enddo
  
    enddo
  
    return 
  end subroutine
  
  
  !###########################################################
  !# Loading background e-, e+, pbar from physical model 
  !# (under construction)
  !###########################################################
!  subroutine load_bkg_phys(WhatGALPROP)
!    implicit none 
!    integer, intent(in) ::  WhatGALPROP
!    !--- dummy varible ---
!    character (LEN = 100) :: string_tmp
!    integer i
!    if (WhatGALPROP<=6) then
!      write(string_tmp,*) WhatGALPROP
!      string_tmp='./dat/ep_ap_bkg_expcon/bkg_prop'//trim(adjustl(string_tmp))//'.txt'
!      !--- fill in GALPROP parameters --- 
!      GALPROP__prop_id=WhatGALPROP
!      GALPROP__esmod=(/ 0.798,0.954,0.870,0.845,0.805,0.789 /)
!      GALPROP__psmod=(/ 0.32,0.34,0.36,0.36,0.36,0.339 /)
!    else
!      write(*,*)'The WhatGALPROP>6 is not supported yet.'
!      stop
!    endif 
! 
!    !--- To read background fluxes map ---
!    open(unit=19,file=string_tmp,status='old')
!    do i=1,nbkp,1
!      !Ek(GeV), bkg_elec(GeV^-1 m^-2 s^-1 sr^-1), bkg_posi, psr_posi, bkg_antip
!      read(19,*) charged_bgnd(1:5,i)
!    enddo
!    close(19)
!    return 
!  end subroutine
  


!###########################################################
! To interpolate dnde by given energy
! The hidden parameters: sets__DMdecay and PYTHIA__dnde
!###########################################################
subroutine Gen_dnde_from_Table(ptype,mx,en,dnde)
  implicit none
  CHARACTER(LEN=*), intent(in) :: ptype 
  real*8, intent(in) :: mx,en
  real*8, intent(out) :: dnde
  integer :: Tab_ID     
  integer i
  integer nrows
  real*8 half_sqrts
  real*8 x(PYTHIA__nrows),y(PYTHIA__nrows)
 
  nrows=PYTHIA__nrows

  if (sets__DMdecay) then 
    half_sqrts=mx*0.5d0
  else 
    half_sqrts=mx
  endif



  if (en>half_sqrts) then
     dnde=0.0d0
     return
  endif


  if (trim(ptype).eq.'gamma_ray') then
    Tab_ID=2
  elseif (trim(ptype).eq.'positron') then
    Tab_ID=3
  elseif (trim(ptype).eq.'antiproton') then
    Tab_ID=4
  else
    write(*,*) 'We are only provide e+/e-, photon, and antiproton spectrum...'
    write(*,*) 'Please use "gamma_ray", "positron", or "antiproton". '
    STOP
  endif


  !--- energy bins ---
  x(1:nrows)=PYTHIA__dnde(1,1:nrows)
  y(1:nrows)=PYTHIA__dnde(Tab_ID,1:nrows)



  ! if input energy smaller than all the data points' energy 
  if (en<x(1) .and. x(1)<x(nrows) ) then 
    dnde=y(1)

  ! if input energy greater than all the data points' energy 
  elseif (en>x(1) .and. x(1)>x(nrows) ) then
    dnde=y(1)

  ! if input energy greater than the last data points' energy
  ! The value of data points' energy is increased.  
  elseif (en>x(nrows) .and. x(1)<x(nrows) ) then
    dnde=y(nrows)

  ! if input energy greater than the last data points' energy
  ! The value of data points' energy is decreased. 
  elseif (en<x(nrows) .and. x(1)>x(nrows) ) then
    dnde=y(nrows)

  else

    do i=2,nrows,1

       if (x(i-1) < en .and. x(i) > en) then 
         call Between_two_pts(x(i-1),y(i-1),x(i),y(i),.true.,en,dnde)
         exit 
       elseif (x(i-1) > en .and. x(i) < en) then 
         call Between_two_pts(x(i-1),y(i-1),x(i),y(i),.true.,en,dnde)
         exit 
       endif

    enddo
  endif

  return
end subroutine



!###########################################################
! To interpolate and average dnde in [en-resol/2, en+resol/2]
!###########################################################
subroutine Gen_AVG_dnde(ptype,resol,mchi,en,dnde)
  implicit none
  CHARACTER(LEN=*), intent(in) :: ptype 
  real*8, intent(in) :: resol        ! energy resolution defined by user. 
  real*8, intent(in) :: mchi,en
  real*8, intent(out) :: dnde
  integer i   
  real*8 X(10),Y(10)
  real*8 area   


  call Gen_linspace(en-resol/2d0,en+resol/2d0,10,.false.,X)

  do i=1,10,1
    call Gen_dnde_from_Table(ptype,mchi,X(i),Y(i))
  enddo 

  call TableArea(X,Y,10,area)

  dnde=area/resol

  return
end subroutine
  








!=======================================================================
end module ReadTable
  
  
