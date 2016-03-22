!###############################################################
!#
!# Note: Module of pythia table interpolation (PPPC4)  
!#
!# Author: Yue-Lin Sming Tsai
!# Email: smingtsai@gmail.com
!# Date: 2013-11-04
!###############################################################
module PYTHIA_PPPC4
  use MathLib
  implicit none

    integer :: PPPC4_Mbin
    integer :: PPPC4_Ebin
    real*8 :: PPPC4_spect(3,62,179)   ! dnde from PPPC4 table products, [mx bins, energy bins].
    real*8 :: PPPC4_K_over_mx(179)
    real*8 :: PPPC4_mx(62)
    real*8 :: PPPC4_dnde(4,179)       ! dnde for LikeDM formate
contains
!--- The dictionary in PPPC4 table ---
! 1. eL       2. eR           3. e            
! 4. \[Mu]L   5. \[Mu]R       6. \[Mu]        
! 7. \[Tau]L  8. \[Tau]R      9. \[Tau]       
! 10. q       11. c           12. b        13. t            
! 14. WL      15. WT          16. W            
! 17. ZL      18. ZT          19. Z            
! 20. g       21. \[Gamma]    22. h            
! 23. \[Nu]e  24. \[Nu]\[Mu]  25. \[Nu]\[Tau]   
! 26. V->e    27. V->\[Mu]    28. V->\[Tau]


!--- Fluxes at production, including EW corrections, taken from Sample.nb ---
! Fluxes are given as Log10[dN/dlog10(x)]  
! (x = K/mx, with K the kinetic energy in GeV)  
! normalized per single DM annihilation. 
! The format is
!
!	dlNdlxIEW[primary -> secondary][mass, lx]
!
!where
!
!	primary = eL, eR, e, \[Mu]L, \[Mu]R, \[Mu], \[Tau]L, \[Tau]R, \[Tau], q, c, b, t, WL, WT, W, ZL, ZT, Z, g, 
!                \[Gamma], h, \[Nu]e, \[Nu]\[Mu], \[Nu]\[Tau], V->e, V->\[Mu], V->\[Tau]
!	secondary = e, p, \[Gamma], d, \[Nu]e, \[Nu]\[Mu], \[Nu]\[Tau]	
!	mass = mx in GeV, on the range mx = 5 GeV -> 100 TeV
!	lx = Log[10, x ], on the range x = 10^-9-> 1 except for the channels V->e, V->\[Mu], V->\[Tau] for which x = 10^-5-> 1
!
  !--- The cols in PPPC4 pythia table ---
  ! 1. mDM       2. Log10[K/mx]   
  ! 3. eL        4. eR           5. e            
  ! 6. \[Mu]L    7. \[Mu]R       8. \[Mu]        
  ! 9. \[Tau]L   10. \[Tau]R     11. \[Tau]       
  ! 12. q(u,d,s) 13. c           14. b        15. t            
  ! 16. WL       17. WT          18. W            
  ! 19. ZL       20. ZT          21. Z            
  ! 22. g        23. \[Gamma]    24. h            
  ! 25. \[Nu]e   26. \[Nu]\[Mu]  27. \[Nu]\[Tau]   
  ! 28. V->e     29. V->\[Mu]    30. V->\[Tau]
  ! http://arxiv.org/pdf/1012.4515v4.pdf

  subroutine ReadPPPC4(BR)
    implicit none   
    real*8, intent(in) :: BR(28)
    real*8  line(30)
    character (LEN = 80) :: buf
    integer i,j,k,ch
    real*8 kin
    !--- dummy varibles ---
    real*8 ans

    PPPC4_spect(:,:,:)=0.0d0
    PPPC4_Ebin=179
    PPPC4_Mbin=62


    do ch=1,3,1  !

      if (ch==1) open(unit=58,file='./dat/PPPC4/AtProduction_gammas.dat',status='old')
      if (ch==2) open(unit=58,file='./dat/PPPC4/AtProduction_positrons.dat',status='old')
      if (ch==3) open(unit=58,file='./dat/PPPC4/AtProduction_antiprotons.dat ',status='old')

      read(unit=58,fmt=*) buf
      do i=1,PPPC4_Mbin,1
        do j=1,PPPC4_Ebin,1

         read(58,*) line(1:30)
         if (j==1 .and. ch==1) PPPC4_mx(i)=line(1)
         if (i==1 .and. ch==1) PPPC4_K_over_mx(j)=10.0d0**line(2)


         kin=PPPC4_K_over_mx(j)*PPPC4_mx(i)
         do k=1,28,1
           ans= line(k+2)/(kin*dlog(10.0d0))
           PPPC4_spect(ch,i,j)=PPPC4_spect(ch,i,j)+BR(k)*ans
         enddo

        enddo
      enddo
      close(58)
    enddo

    return
  end subroutine



  subroutine dNdE_From_PPPC4(mx_in,sqrts)
    implicit none  
    real*8, intent(in) :: mx_in
    real*8,optional,intent(in) :: sqrts
    real*8 mx,mx_lt_5GeV
    integer lab(2)
    !--- dummy varibles ---
    real*8 ans,x(2),y(2)
    integer i,j

    !--- decaying DM or annihilation DM ---           
    ! since PPPC4 is only for annihilation, the sqrt(s)/2=mx is required for PPPC4 mx entry.
    ! This is not ture for decaying DM. This entry is sqrt(s)/2=mx/2
    if (present(sqrts)) then
      mx=sqrts/2d0
    else
      mx=mx_in   
    endif

    mx_lt_5GeV=mx
    if (mx<5d0) mx=5d0 
      

    !--- find out the two points near mx ---
    do i=1,PPPC4_Mbin-1,1



      if (PPPC4_mx(i)==mx) then 

        PPPC4_dnde(1,:)=PPPC4_K_over_mx(:)*mx
        PPPC4_dnde(2,:)=PPPC4_spect(1,i,:)
        PPPC4_dnde(3,:)=PPPC4_spect(2,i,:)
        PPPC4_dnde(4,:)=PPPC4_spect(3,i,:)

        ! Make the correction for mx< 5 GeV
        if (mx_lt_5GeV<5d0)  then 

           PPPC4_dnde(1,:)=PPPC4_K_over_mx(:)*mx_lt_5GeV
           PPPC4_dnde(2,:)=PPPC4_spect(1,i,:)*mx/mx_lt_5GeV
           PPPC4_dnde(3,:)=PPPC4_spect(2,i,:)*mx/mx_lt_5GeV
           PPPC4_dnde(4,:)=PPPC4_spect(3,i,:)*mx/mx_lt_5GeV

        endif

        return

      elseif (PPPC4_mx(i+1)>=mx .and. PPPC4_mx(i)<=mx) then 
        lab(1)=i
        lab(2)=i+1
      endif
    enddo


    !--- j=1,2,3 are gammas, positron, antiproton---
    do j=1,3,1

      do i=1,PPPC4_Ebin,1
        x(1)=PPPC4_mx(lab(1))
        x(2)=PPPC4_mx(lab(2))
        y(1)=PPPC4_spect(j,lab(1),i)
        y(2)=PPPC4_spect(j,lab(2),i)
        call Between_two_pts(x(1),y(1),x(2),y(2),.true.,mx,ans)

        if (j==1) PPPC4_dnde(1,i)=PPPC4_K_over_mx(i)*mx
        PPPC4_dnde(j+1,i)=ans

        if (mx_lt_5GeV<5d0 .and. j==1) PPPC4_dnde(1,i)=PPPC4_K_over_mx(i)*mx_lt_5GeV
        if (mx_lt_5GeV<5d0) PPPC4_dnde(j+1,i)=ans*mx/mx_lt_5GeV

      enddo

    enddo

    return
  end subroutine




!'''
!  Try to find the Ecm for 2 particles with different masses.  
!  X, X > A1, B2
!'''
  subroutine Boost_Ecm_XX_to_AB(mx,m1,m2, Ecm1,Ecm2)
    implicit none
    real*8, intent(in) :: mx,m1,m2
    real*8, intent(out) :: Ecm1,Ecm2
    real*8 delta,Ecm

    Ecm=2d0*mx
    ! moment is the same at the cm frame
    delta=(m1**2-m2**2)/(2d0*Ecm)
    Ecm1=Ecm*0.5d0+delta
    Ecm2=Ecm*0.5d0-delta
    return
  end subroutine





!=======================================================================
end module PYTHIA_PPPC4
