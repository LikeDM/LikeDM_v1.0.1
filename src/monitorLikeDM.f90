
subroutine print_debug_info(mx,sv,aep,bep,aap,bap,chisq)
  use MathLib
  use ReadTable
  use charge_lepton
  use charge_antip
  use dSphs_gamma

  implicit none

  real*8, intent(in) ::  mx,sv,chisq(7)  
  real*8, intent(in) ::  aep(3),bep(3),aap,bap 
  integer i
  real*8 Fep_mod(100),Fap_mod(100)




  write(*,*)''
  write(*,*)'************************************************************'
  write(*,*)'--------- charged particle result >>>>>'
  if (sets__use_ep) then

      write(*,*)'AMS02efr:',chisq(3)
      write(*,*)'AMS02e+:',chisq(4)
      write(*,*)'AMS02e-:',chisq(5)
      write(*,*)'AMS02e+e-:',chisq(6)
      write(*,*)'AMS02_total_ep:',chisq(3)+chisq(4)+chisq(5)+chisq(6)
      write(*,*)''
  endif


  if (sets__use_ap) then
      write(*,*)'PAMELA_pbar:',chisq(7)
  endif
  write(*,*)'--------- charged particle result <<<<<'
  write(*,*)'************************************************************'
  write(*,*)''
  !=======================================================================
  !++++++++++ output debug info +++++++++++++++

  if (sets__seebug>0) then !++ 0 for nothing
    write(*,*)''
    if (sets__DMdecay .and. sets__seebug>0) then 
      write(*,'(A,1G12.4,A,1G12.4)')' mx (GeV) = ',mx,' tau (s) = ',1d0/sv
    else
      write(*,'(A,1G12.4,A,1G12.4)')' mx (GeV) = ',mx,' sigmav (cm^3 s^-1) =',sv
    endif

    if (sets__seebug==1) then !++ 1 for inputs
      write(*,*)'************************************************************'
      write(*,*)'seebug=',sets__seebug,':  LikeDM inputs >>>>>'
      write(*,*)'sets__DMdecay=',sets__DMdecay
      write(*,*)'sets__seebug=',sets__seebug
      write(*,*)'sets__halo=',sets__halo
      write(*,*)'sets__use_dSphs=',sets__use_dSphs
      write(*,*)'sets__use_ep=',sets__use_ep
      write(*,*)'sets__use_ap=',sets__use_ap
      !write(*,*)'sets__use_DD=',sets__use_DD
      write(*,*)'smod__ep :',smod__ep
      write(*,*)'smod__ap :',smod__ap
    endif

    if (sets__seebug==2) then !++ 2 for (alpha,beta)
      write(*,*)'************************************************************'
      write(*,*)'seebug=',sets__seebug,':  fitting results >>>>>'
      write(*,*)'--------- e+/e- >>>>>'
      write(*,'(A,3G12.4)')'alpha[1:3]=',aep(:)
      write(*,'(A,3G12.4)')'beta[1:3]=',bep(:)
      write(*,'(A,1E12.4)')' chisq/-2ln(likelihood)= ',chisq(3)+chisq(4)+chisq(5)+chisq(6)
    endif

    if (sets__seebug==2) then
      write(*,*)''
      write(*,*)'--------- antiproton >>>>>'
      write(*,'(A,G12.4,A,G12.4)')'alpha=',aap,' beta=',bap
      write(*,'(A,1E12.4)')' chisq/-2ln(likelihood)= ',chisq(7)
    endif


    if (sets__seebug==3) then !++ 3 for input dNdE
      write(*,*)'************************************************************'
      write(*,*)'seebug=',sets__seebug,':  dNdE >>>>>'
      write(*,*)' >>> E(GeV), dNdE(gamma), dNdE(e+), dNdE(pbar) <<< ' 
      write(*,*)'' 
      do i=1,PYTHIA__nrows,1
        write(*,'(4G12.4)') PYTHIA__dnde(1:4,i)
      enddo
    endif

 
    if (sets__seebug==4) then !++ 4 for propagated fluxes of e+ and pbar
      write(*,*)'************************************************************'
      write(*,*)'seebug=',sets__seebug,':  propagated spectra >>>>>'
      write(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      write(*,*)'  positron and antiproton flux AFTER propagation  '
      write(*,*)'              BEFORE solar Modulation             ' 
      write(*,*)'         En(e+), Phi(e+), En(ap), Phi(ap)         '
      write(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      do i=1,100,1
        write(*,'(4G12.4)') PreMod__DMpred(1,1:2,i),PreMod__DMpred(2,1:2,i)
      enddo

      write(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      write(*,*)'  positron and antiproton flux AFTER propagation  '
      write(*,*)'              AFTER solar Modulation              ' 
      write(*,'(A,G12.4,A,G12.4)')'   PHI_e^+=',smod__ep,' PHI_pbar=',smod__ap
      write(*,*)'         En(e+), Phi(e+), En(ap), Phi(ap)         '
      write(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      call AfterModul('e+',100,PreMod__DMpred(1,1,:),PreMod__DMpred(1,2,:),Fep_mod)
      call AfterModul('pbar',100,PreMod__DMpred(2,1,:),PreMod__DMpred(2,2,:),Fap_mod)
      do i=1,100,1
        write(*,'(4G12.4)') PreMod__DMpred(1,1,i),Fep_mod(i),PreMod__DMpred(2,1,i),Fap_mod(i)
      enddo
    endif

  endif


  return
end subroutine print_debug_info
