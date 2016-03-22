!###############################################################
!#
!# Note: Module of Math lib 
!#
!# Author: Yue-Lin Sming Tsai 
!# Email: smingtsai@gmail.com                                                                                   
!# Date: 2012-11-15
!# Date: 2013-10-21
!###############################################################


module MathLib

  implicit none
      
contains



!###########################################################
!#  Simple interpolation between two points (array)
!###########################################################
subroutine linInt2p(X,Y,n,ORlog,x_in,ans,DoSort_,Extropl_)                                                                          
  implicit none
  integer, intent(in) :: n
  logical, intent(in) :: ORlog
  real*8, intent(in) :: X(n),Y(n),x_in
  real*8, intent(out) :: ans
  logical,optional,intent(in) :: DoSort_,extropl_ 
  logical :: DoSort,Extropl
  real*8 x0,y0,x1,y1,xin  
  real*8 P(2,n),A(2,n)                                                                 
  integer i                                                                              

  xin=x_in

  if (present(DoSort_)) then
    DoSort=DoSort_
  else
     DoSort=.true.   
  endif

  if (present(Extropl_)) then
    Extropl=Extropl_
  else
     Extropl=.true.   
  endif



  if (DoSort) then
    !first to sort the X from 
    P(1,:)=X(:)
    P(2,:)=Y(:)
    call sort_small_big(2,n,P,A)
  else
    A(1,:)=X(:)
    A(2,:)=Y(:)
  endif


  if (A(1,1)>=xin) then 
    if (Extropl) then
        x0=A(1,1)
        y0=A(2,1)
        x1=A(1,2)
        y1=A(2,2)
    else
        ans=A(2,1)
        return
    endif


  elseif (A(1,n)<=xin) then

    if (Extropl) then
        x0=A(1,n)
        y0=A(2,n)
        x1=A(1,n-1)
        y1=A(2,n-1)
    else
        ans=A(2,n)
        return
    endif



  else
    do i=1,n-1,1
      if (A(1,i)<=xin .and. A(1,i+1)>xin ) then
        x0=A(1,i)
        y0=A(2,i)
        x1=A(1,i+1)
        y1=A(2,i+1)
        exit
      endif

    enddo
  endif

  call Between_two_pts(x0,y0,x1,y1,ORlog,xin,ans)

  return
end subroutine


!###########################################################
!#  Simple interpolation between two points (giving two point)
!###########################################################

subroutine Between_two_pts(a0,b0,a1,b1,ORlog,ain,ans)                                                                          
  implicit none                                                                              
  real*8 x0,y0,x1,y1,xin  
  real*8, intent(in) :: a0,b0,a1,b1,ain                                                                 
  logical, intent(in) :: ORlog
  logical ORlog_
  real*8, intent(out) :: ans
  logical use_absx,use_absy

  use_absx=.false.
  use_absy=.false.
  ORlog_=ORlog


  if (ORlog) then

    if (a0==0.0d0) then 
      x0=1.0d-80
    else
      x0=a0
    endif

    if (b0==0.0d0) then 
      y0=1.0d-80
    else
      y0=b0
    endif

    if (a1==0.0d0) then 
      x1=1.0d-80
    else
      x1=a1
    endif

    if (b1==0.0d0) then 
      y1=1.0d-80
    else
      y1=b1
    endif

    ! to decided whether is to use abs(x) and abs(y) in interpolation  
    if (a0<0d0 .and. a1<0d0 .and. ain<0d0) then 
      use_absx=.true. 
      x0=dabs(x0)
      x1=dabs(x1)
      xin=dabs(ain) 
    endif

    if (b0<0d0 .and. b1<0d0) then 
      use_absy=.true. 
      y0=dabs(y0)
      y1=dabs(y1)
    endif

    if (a0*a1<0d0 .or. b0*b1<0d0 ) then 
      xin= ain
      x0= a0
      y0= b0
      x1= a1
      y1= b1
      ORlog_=.false. 
    else
      xin= dlog10(ain)
      x0= dlog10(x0)
      y0= dlog10(y0)
      x1= dlog10(x1)
      y1= dlog10(y1)

    endif


  else
    xin= ain
    x0= a0
    y0= b0
    x1= a1
    y1= b1

  endif

  ans=((xin-x0)/(x1-x0)) *  (y1-y0) + y0


  if (ORlog_) then
        ans=10.0d0**(ans)
  endif

  if (use_absy) ans=-ans

  return
end subroutine






!#########################################################################
!#
!#  Subroutine to order according to first cols
!#
!#########################################################################

!%%% sort_chisq  %%%
subroutine sort_small_big(ncols,nrows,P,A)
  implicit none

  integer, intent(in) :: ncols,nrows ! number variables, number of nrows 
  real*8, intent(in) :: P(ncols,nrows) !input
  real*8, intent(out) :: A(ncols,nrows) !output
  real*8 temp
  integer i,ione(nrows),k,smallest,j
  real*8 P_(ncols,nrows) !input


  P_(:,:)=P(:,:)
  ! declare labels of lines
  do i=1,nrows,1
    ione(i)=i
  enddo


  do i = 1,nrows-1,1
     smallest = i
     do j=i+1,nrows,1
       if (P_(1,j) <= P_(1,smallest)) then
         Smallest = j
       endif
     enddo !j 

     !change A[I] and A[smallest]
     temp=P_(1,i)
     P_(1,i)=P_(1,smallest)
     P_(1,smallest)=temp

     k=ione(i)
     ione(i)=ione(smallest)
     ione(smallest)=k

  enddo !finished ordering according to chisq

  !write out
  do i=1,nrows,1  

    !check the errors
    if (P_(1,i)>P_(1,i+1) .and. i<nrows) then
      write(*,*)'wrong?'
      STOP
    endif

    ! j is the address according to chisq
    j=ione(i)
    A(1,i)=P_(1,i)
    A(2:ncols,i)=P_(2:ncols,j)

  enddo   ! i

  return
end subroutine







!###########################################################
!#  Simple 2D interpolation by reusing 1D subroutine 
!###########################################################
  subroutine interpXYtoZ(X,Y,Z,nx,ny,use_log,xo,yo,ans)
    implicit none
    integer, intent(in) :: nx,ny                                                                             
    real*8, intent(in) :: X(nx),Y(ny),Z(nx,ny)
    logical, intent(in) :: use_log
    real*8, intent(in) :: xo,yo
    real*8, intent(out) :: ans
    real*8 z1,z2
    integer i,i1,i2 
    i1=0
    i2=0
    do i=1,nx-1,1
      if (X(i)<=xo .and. X(i+1)>xo ) then
        i1=i
        i2=i+1
        exit
      endif    
    enddo 

    if (i1==0 .or. i2==0) then 
      write(*,*)'interpXYtoZ: your table has problem to locate xo'
      STOP
    endif

    call linInt2p(Y,Z(i1,:),ny,use_log,yo,z1)
    call linInt2p(Y,Z(i2,:),ny,use_log,yo,z2)
 
    call Between_two_pts(X(i1),z1,X(i2),z2,use_log,xo,ans)

    return
  end subroutine






!###########################################################
!#  Compute a step size in the log scale
!###########################################################
subroutine factorCalc(facor,uplim,lowlim,n)                                                                          
  implicit none 
  integer, intent(in) :: n                                                                             
  real*8, intent(in) :: uplim,lowlim
  real*8, intent(out) :: facor                                                                    
  real*8 logx0

  logx0=dlog10(uplim/lowlim)/dble(n-1)
  facor =10.0d0**(logx0)
  return
end subroutine




!#########################################################################
!#
!#  To count how many rows in a given files
!#
!#########################################################################
subroutine count_nrows(file_name,nrows)
  implicit none
  character (LEN = 80), intent(in) :: file_name
  integer, intent(out) :: nrows
  integer i
  real*8 aa
  open(unit=32,file=file_name,status='old')
  i=1
  do while(.true.)
    read(unit=32,fmt=*,end=117) aa
    i = i + 1
  enddo 
117   close(32)

  nrows=i
  return
end subroutine 




  subroutine TableArea(X,Y,n,area)
    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: X(n),Y(n)
    real*8, intent(out) ::  area
    real*8 ans
    integer i
    area=0d0
    do i=1,n-1,1
      ans=(Y(i)+Y(i+1))*dabs(X(i+1)-X(i))*0.5d0
      area=area+ans
    enddo 

    return
  end subroutine 






!##################################################################
  !--- Compute Likelihood ---
  subroutine PoissonProb(obs,bkg,s,ans)
!##################################################################
   implicit none
   real*8, intent(in) :: obs,bkg,s ! input
   real*8, intent(out) :: ans ! Output
   !=========================
   integer i,o
   real*8 sumlogx,lnPoisson

   ! ans= (s+b)^o exp(-(s+b)) / o!
   ! In order to save the diverge problem, we use log...  
   ! 1) o!->ln(o!)
   o=dint(obs)
   sumlogx=0.0d0
   do i=1,o,1
    sumlogx=sumlogx+dlog(dble(i))
   enddo

   ! 2) ln Poisson = O*ln(s+b) -(s+b) -ln(o!)
   lnPoisson= obs*dlog(bkg+s) - (bkg+s) - sumlogx

   ! 3) Poisson=exp (lnPoisson)
   ans=dexp(lnPoisson)

   return
  end subroutine 



!##################################################################
 subroutine GauProb(cv,theo,sigma,tau,prob)
!##################################################################
   implicit none
   real*8, intent(in) :: cv,theo,sigma,tau
   real*8, intent(out) :: prob
   real*8 errsq
   real*8 chisq
   errsq= sigma**2 + tau**2
   chisq=(cv-theo)**2 /errsq
   prob=dexp(-0.5d0*chisq)
   return
 end subroutine


!##################################################################
  subroutine Gau_prior(f0,f1,cv,width,ans)
!##################################################################
    implicit none
    real*8, intent(in) :: f0,f1,cv,width
    real*8, intent(out) :: ans
    real*8 chi
    chi=dsqrt(dlog(f0)*(-2d0))*dcos( f1* dasin(1.0d0)*4d0)
    ans=cv+chi*width

    return
  end subroutine


  subroutine Gen_linspace(vmin,vmax,num,use_log,ans)
    implicit none
    integer, intent(in) :: num       ! input             
    logical, intent(in) :: use_log   ! input
    real*8, intent(in) :: vmin,vmax  ! input
    real*8, intent(out) :: ans(num)   ! output
    !--- dummy varible ---
    real*8 fac
    integer i
    if (use_log) then 
      call factorCalc(fac,vmax,vmin,num)
    else
      fac=(vmax-vmin)/dble(num-1)
    endif
    ans(1)=vmin
    do i=2,num,1
      if (use_log) then
        ans(i)=ans(i-1)*fac
      else
        ans(i)=ans(i-1)+fac
      endif
    enddo


    return
  end subroutine 





  !-------------------------------------------------------------------
  !# Numerov algorithm to solve y"+f(x)*y=0
  !-------------------------------------------------------------------
  subroutine numerov(n,x,f,y1,y2,y,yp)
    implicit none               
    integer, intent(in) :: n    ! total n grid points 
    ! y"+f(x)*y=0
    real*8, intent(in) :: x(n),y1,y2
    real*8, intent(in) :: f(n)  ! tabled in n grid 
    real*8, intent(out) :: y(n),yp(n)
    integer i
    real*8 h
    real*8 b

    h=x(n)-x(n-1)
    b=(h**2)/12d0                         

    y(1)=y1
    y(2)=y2
    do i=2,n-1,1  
      y(i+1) = 2d0 *y(i)*(1d0-5d0*b*f(i)) 

      y(i+1) = y(i+1)- (1d0 +b*f(i-1))*y(i-1)
      y(i+1) = y(i+1)/(1d0+b*f(i+1))

    enddo


    !--- Generally speaking the symmetric derivatives are preferred over 
    !    the forward or backward derivatives.  
    !    Another useful approximation to the derivative is the "5 point formula"
    
    do i=3,n-2,1      
      yp(i)=y(i-2)-8d0*y(i-1)+8*y(i+1)-y(i+2)
      yp(i)=yp(i)/12d0/h
    enddo
    !--- extroplate to i=1,2,n-1,n ---
    call Between_two_pts(x(3),yp(3),x(4),yp(4),.false.,x(1),yp(1))
    call Between_two_pts(x(3),yp(3),x(4),yp(4),.false.,x(2),yp(2)) 

    call Between_two_pts(x(n-3),yp(n-3),x(n-2),yp(n-2),.false.,x(n-1),yp(n-1))
    call Between_two_pts(x(n-3),yp(n-3),x(n-2),yp(n-2),.false.,x(n),yp(n))

    return

  end subroutine








subroutine Jacobi(ain,xin,n,abserr,a,x)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: ain(n,n),xin(n,n)
  double precision, intent(out) :: a(n,n),x(n,n)
  integer i, j, k
  double precision abserr, b2, bar
  double precision beta, coeff, c, s, cs, sc

  a(:,:)=ain(:,:)
  x(:,:)=xin(:,:)
! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
  x = 0.0
  do i=1,n
    x(i,i) = 1.0
  enddo

! find the sum of all off-diagonal elements (squared)
  b2 = 0.0
  do i=1,n
    do j=1,n
      if (i.ne.j) b2 = b2 + a(i,j)**2
    end do
  end do

  if (b2 <= abserr) return

  ! average for off-diagonal elements /2
  bar = 0.5*b2/float(n*n)

  do while (b2.gt.abserr)
    do i=1,n-1
      do j=i+1,n
        if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
        b2 = b2 - 2.0*a(j,i)**2
        bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
        beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
        coeff = 0.5*beta/sqrt(1.0+beta**2)
        s = dsqrt(max(0.5+coeff,0.0))
        c = dsqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
        do k=1,n
          cs =  c*a(i,k)+s*a(j,k)
          sc = -s*a(i,k)+c*a(j,k)
          a(i,k) = cs
          a(j,k) = sc
        end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
        do k=1,n
          cs =  c*a(k,i)+s*a(k,j)
          sc = -s*a(k,i)+c*a(k,j)
          a(k,i) = cs
          a(k,j) = sc
          cs =  c*x(k,i)+s*x(k,j)
          sc = -s*x(k,i)+c*x(k,j)
          x(k,i) = cs
          x(k,j) = sc
        enddo
      enddo
    enddo
  enddo
  return
end subroutine






!=======================================================================
end module MathLib













