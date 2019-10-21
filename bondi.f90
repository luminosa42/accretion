! "Analytic" solution for Bondi accretion
! PR 8/14/13

program bondi

!==============================================================================

implicit none

integer, parameter :: N = 1000  !number of data points
real, parameter    :: l = 3.  

integer :: i
real    :: dx, logx
real	:: lambda, x, alpha, v, gamma, rs, xmin, xmax, vmin, vmax
real    :: alpha_old, v_old, dlnadlnx, dlnvdlnx, zbrent, f, H, g
real    :: fminloc, xsonic
external :: f, g

common gamma, lambda, x

!==============================================================================
dx = 2*l/N

! If given a single gamma
!gamma = 5./3.
gamma = 1.1
if (gamma == 1.0) then
  lambda = 0.25*exp(1.5)
else if (gamma == 5./3.) then
  lambda = 0.25
else
  lambda = 0.5**((gamma+1.)/(2.*(gamma-1.))) * &
           ((5.-3.*gamma)/4.)**(-(5.-3.*gamma)/(2.*(gamma-1.)))
endif
!write(*,'(A,F13.7)') '# lambda_c = ', lambda

rs = xsonic()
!write(*,'(A,F13.7)') '# x_sonic = ', rs

!logarithmic distribution
x = 10.**((N+1)*dx/2. )
vmin = 1.e-10
vmax = fminloc()
v_old = zbrent(f, vmin, vmax, 1.e-6)
alpha_old = lambda/(x*x*v_old)

do i = N, 1, -1
  logx = -l + i*dx
  x = 10.**logx
  !print *, logx
  vmin = 1.e-10
  vmax = fminloc()
  v = zbrent(f, vmin, vmax, 1.e-6)
  print *, v 
  if (v < v_old) then
    vmin = fminloc()
    vmax = 1.e5
    v = zbrent(f, vmin, vmax, 1.e-6)
  endif
  alpha = lambda/(x*x*v)
  dlnadlnx = alog(alpha_old/alpha) / alog((x+dx)/x)
  dlnvdlnx = alog(v_old/v) / alog((x+dx)/x)
  alpha_old = alpha
  v_old = v
  write(*,'(5(ES13.6,:,2X))') x, alpha, v, dlnadlnx, dlnvdlnx
enddo

!==============================================================================

end program bondi

!==============================================================================

real function xsonic()
implicit none
real :: lambda, gamma
real :: xsonic
common gamma, lambda

if (gamma == 1.0) then
  xsonic = 0.5
else
  !xsonic = (1./(2.*lambda**(gamma-1.)))**(1./(3.-2.*gamma))
  xsonic = (5. - 3. * gamma)/4.
endif

return
end function xsonic

real function f(v)
implicit none
real :: v, x, lambda, gamma
real :: H
common gamma, lambda, x
f = 0.5*v**2. + H(lambda/(x*x*v)) - 1./x
return
end function f

real function fminloc()
implicit none
real :: x, lambda, gamma
common gamma, lambda, x
if (gamma == 1.0) then
  fminloc = 1.
else
  fminloc = (lambda/x**2)**((gamma-1.)/(gamma+1.))
endif
return
end function fminloc


real function g(x)
implicit none
real :: x, lambda, gamma
real :: H
common gamma, lambda
g = 0.5 + H(lambda/x**2) - 1./x
return
end function g


real function H(alpha)
implicit none
real :: alpha, gamma, lambda, x
common gamma, lambda, x

if (gamma == 1.0) then
  H = log(alpha)
else
  H = 1./(gamma-1.)*(alpha**(gamma-1.) - 1.)
endif

return
end function H



      FUNCTION zbrent(func,x1,x2,tol)
      INTEGER ITMAX
      REAL zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))stop &
     &'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      stop 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
