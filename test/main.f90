module shared
  implicit none
  double precision :: rad
end module shared



program main
  use, intrinsic :: iso_fortran_env, only : wp=>real64, ip=>int64
  use vint
  use shared
  implicit none
  integer, parameter :: ndim = 18
  real(wp), dimension(2*ndim) :: region
  integer(ip) :: init,ncall,itmx,nprn
  real(wp) :: tgral,sd,chi2a
  interface
    double precision function func(pt,wgt)
      implicit none
      double precision, intent(in) :: pt(:), wgt
    end function func
  end interface
  double precision, parameter :: PI = 4d0*atan(1d0)
  rad = 1d0
  ncall = 20000000_ip
  itmx = 1
  nprn = -1
  region(1:ndim) = -rad
  region(1+ndim:2*ndim) = +rad
  init = -1
  call vegas(region,func,init,ncall,itmx*10,nprn,tgral,sd,chi2a)
  init = +1
  call vegas(region,func,init,ncall*10,itmx,nprn,tgral,sd,chi2a)
  write(*,*) "result: ",tgral
  write(*,*) "volume: ",PI**(ndim*0.5d0) / gamma(ndim*0.5d0+1d0) * rad**ndim
  write(*,*) "stddev: ",sd
end program main



double precision function func(pt,wgt)
  use shared
  implicit none
  double precision, intent(in) :: pt(:), wgt
  double precision :: temp
  temp = wgt
  func = 0d0
  if(sum(pt(:)**2)<=rad*rad) func = 1d0
end function func
