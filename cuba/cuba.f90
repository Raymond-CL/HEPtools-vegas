module shared
  implicit none
  double precision :: rad
end module shared

program main
  use shared
  implicit none
  integer :: ndim,ncomp,integrand,userdata,nvec
  external :: integrand
  double precision :: epsrel,epsabs
  integer :: flags,seed,mineval,maxeval
  integer :: nstart,nincrease,nbatch,gridno
  character(len=80) :: statefile
  integer :: spin
  integer :: neval,fail
  double precision :: integral(1),error(1),prob(1)
  double precision, parameter :: PI = 4d0*atan(1d0)  

  rad = 5d0

  ndim = 3
  ncomp = 1
  userdata = 0
  nvec = 1
  epsrel = 1d-3
  epsabs = 0d0
  flags = 0
  seed = 0
  mineval = 10000
  maxeval = 200000000
  nstart = 1000
  nincrease = 1000
  nbatch = 1000
  gridno = 0
  statefile = ''
  spin = 0

  call vegas(ndim, ncomp, integrand, userdata, nvec, &
    epsrel, epsabs, flags, seed, mineval, maxeval, &
    nstart, nincrease, nbatch, gridno, statefile, spin, &
    neval, fail, integral, error, prob)

  write(*,*) 'neval = ', neval
  write(*,*) 'fail = ', fail
  write(*,*) 'integral = ', integral
  write(*,*) 'actual   = ', PI**(ndim/2d0) / gamma(ndim/2d0+1d0) * rad**ndim
  write(*,*) 'error = ', error
  write(*,*) 'prob = ', prob

end program main



integer function integrand(ndim, x, ncomp, f)
  use shared
  implicit none
  integer, intent(in) :: ndim, ncomp
  double precision, intent(in) :: x(ndim)
  double precision, intent(out) :: f(ncomp)
  integer :: i
  double precision :: dx(ndim),upper(ndim),lower(ndim)
  double precision :: jacobian
  double precision :: r

  upper = +rad
  lower = -rad

  jacobian = 1d0
  do i=1,ndim
    dx(i) = lower(i) + (upper(i) - lower(i)) * x(i)
    jacobian = jacobian * (upper(i) - lower(i))
  enddo

  r = sqrt(sum(dx(:)**2))

  if(r < rad) then
    f(1) = 1d0 * jacobian
  else
    f(1) = 0d0
  endif

  integrand = 0

end function integrand
