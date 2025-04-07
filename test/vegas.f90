module vint

  ! VEGAS multi-dimensional numerical integration module in 64-bit Fortran
  ! ------------------------------------------------------------------------
  ! VEGAS subroutine is inspired from:
  ! - original code by P. Lepage [J. Comput. Phys. 27 (1978) 192]
  ! - VEGAS subroutine in Numerical Recipe (NR) by Press et.al. (not open source!)
  ! - VEGAS subroutine in GNU Scientific Library (GSL)
  ! - VEGAS subroutine in Cuba library by T. Hahn [hep-ph/0404043]
  ! VEGAS uses Random Number Generator (RNG):
  ! - legacy code uses RANDA from cernlib
  ! - NR uses ran1 (low period)
  ! - this version uses Mersenne twister mt19937 64-bit code modified from:
  !   Remi Piatek, original from Takuji Nishimura and Makoto Matsumoto

  use, intrinsic :: iso_fortran_env, only : r8=>real64, i8=>int64
  implicit none
  private
  ! RNG parameters (do not alter)
  integer(i8), parameter :: nn       = 312_i8
  integer(i8), parameter :: mm       = 156_i8
  integer(i8), parameter :: seed_def = 5489_i8
  integer(i8), parameter :: matrix_a = -5403634167711393303_i8
  integer(i8), parameter :: um       = -2147483648_i8 ! most significant 33 bits
  integer(i8), parameter :: lm       =  2147483647_i8 ! least significant 31 bits
  real(r8),    parameter :: pi252    = 1._r8/(2._r8**52)
  ! RNG saved states
  integer(i8) :: mt(nn)       ! array for the state vector
  integer     :: mti = nn+1   ! mti==nn+1 means mt(nn) is not initialized
  ! public subroutine
  public :: vegas!, init_genrand64, genrand64_real ! (uncomment if you want to use the RNG)

contains

  ! initialize RNG with a seed (optional)
  subroutine init_genrand64(seed)
    integer(i8), intent(in) :: seed
    integer :: i
    mt(1) = seed
    do i = 1, nn-1
      mt(i+1) = 6364136223846793005_i8 * ieor(mt(i), ishft(mt(i), -62)) + i
    end do
    mti = nn
  end subroutine init_genrand64

  ! integer RNG
  integer(r8) function genrand64_int64()
    implicit none
    integer(i8) :: mag01(0:1) = (/0_i8, matrix_a/)
    integer(i8) :: x
    integer     :: i
    if(mti >= nn) then ! generate nn words at one time
      ! if init_genrand64() has not been called, a default initial seed is used
      if(mti == nn+1) call init_genrand64(seed_def)
      do i = 1, nn-mm
        x = ior(iand(mt(i),um), iand(mt(i+1), lm))
        mt(i) = ieor(ieor(mt(i+mm), ishft(x, -1)), mag01(iand(x, 1_i8)))
      end do
      do i = nn-mm+1, nn-1
        x = ior(iand(mt(i), um), iand(mt(i+1), lm))
        mt(i) = ieor(ieor(mt(i+mm-nn), ishft(x, -1)), mag01(iand(x, 1_i8)))
      end do
      x = ior(iand(mt(nn), um), iand(mt(1), lm))
      mt(nn) = ieor(ieor(mt(mm), ishft(x, -1)), mag01(iand(x, 1_i8)))
      mti = 0
    end if
    mti = mti + 1
    x = mt(mti)
    x = ieor(x, iand(ishft(x,-29), 6148914691236517205_i8))
    x = ieor(x, iand(ishft(x, 17), 8202884508482404352_i8))
    x = ieor(x, iand(ishft(x, 37),   -2270628950310912_i8))
    x = ieor(x, ishft(x, -43))
    genrand64_int64 = x
  end function genrand64_int64

  ! real RNG used in vegas
  real(r8) function genrand64_real()
    genrand64_real = real(ishft(genrand64_int64(), -12), kind=r8)
    genrand64_real = (genrand64_real + 0.5_r8) * pi252
  end function genrand64_real

  ! VEGAS subroutine
  subroutine vegas(region,func,init,ncall,itmx,nprn,tgral,sd,chi2a)
    real(r8), dimension(:), intent(in) :: region
    interface
      double precision function func(pt,wgt)
        double precision, intent(in) :: pt(:), wgt
      end function func
    end interface
    integer(i8), intent(in) :: init,ncall,itmx,nprn
    real(r8), intent(out) :: tgral,sd,chi2a
    real(r8), parameter :: alph=1.5_r8,tiny=1.0e-30_r8
    integer(i8), parameter :: mxdim=20,ndmx=50
    integer(i8), save :: i,it,j,k,mds,nd,ndim,ndo,ng,npg
    integer(i8), dimension(mxdim), save :: ia,kg
    real(r8), save :: calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo,harvest
    real(r8), dimension(ndmx,mxdim), save :: d,di,xi
    real(r8), dimension(mxdim), save :: dt,dx,x
    real(r8), dimension(ndmx), save :: r,xin
    real(r8), save :: schi,si,swgt
    ndim=size(region)/2
    if (init <= 0) then
      mds=1
      ndo=1
      xi(1,:)=1.0_r8
    end if
    if (init <= 1) then
      si=0.0_r8
      swgt=0.0_r8
      schi=0.0_r8
    end if
    if (init <= 2) then
      nd=ndmx
      ng=1
      if (mds /= 0) then
        ng=int((ncall/2.0_r8+0.25_r8)**(1.0_r8/ndim))
        mds=1
        if ((2*ng-ndmx) >= 0) then
          mds=-1
          npg=ng/ndmx+1
          nd=ng/npg
          ng=npg*nd
        end if
      end if
      k=ng**ndim
      npg=max(ncall/k,2)
      calls=real(npg,r8)*real(k,r8)
      dxg=1.0_r8/ng
      dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1.0_r8)
      xnd=nd
      dxg=dxg*xnd
      dx(1:ndim)=region(1+ndim:2*ndim)-region(1:ndim)
      xjac=1.0_r8/calls*product(dx(1:ndim))
      if (nd /= ndo) then
        r(1:max(nd,ndo))=1.0
        do j=1,ndim
          call rebin(ndo/xnd,nd,r,xin,xi(:,j))
        end do
        ndo=nd
      end if
      if (nprn >= 0) write(*,200) ndim,calls,it,itmx,nprn,&
        alph,mds,nd,(j,region(j),j,region(j+ndim),j=1,ndim)
    end if
    do it=1,itmx
      ti=0.0
      tsi=0.0
      kg(:)=1
      d(1:nd,:)=0.0
      di(1:nd,:)=0.0
      iterate: do
        fb=0.0
        f2b=0.0
        do k=1,npg
          wgt=xjac
          do j=1,ndim
            harvest = genrand64_real() ! get random number from RNG
            xn=(kg(j)-harvest)*dxg+1.0_r8
            ia(j)=max(min(int(xn),ndmx),1)
            if (ia(j) > 1) then
              xo=xi(ia(j),j)-xi(ia(j)-1,j)
              rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
            else
              xo=xi(ia(j),j)
              rc=(xn-ia(j))*xo
            end if
            x(j)=region(j)+rc*dx(j)
            wgt=wgt*xo*xnd
          end do
          f=wgt*func(x(1:ndim),wgt)
          f2=f*f
          fb=fb+f
          f2b=f2b+f2
          do j=1,ndim
            di(ia(j),j)=di(ia(j),j)+f
            if (mds >= 0) d(ia(j),j)=d(ia(j),j)+f2
          end do
        end do
        f2b=sqrt(f2b*npg)
        f2b=(f2b-fb)*(f2b+fb)
        if (f2b <= 0.0) f2b=tiny
        ti=ti+fb
        tsi=tsi+f2b
        if (mds < 0) then
          do j=1,ndim
            d(ia(j),j)=d(ia(j),j)+f2b
          end do
        end if
        do k=ndim,1,-1
          kg(k)=mod(kg(k),ng)+1
          if (kg(k) /= 1) cycle iterate
        end do
        exit iterate
      end do iterate
      tsi=tsi*dv2g
      wgt=1.0_r8/tsi
      si=si+real(wgt,r8)*real(ti,r8)
      schi=schi+real(wgt,r8)*real(ti,r8)**2
      swgt=swgt+real(wgt,r8)
      tgral=si/swgt
      chi2a=max((schi-si*tgral)/(it-0.99_r8),0.0_r8)
      sd=sqrt(1.0_r8/swgt)
      tsi=sqrt(tsi)
      if (nprn >= 0) then
        write(*,201) it,ti,tsi,tgral,sd,chi2a
        if (nprn /= 0) then
          do j=1,ndim
            write(*,202) j,(xi(i,j),di(i,j),&
              i=1+nprn/2,nd,nprn)
          end do
        end if
      end if
      do j=1,ndim
        xo=d(1,j)
        xn=d(2,j)
        d(1,j)=(xo+xn)/2.0_r8
        dt(j)=d(1,j)
        do i=2,nd-1
          rc=xo+xn
          xo=xn
          xn=d(i+1,j)
          d(i,j)=(rc+xn)/3.0_r8
          dt(j)=dt(j)+d(i,j)
        end do
        d(nd,j)=(xo+xn)/2.0_r8
        dt(j)=dt(j)+d(nd,j)
      end do
      where (d(1:nd,:) < tiny) d(1:nd,:)=tiny
      do j=1,ndim
        r(1:nd)=((1.0_r8-d(1:nd,j)/dt(j))/(log(dt(j))-log(d(1:nd,j))))**alph
        rc=sum(r(1:nd))
        call rebin(rc/xnd,nd,r,xin,xi(:,j))
      end do
    end do
200 format(/' input parameters for vegas:  ndim=',i3,'  ncall=',f8.0&
      /28x,'  it=',i5,'  itmx=',i5&
      /28x,'  nprn=',i3,'  alph=',f5.2/28x,'  mds=',i3,'   nd=',i4&
      /(30x,'xl(',i2,')= ',g11.4,' xu(',i2,')= ',g11.4))
201 format(/' iteration no.',i3,': ','integral =',g14.7,' +/- ',g9.2,&
      /' all iterations:   integral =',g14.7,' +/- ',g9.2,&
      ' chi**2/it''n =',g9.2)
202 format(/' data for axis ',i2/'    x       delta i       ',&
      '   x       delta i       ','    x       delta i       ',&
      /(1x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
  end subroutine vegas

  ! rebin subroutine for VEGAS
  subroutine rebin(rc,nd,r,xin,xi)
    real(r8), intent(in) :: rc
    integer(i8), intent(in) :: nd
    real(r8), dimension(:), intent(in) :: r
    real(r8), dimension(:), intent(out) :: xin
    real(r8), dimension(:), intent(inout) :: xi
    integer(i8) :: i,k
    real(r8) :: dr,xn,xo
    k=0_r8
    xo=0.0_r8
    dr=0.0_r8
    do i=1,nd-1
      do
        if (rc <= dr) exit
        k=k+1
        dr=dr+r(k)
      end do
      if (k > 1) xo=xi(k-1)
      xn=xi(k)
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr/r(k)
    end do
    xi(1:nd-1)=xin(1:nd-1)
    xi(nd)=1.0_r8
  end subroutine rebin

end module vint
