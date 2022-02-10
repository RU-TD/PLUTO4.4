module prizmo_tdust
  use prizmo_commons
  integer,parameter::ntdust=100
  real*8::xtdust_em(ntdust), ytdust_em(ntdust), rtdust_em(ntdust-1)
  real*8::ftdust_abs(nphoto), collision_factor
  real*8::xtdust_min, xtdust_fact
contains

  ! *******************
  ! get dust temperature, K
  function get_tdust(n, Tgas, jflux) result(Tdust)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas, jflux(nphoto)
    real*8::Tdust, Av, ngas, mu

    ! compute Tdust with bisection only if flag is set
    if(tdust_mode == 0) then
       ngas = get_ntot(n(:))
       mu = get_mu(n(:))

       ! compute dust temperature using bisection
       Tdust = get_tdust_bisect(Tgas, mu, ngas, jflux(:))

       if(Tdust < 0d0) then
          print *, "ERROR: negative Tdust (prizmo_tdust.f90)"
          stop
       end if
    else if(tdust_mode == 1) then
       ! constant dust temperature, K
       Tdust = 2d1
    else if(tdust_mode == 2) then
      ! use Hocuk+2017, eq.8, https://arxiv.org/pdf/1704.02763.pdf
      Av = min(max(variable_Av, 1e-2), 4e2)
      Tdust = (11. + 5.7 * tanh(0.61 - log10(Av))) * variable_G0**(1e0 / 5.9)
    else if(tdust_mode == 3) then
      Tdust = Tgas
    end if

  end function get_tdust

  ! *******************
  subroutine load_tdust_data()
    use prizmo_commons
    use prizmo_utils
    implicit none
    integer::unit_tdust, i
    real*8::alpha, mu, dummy

    ! load emission as a function of tdust (precomputed for bisection)
    print *, "loading dust emission precomputed data..."
    open(newunit=unit_tdust, file=runtime_folder//"dust_emission.dat", status="old")
    do i=1,ntdust
       read(unit_tdust, *) xtdust_em(i), ytdust_em(i)
    end do
    close(unit_tdust)

    xtdust_em(:) = log10(xtdust_em(:))
    ytdust_em(:) = log10(ytdust_em(:))

    ! compute ratio for later use in interpolation
    do i=1, ntdust-1
       rtdust_em(i) = (ytdust_em(i+1) - ytdust_em(i)) / (xtdust_em(i+1) - xtdust_em(i))
    end do

    xtdust_min = xtdust_em(1)
    xtdust_fact = (ntdust - 1) / (xtdust_em(ntdust) - xtdust_em(1))

    ! load absorption coefficients to be multiplied by the multi-frequency intensity
    print *, "loading dust absorption precomputed data..."
    open(newunit=unit_tdust, file=runtime_folder//"dust_absorption.dat", status="old")
    do i=1,nphoto
       read(unit_tdust, *) dummy, ftdust_abs(i)
    end do
    close(unit_tdust)

    alpha = 0.5  ! thermal exchange efficiency factor
    mu = 2.34  ! mean molecular weight

    ! compute dust collisional pre-factor, eV/s/cm3
    ! NOTE: first kboltzmann in erg, second in eV for matching emission and absorption units
    collision_factor = 2d0 * sqrt(8d0 * kboltzmann / pi / proton_mass) &
         * kboltzmann_eV * pi * dust_fact3 * alpha

  end subroutine load_tdust_data

  ! *******************
  ! compute absorption term for bisection
  function get_absorption_term(jflux) result(af)
    use prizmo_commons
    implicit none
    real*8,intent(in)::jflux(nphoto)
    real*8::af

    af = 2d0 * sum(jflux(:) * ftdust_abs(:))

  end function get_absorption_term

  ! *******************
  ! get dust temperature using bisection
  function get_tdust_bisect(Tgas, mu, ngas, jflux) result(x)
    use prizmo_commons
    implicit none
    real*8,intent(in)::Tgas, ngas, jflux(nphoto), mu
    real*8::af, xmin, xmax, x, f, fmin, fmax
    integer::i, istep, unit

    xmin = 1d0 ! initial left temperature
    xmax = 3.99999d3 ! initial right temperature

    ! absorption term
    af = get_absorption_term(jflux(:))

    ! compute initial points
    fmin = f_bisect(xmin, Tgas, mu,  ngas, af)
    fmax = f_bisect(xmax, Tgas, mu, ngas, af)

    ! check if bisection has roots
    if(fmin * fmax > 0d0) then

       print *, "ERROR: can't find root in tdust bisection!"
       print *, xmin, fmin
       print *, xmax, fmax

       ! write function to find roots in the interval in debug file
       istep = 100
       open(newunit=unit, file="bisection_debug.dat", status="replace")
       do i=1,istep
          x = 1e1**((log10(xmax) - log10(xmin)) * (i - 1) / (istep - 1) + log10(xmin))
          write(unit, '(99E17.8e3)') x, f_bisect(x, Tgas, mu, ngas, af), &
               femission(x), collision_factor * sqrt(Tgas / mu) * (Tgas - x) * ngas, af
       end do
       close(unit)
       print *, "debug in bisection_debug.dat"
       print *, "file format is: Tdust, f_bisect, f_emission, f_thermal, f_absorption"
       stop
    end if

    ! loop for bisection
    do
       ! average point
       x = (xmin + xmax) / 2d0
       ! evaluate function at average point
       f = f_bisect(x, Tgas, mu, ngas, af)
       ! determine new position according to signs
       if(f * fmax > 0d0) then
          xmax = x
       else
          xmin = x
       end if
       ! loop until some accuracy is found
       if (abs(xmin - xmax) < 1d-12) exit
    end do

  end function get_tdust_bisect

  ! *******************
  function f_bisect(x, Tgas, mu, ngas, af) result(f)
    use prizmo_commons
    implicit none
    real*8,intent(in)::x, ngas, Tgas, af, mu
    real*8::f

    ! compute function for bisection
    f = femission(x) - collision_factor * sqrt(Tgas / mu) * (Tgas - x) * ngas - af

  end function f_bisect

  ! *****************
  ! compute emission function from fit
  function femission(xin) result(f)
    implicit none
    real*8,intent(in)::xin
    real*8::f, x1, f1, x
    integer::idx

    x = log10(xin)

    idx = floor((x - xtdust_min) * xtdust_fact + 1)

    x1 = xtdust_em(idx)
    f1 = ytdust_em(idx)

    f = rtdust_em(idx) * (x - x1) + f1
    f = 1d1**f

  end function femission

end module prizmo_tdust
