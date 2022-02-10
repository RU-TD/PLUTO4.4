module prizmo_photo
contains

  ! *************************
  function get_draine_flux() result(jflux)
    use prizmo_commons, only: nphoto, energy_grid
    implicit none
    real*8::jflux(nphoto)
    integer::i

    ! loop on photo bins
    do i=1,nphoto
      jflux(i) = f_draine(energy_grid(i))
    end do

  end function get_draine_flux

  ! *************************
  function get_black_body_flux(Tbb) result(jflux)
    use prizmo_commons, only: nphoto, energy_grid
    implicit none
    real*8,intent(in)::Tbb
    real*8::jflux(nphoto)
    integer::i

    ! loop on photo bins
    do i=1,nphoto
      jflux(i) = f_black_body(energy_grid(i), Tbb)
    end do

  end function get_black_body_flux

  ! *************************
  function f_draine(energy_eV) result(f)
    use prizmo_commons, only: hplanck_eV
    implicit none
    real*8,intent(in)::energy_eV
    real*8::f

    ! erg/cm2
    f = 0d0

    ! Draine flux non-zero only within 6, 13.6 eV
    if((energy_eV > 13.6d0) .or. (energy_eV < 6d0)) then
      return
    end if

    ! return Draine flux
    f = (1.658d6*energy_ev - 2.152d5*energy_ev**2 + 6.919d3*energy_ev**3) &
      * energy_ev * hplanck_eV

  end function f_draine

  ! ***************************
  ! eV/cm2
  function f_black_body(energy_ev, Tbb) result(f)
    use prizmo_commons, only: kboltzmann_ev, stefan_boltzmann_ev, &
    clight, hplanck_eV
    implicit none
    real*8,intent(in)::energy_ev, Tbb
    real*8::f, xexp

    ! exponential argument
    xexp = energy_ev / kboltzmann_eV / Tbb

    f = 0d0

    ! large and small exponents are excluded to avoid under/overflow
    if((xexp > 3d2) .or. (xexp < 1d-10)) then
      return
    end if

    f = 2d0 * energy_ev**3 / hplanck_eV**2 / clight**2 &
    / (exp(xexp) - 1d0)

    ! normalize
    !f = f * integral_draine() / (stefan_boltzmann_eV * Tbb**4)

  end function f_black_body

  ! *************************
  ! xray flux in the range [6, 1183.5 eV], flux in eV / cm2
  function f_xray(energy_eV) result(f)
    use prizmo_commons
    implicit none
    real*8,intent(in)::energy_eV
    real*8::f, loge

    ! eV / cm2
    f = 0d0

    ! range is because of fit validity
    if((energy_eV > 1.1835d4) .or. (energy_eV < 6d0)) then
      return
    end if

    ! fit is in log space
    loge = log10(energy_eV)

    ! return Draine flux
    f = 16.54963901 - 29.05076878 * loge + 20.25163934 * loge**2 &
        - 6.74466784 * loge**3 + 1.09544828 * loge**4 - 0.07121136 * loge**5

    ! 0.1663 factor is to have \xi = 1
    f = 1d1**f / energy_eV / 4d0 / pi / 0.1663 / ev2erg * hplanck_eV

  end function f_xray

end module prizmo_photo
