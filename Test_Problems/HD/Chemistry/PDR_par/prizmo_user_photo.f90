module prizmo_user_photo
contains

  ! *************************
  function prizmo_get_energy() result(energy)
    use prizmo_commons, only: nphoto, energy_grid
    real*8::energy(nphoto)

    energy(:) = energy_grid(:)

  end function prizmo_get_energy

  ! *************************
  ! load flux from file, with format
  ! energy/[eV], flux/[eV/cm2/s/Hz/sr]
  function prizmo_load_flux(fname) result(jflux)
    use prizmo_commons
    implicit none
    character(len=*),intent(in)::fname
    real*8::jflux(nphoto), energy
    integer::i, unit

    open(newunit=unit, file=trim(fname), status="old")
    do i=1,nphoto
      read(unit, *) energy, jflux(i)
      if(abs(energy_grid(i) - energy) / energy_grid(i) > 1d-4) then
        print *, "ERROR: when loading flux from "//trim(fname)//" energy grid is not matching!"
        print *, "found, eV", energy
        print *, "expected, eV", energy_grid(i)
        print *, "bin number", i
        stop
      end if
    end do
    close(unit)

  end function prizmo_load_flux

  ! *************************
  function prizmo_get_draine_flux() result(jflux)
    use prizmo_commons, only: nphoto, energy_grid
    implicit none
    real*8::jflux(nphoto)
    integer::i

    ! loop on photo bins
    do i=1,nphoto
      jflux(i) = prizmo_f_draine(energy_grid(i))
    end do

  end function prizmo_get_draine_flux

  ! *************************
  function prizmo_get_black_body_flux(Tbb) result(jflux)
    use prizmo_commons, only: nphoto, energy_grid
    implicit none
    real*8,intent(in)::Tbb
    real*8::jflux(nphoto)
    integer::i

    ! loop on photo bins
    do i=1,nphoto
      jflux(i) = prizmo_f_black_body(energy_grid(i), Tbb)
    end do

  end function prizmo_get_black_body_flux

  ! *************************
  function prizmo_get_xray_flux() result(jflux)
    use prizmo_commons, only: nphoto, energy_grid
    implicit none
    real*8::jflux(nphoto)
    integer::i

    ! loop on photo bins
    do i=1,nphoto
      jflux(i) = prizmo_f_xray(energy_grid(i))
    end do

  end function prizmo_get_xray_flux

  ! *************************
  function prizmo_f_draine(energy_eV) result(f)
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

  end function prizmo_f_draine

  ! ***************************
  ! eV/cm2
  function prizmo_f_black_body(energy_ev, Tbb) result(f)
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

  end function prizmo_f_black_body

  ! *************************
  ! xray flux in the range [6, 1183.5 eV], flux in eV / cm2
  function prizmo_f_xray(energy_eV) result(f)
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

  end function prizmo_f_xray

  ! ***************************
  ! get ratio between current J(E)/E integral in the range [6,13.6] eV
  ! and the corresponding Draine flux integral in the same range
  function prizmo_get_Gnot(jflux) result(Gnot)
    use prizmo_commons
    use prizmo_rates_photo
    implicit none
    real*8,intent(in)::jflux(nphoto)
    real*8::Gnot

    Gnot = get_Gnot(jflux(:))

  end function prizmo_get_Gnot

  ! **********************
  subroutine prizmo_attenuate(flux, n, Tgas, lenght)
    use prizmo_commons
    use prizmo_attenuation
    implicit none
    real*8,intent(in)::lenght, n(nmols), Tgas
    real*8,intent(inout)::flux(nphoto)

    call attenuate(flux(:), n(:), Tgas, lenght)

  end subroutine prizmo_attenuate

  ! ***********************
  function prizmo_get_tau(n, Tgas, lenght) result(tau)
    use prizmo_commons
    use prizmo_attenuation
    implicit none
    real*8,intent(in)::lenght, n(nmols), Tgas
    real*8::tau(nphoto)

    tau(:) = get_tau(n(:), Tgas, lenght)

  end function prizmo_get_tau

end module prizmo_user_photo
