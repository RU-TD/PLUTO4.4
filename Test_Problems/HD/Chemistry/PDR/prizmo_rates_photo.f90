module prizmo_rates_photo
contains

  ! **************************
  subroutine compute_rates_photo(n, Tgas, jflux)
    use prizmo_commons
    use prizmo_self_shielding
    implicit none
    real*8,intent(in)::n(nmols)
    real*8,intent(in)::Tgas
    real*8,intent(in)::jflux(nphoto)

    !!BEGIN_RATES_PHOTO
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    !H2 -> H + H
    kall(258) = (integrate_rate_photo(jflux(:), 258)) * get_self_shielding_H2(variable_NH2_incoming, Tgas)

    !CO -> C + O
    kall(259) = (integrate_rate_photo(jflux(:), 259)) * get_self_shielding_CO(variable_NH2_incoming, variable_NCO_incoming)

    !H2+ -> H+ + H
    kall(260) = integrate_rate_photo(jflux(:), 260)

    !H3+ -> H2+ + H
    kall(261) = (integrate_rate_photo(jflux(:), 261)) * 0.5

    !H3+ -> H2 + H+
    kall(262) = (integrate_rate_photo(jflux(:), 262)) * 0.5

    !C -> C+ + E
    kall(263) = integrate_rate_photo(jflux(:), 263)

    !CH -> CH+ + E
    kall(264) = integrate_rate_photo(jflux(:), 264)

    !CH -> C + H
    kall(265) = integrate_rate_photo(jflux(:), 265)

    !CH+ -> C+ + H
    kall(266) = integrate_rate_photo(jflux(:), 266)

    !CH2 -> CH2+ + E
    ! xsecs always < 1.000000e-40 cm2, skip

    !CH2 -> CH + H
    kall(268) = integrate_rate_photo(jflux(:), 268)

    !CH2+ -> CH+ + H
    kall(269) = integrate_rate_photo(jflux(:), 269)

    !CH3 -> CH3+ + E
    kall(270) = integrate_rate_photo(jflux(:), 270)

    !CH3 -> CH2 + H
    kall(271) = integrate_rate_photo(jflux(:), 271)

    !CH3 -> CH + H2
    kall(272) = integrate_rate_photo(jflux(:), 272)

    !CH4 -> CH3 + H
    kall(273) = integrate_rate_photo(jflux(:), 273)

    !CH4 -> CH2 + H2
    kall(274) = integrate_rate_photo(jflux(:), 274)

    !CH4 -> CH + H2 + H
    kall(275) = integrate_rate_photo(jflux(:), 275)

    !OH -> OH+ + E
    ! xsecs always < 1.000000e-40 cm2, skip

    !OH -> O + H
    kall(277) = integrate_rate_photo(jflux(:), 277)

    !OH+ -> O + H+
    kall(278) = integrate_rate_photo(jflux(:), 278)

    !H2O -> OH + H
    kall(279) = integrate_rate_photo(jflux(:), 279)

    !H2O -> H2O+ + E
    kall(280) = integrate_rate_photo(jflux(:), 280)

    !O2 -> O + O
    kall(281) = integrate_rate_photo(jflux(:), 281)

    !O2 -> O2+ + E
    kall(282) = integrate_rate_photo(jflux(:), 282)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_RATES_PHOTO

  end subroutine compute_rates_photo

  ! *************************
  ! compute rate reaction, 1/s
  function integrate_rate_photo(jflux, idx) result(k)
    use prizmo_commons
    implicit none
    real*8,intent(in)::jflux(nphoto)
    integer,intent(in)::idx
    integer::i
    real*8::k, f0, f1, x0, x1

    k = sum(Jflux(:) * xsecs_trapz(:, idx)) / 2d0
    k = k / hplanck_eV

    return

    !k = 0d0  ! eV
    !f0 = xsecs(1, idx) * jflux(1) * inv_energy_grid(1)
    !x0 = energy_grid(1)  ! eV
    !do i=2,nphoto
    !   f1 = xsecs(i, idx) * jflux(i) * inv_energy_grid(i)
    !   x1 = energy_grid(i)
    !   k = k + 0.5 * (f0 + f1) * (x1 - x0)
    !   f0 = f1
    !   x0 = x1
    !end do

    ! convert to 1/s
    !k = k / hplanck_eV

  end function integrate_rate_photo


  ! **************************
  ! G0 as the ratio betwenn the current integral of jflux
  ! over the Draine energy interval scaled by the integral of the
  ! corresponding Draine flux (6, 13.6eV)
  ! NOTE: this is for rates, hence J(E) / E is the integration kernel
  function get_Gnot_rates(jflux) result(Gnot)
    use prizmo_commons
    implicit none
    real*8,intent(in)::jflux(nphoto)
    real*8::Gnot, f1, f2, intf
    real*8,parameter::intf0=1.547d7*hplanck_eV ! J(E)/E intgeral, eV/cm2
    integer::i

    intf = 0d0  ! eV/cm2
    do i=imin_fDraine,imax_fDraine
      f1 = jflux(i) / energy_grid(i)
      f2 = jflux(i+1) / energy_grid(i+1)
      intf = intf + (f1 + f2) / 2d0 * (energy_grid(i+1) - energy_grid(i))
    end do

    Gnot = intf / intf0

  end function get_Gnot_rates

  ! **************************
  ! G0 as the ratio betwenn the current integral of jflux
  ! over the Draine energy interval scaled by the integral of the
  ! corresponding Draine flux (6, 13.6eV)
  ! NOTE: this is classic definition, hence J(E) is the integration kernel
  function get_Gnot(jflux) result(Gnot)
    use prizmo_commons
    implicit none
    real*8,intent(in)::jflux(nphoto)
    real*8::Gnot, f1, f2, intf
    real*8,parameter::intf0=1.33113d8*hplanck_eV ! J(E) intgeral, eV/cm2
    integer::i

    intf = 0d0  ! eV/cm2
    do i=imin_fDraine,imax_fDraine
      f1 = jflux(i)
      f2 = jflux(i+1)
      intf = intf + (f1 + f2) / 2d0 * (energy_grid(i+1) - energy_grid(i))
    end do

    Gnot = intf / intf0

  end function get_Gnot

  ! **************************
  ! H2 rate as the ratio betwenn the current integral of jflux
  ! over the Draine energy interval scaled by the integral of the
  ! corresponding Draine flux
  function get_H2_photodissociation_thin(jflux) result(kH2)
    use prizmo_commons
    implicit none
    real*8,intent(in)::jflux(nphoto)
    real*8::kH2

    ! H2 photodissociation rate reference is from Leiden database, ISRF
    ! Gnot factor is normalized with a J(E) / E kernel
    ! https://home.strw.leidenuniv.nl/~ewine/photo/display_h2_ca2bf3f6b7e18a508253e9521510a4b5.html
    kH2 = 5.68d-11 * get_Gnot_rates(jflux(:))  ! 1/s

  end function get_H2_photodissociation_thin

end module prizmo_rates_photo
