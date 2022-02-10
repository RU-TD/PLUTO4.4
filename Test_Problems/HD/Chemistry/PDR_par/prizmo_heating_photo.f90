module prizmo_heating_photo
  use prizmo_commons
  ! heating threshold energy, eV
  real*8::xsecs_heating_threshold(nrea)
  real*8::pre_integral_heat(nrea)
contains

  ! **************************
  ! photo heating based on Av, including H2 dissociation, erg/s/cm3
  function heating_photo(n, Tgas, jflux) result(heat)
    use prizmo_commons
    use prizmo_self_shielding
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas, jflux(nphoto)
    real*8::heat, dd, ncrn, ncrd1, ncrd2, yH, yH2
    real*8::ncr, hf, Rdiss, heat_eV, nH2, f_ntot, xe

    nH2 = max(n(idx_H2), 0d0)

    dd = get_Hnuclei(n(:))
    ncrn = 1d6 / sqrt(Tgas)
    ncrd1 = 1.6 * exp(-(4d2 / Tgas)**2) + 1d-40
    ncrd2 = 1.4 * exp(- 1.2d4 / (Tgas + 1.2d3)) + 1d-40

    yH = n(idx_H) / dd + 1d-40
    yH2 = nH2 / dd + 1d-40

    ncr = ncrn / (ncrd1 * yH + ncrd2 * yH2)
    hf = dd / (dd + ncr)

    ! default if H2 photodissociation is not present
    Rdiss = 0d0

    !!BEGIN_H2_PHOTODISS_RATE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    Rdiss = kall(258)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_H2_PHOTODISS_RATE

    ! H2 dissociation heating, erg/s/cm3
    heat = (6.4d-13 + 2.7d-11 * hf) * Rdiss * nH2

    ! ntot threshold, cm-3
    f_ntot = 1d-10 * get_ntot(n(:))

    ! heating in eV for multiline frequency
    heat_eV = 0d0

    !!BEGIN_ENERGY_EXCESS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! CO -> C + O
    heat_eV = heat_eV + pre_integral_heat(259) * max(n(idx_CO), 0d0)
    ! H2+ -> H+ + H
    heat_eV = heat_eV + pre_integral_heat(260) * max(n(idx_H2j), 0d0)
    ! H3+ -> H2+ + H
    heat_eV = heat_eV + pre_integral_heat(261) * max(n(idx_H3j), 0d0)
    ! H3+ -> H2 + H+
    heat_eV = heat_eV + pre_integral_heat(262) * max(n(idx_H3j), 0d0)
    ! C -> C+ + E
    heat_eV = heat_eV + pre_integral_heat(263) * max(n(idx_C), 0d0)
    ! CH -> CH+ + E
    heat_eV = heat_eV + pre_integral_heat(264) * max(n(idx_CH), 0d0)
    ! CH -> C + H
    heat_eV = heat_eV + pre_integral_heat(265) * max(n(idx_CH), 0d0)
    ! CH+ -> C+ + H
    heat_eV = heat_eV + pre_integral_heat(266) * max(n(idx_CHj), 0d0)
    ! CH2 -> CH2+ + E
    heat_eV = heat_eV + pre_integral_heat(267) * max(n(idx_CH2), 0d0)
    ! CH2 -> CH + H
    heat_eV = heat_eV + pre_integral_heat(268) * max(n(idx_CH2), 0d0)
    ! CH2+ -> CH+ + H
    heat_eV = heat_eV + pre_integral_heat(269) * max(n(idx_CH2j), 0d0)
    ! CH3 -> CH3+ + E
    heat_eV = heat_eV + pre_integral_heat(270) * max(n(idx_CH3), 0d0)
    ! CH3 -> CH2 + H
    heat_eV = heat_eV + pre_integral_heat(271) * max(n(idx_CH3), 0d0)
    ! CH3 -> CH + H2
    heat_eV = heat_eV + pre_integral_heat(272) * max(n(idx_CH3), 0d0)
    ! CH4 -> CH3 + H
    heat_eV = heat_eV + pre_integral_heat(273) * max(n(idx_CH4), 0d0)
    ! CH4 -> CH2 + H2
    heat_eV = heat_eV + pre_integral_heat(274) * max(n(idx_CH4), 0d0)
    ! CH4 -> CH + H2 + H
    heat_eV = heat_eV + pre_integral_heat(275) * max(n(idx_CH4), 0d0)
    ! OH -> OH+ + E
    heat_eV = heat_eV + pre_integral_heat(276) * max(n(idx_OH), 0d0)
    ! OH -> O + H
    heat_eV = heat_eV + pre_integral_heat(277) * max(n(idx_OH), 0d0)
    ! OH+ -> O + H+
    heat_eV = heat_eV + pre_integral_heat(278) * max(n(idx_OHj), 0d0)
    ! H2O -> OH + H
    heat_eV = heat_eV + pre_integral_heat(279) * max(n(idx_H2O), 0d0)
    ! H2O -> H2O+ + E
    heat_eV = heat_eV + pre_integral_heat(280) * max(n(idx_H2O), 0d0)
    ! O2 -> O + O
    heat_eV = heat_eV + pre_integral_heat(281) * max(n(idx_O2), 0d0)
    ! O2 -> O2+ + E
    heat_eV = heat_eV + pre_integral_heat(282) * max(n(idx_O2), 0d0)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_ENERGY_EXCESS

    ! effciency following Xu & McCray 1991
    xe = get_electrons(n(:)) / get_ntot(n(:))
    !if(xe < 0.95) then
    !  heat_eV = heat_eV * 0.9971 * (1d0 - (1d0 - xe**0.2663)**1.3163)
    !end if

    ! effciency with a continuos function, but with a 10% error on Xu & McCray 1991
    !FIXME heat_eV = heat_eV * 1d1**(log10(xe + 1d-40) * 0.25452)

    ! add photoheating to H2 photodissociation heating
    heat = heat + heat_eV * ev2erg

  end function heating_photo

  ! ****************************
  ! pre integrate heating factors on the current radiation flux
  subroutine pre_integrate_photoheating(jflux)
    use prizmo_commons
    use prizmo_self_shielding
    implicit none
    real*8,intent(in)::jflux(nphoto)
    integer::i

    ! set default value
    pre_integral_heat(:) = 0d0

    !!BEGIN_PHOTOHEATING_PREINTEGRATE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! CO -> C + O
    pre_integral_heat(259) = (integrate_rate_photo_heating(jflux(:), 259)) * get_self_shielding_CO(variable_NH2_incoming, variable_NCO_incoming)

    ! H2+ -> H+ + H
    pre_integral_heat(260) = integrate_rate_photo_heating(jflux(:), 260)

    ! H3+ -> H2+ + H
    pre_integral_heat(261) = (integrate_rate_photo_heating(jflux(:), 261)) * 0.5

    ! H3+ -> H2 + H+
    pre_integral_heat(262) = (integrate_rate_photo_heating(jflux(:), 262)) * 0.5

    ! C -> C+ + E
    pre_integral_heat(263) = integrate_rate_photo_heating(jflux(:), 263)

    ! CH -> CH+ + E
    pre_integral_heat(264) = integrate_rate_photo_heating(jflux(:), 264)

    ! CH -> C + H
    pre_integral_heat(265) = integrate_rate_photo_heating(jflux(:), 265)

    ! CH+ -> C+ + H
    pre_integral_heat(266) = integrate_rate_photo_heating(jflux(:), 266)

    ! CH2 -> CH2+ + E
    pre_integral_heat(267) = integrate_rate_photo_heating(jflux(:), 267)

    ! CH2 -> CH + H
    pre_integral_heat(268) = integrate_rate_photo_heating(jflux(:), 268)

    ! CH2+ -> CH+ + H
    pre_integral_heat(269) = integrate_rate_photo_heating(jflux(:), 269)

    ! CH3 -> CH3+ + E
    pre_integral_heat(270) = integrate_rate_photo_heating(jflux(:), 270)

    ! CH3 -> CH2 + H
    pre_integral_heat(271) = integrate_rate_photo_heating(jflux(:), 271)

    ! CH3 -> CH + H2
    pre_integral_heat(272) = integrate_rate_photo_heating(jflux(:), 272)

    ! CH4 -> CH3 + H
    pre_integral_heat(273) = integrate_rate_photo_heating(jflux(:), 273)

    ! CH4 -> CH2 + H2
    pre_integral_heat(274) = integrate_rate_photo_heating(jflux(:), 274)

    ! CH4 -> CH + H2 + H
    pre_integral_heat(275) = integrate_rate_photo_heating(jflux(:), 275)

    ! OH -> OH+ + E
    pre_integral_heat(276) = integrate_rate_photo_heating(jflux(:), 276)

    ! OH -> O + H
    pre_integral_heat(277) = integrate_rate_photo_heating(jflux(:), 277)

    ! OH+ -> O + H+
    pre_integral_heat(278) = integrate_rate_photo_heating(jflux(:), 278)

    ! H2O -> OH + H
    pre_integral_heat(279) = integrate_rate_photo_heating(jflux(:), 279)

    ! H2O -> H2O+ + E
    pre_integral_heat(280) = integrate_rate_photo_heating(jflux(:), 280)

    ! O2 -> O + O
    pre_integral_heat(281) = integrate_rate_photo_heating(jflux(:), 281)

    ! O2 -> O2+ + E
    pre_integral_heat(282) = integrate_rate_photo_heating(jflux(:), 282)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_PHOTOHEATING_PREINTEGRATE

    ! DEBUG: uncomment here below if needed
    !do i=1,nrea
    !  if(pre_integral_heat(i) > 0d0) print *, i, pre_integral_heat(i)
    !end do

  end subroutine pre_integrate_photoheating

  ! *************************
  ! compute rate heating, erg/s
  function integrate_rate_photo_heating(jflux, idx) result(k)
    use prizmo_commons
    implicit none
    real*8,intent(in)::jflux(nphoto)
    integer,intent(in)::idx
    real*8::k, kh

    ! xsecs integral arrays are loaded in prizmo_xsecs.f90
    k = sum(Jflux(:) * xsecs_trapz(:, idx)) / 2d0
    kh = sum(Jflux(:) * xsecs_trapz_heat(:, idx)) / 2d0

    ! xsecs_heating_threshold (eV) is loaded in this module (check next functions)
    k = (kh - xsecs_heating_threshold(idx) * k) / hplanck_eV

    if(k < 0d0) then
       print *, "WARNING: negative photo heating integral, idx, rate", idx, k
       k = max(k, 0d0)
    end if

  end function integrate_rate_photo_heating

  ! ************************
  ! load thresholds (eV) from file for heating
  subroutine load_heating_thresholds(fname)
    use prizmo_commons
    implicit none
    integer::i, unit
    character(len=*),intent(in)::fname

    open(newunit=unit, file=trim(fname), status="old")
    do i=1,nrea
       read(unit, *) xsecs_heating_threshold(i)
    end do
    close(unit)

  end subroutine load_heating_thresholds

end module prizmo_heating_photo
