module prizmo_attenuation
  use prizmo_commons
  real*8::dust_kappa(nphoto)
contains

  ! ****************************
  ! attenuate the flux in a cell with a given length and chemistry
  subroutine attenuate(flux, n, Tgas, length)
    use prizmo_commons
    use prizmo_utils
    use prizmo_photo
    use prizmo_tdust
    implicit none
    real*8,intent(in)::length, n(nmols), Tgas
    real*8,intent(inout)::flux(nphoto)
    real*8::tau(nphoto), Tdust, jem(nphoto)
    integer::i

    tau(:) = get_tau(n(:), Tgas, length)

    !Tdust = get_Tdust(n, Tgas, flux(:))

    !jem(:) = get_black_body_flux(Tdust) * 2d0 * pi

    flux(:) = flux(:) * exp(-tau(:))
    !flux(:) = flux(:) - flux(:) * tau(:) + jem(:)

    !do i=1, nphoto
    !    flux(i) = max(flux(i), 0d0)
    !end do

  end subroutine attenuate

  ! ****************************
  ! get optical depth, dimensionless
  function get_tau(n, Tgas, length) result(tau)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::length, n(nmols), Tgas
    real*8::tau(nphoto), inv_ntot

    tau(:) = 0d0  ! 1/cm

    inv_ntot = 1d0 / get_ntot(n(:))

    !!BEGIN_ATTENUATION
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! H2+ -> H+ + H
    if(n(idx_H2j) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 260) * n(idx_H2j)
    end if

    ! H3+ -> H2+ + H
    if(n(idx_H3j) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 261) * n(idx_H3j)
    end if

    ! H3+ -> H2 + H+
    if(n(idx_H3j) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 262) * n(idx_H3j)
    end if

    ! C -> C+ + E
    if(n(idx_C) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 263) * n(idx_C)
    end if

    ! CH -> CH+ + E
    if(n(idx_CH) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 264) * n(idx_CH)
    end if

    ! CH -> C + H
    if(n(idx_CH) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 265) * n(idx_CH)
    end if

    ! CH+ -> C+ + H
    if(n(idx_CHj) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 266) * n(idx_CHj)
    end if

    ! CH2 -> CH2+ + E
    ! xsecs < 1e-40 cm2, skip
    ! CH2 -> CH + H
    if(n(idx_CH2) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 268) * n(idx_CH2)
    end if

    ! CH2+ -> CH+ + H
    if(n(idx_CH2j) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 269) * n(idx_CH2j)
    end if

    ! CH3 -> CH3+ + E
    if(n(idx_CH3) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 270) * n(idx_CH3)
    end if

    ! CH3 -> CH2 + H
    if(n(idx_CH3) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 271) * n(idx_CH3)
    end if

    ! CH3 -> CH + H2
    if(n(idx_CH3) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 272) * n(idx_CH3)
    end if

    ! CH4 -> CH3 + H
    if(n(idx_CH4) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 273) * n(idx_CH4)
    end if

    ! CH4 -> CH2 + H2
    if(n(idx_CH4) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 274) * n(idx_CH4)
    end if

    ! CH4 -> CH + H2 + H
    if(n(idx_CH4) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 275) * n(idx_CH4)
    end if

    ! OH -> OH+ + E
    ! xsecs < 1e-40 cm2, skip
    ! OH -> O + H
    if(n(idx_OH) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 277) * n(idx_OH)
    end if

    ! OH+ -> O + H+
    if(n(idx_OHj) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 278) * n(idx_OHj)
    end if

    ! H2O -> OH + H
    if(n(idx_H2O) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 279) * n(idx_H2O)
    end if

    ! H2O -> H2O+ + E
    if(n(idx_H2O) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 280) * n(idx_H2O)
    end if

    ! O2 -> O + O
    if(n(idx_O2) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 281) * n(idx_O2)
    end if

    ! O2 -> O2+ + E
    if(n(idx_O2) * inv_ntot > 1d-14) then
        tau(:) = tau(:) + xsecs(:, 282) * n(idx_O2)
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_ATTENUATION

    ! add dust opacity (note: 1/cm)
    tau(:) = tau(:) + d2g * get_rho(n(:)) * dust_kappa(:)

    ! multiply by the segment length to have dimensionless tau
    tau(:) = tau(:) * length

  end	function get_tau

  ! ***************************
  ! load dust opacity from file to dust_kappa common variable, cm2/g
  subroutine load_dust_opacity(fname)
    use prizmo_commons
    implicit none
    character(len=*),intent(in)::fname
    integer::unit, i
    real*8::energy

    print *, "loading opacity file "//trim(fname)//"..."
    open(newunit=unit,file=fname, status='old')
    do i=1,nphoto
       read(unit, *) energy, dust_kappa(i)
    end do
    close(unit)

  end subroutine load_dust_opacity

end module prizmo_attenuation
