module prizmo_heating_PAH
contains

  ! **************************
  function heating_PAH(n, Tgas, G0, Av) result(heat)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas, G0, Av
    real*8::heat, dd, xpah, epah, ne
    real*8,parameter::fpah=0.2

    ! check negative electrons
    ne = max(n(idx_E), 1d-40)

    dd = get_Hnuclei(n(:))
    xpah = G0 * sqrt(Tgas) / ne
    epah = 0.0487 / (1d0 + 4d-3 * xpah**0.73)

    heat = 1d-24 * fpah * dd * epah * G0 * exp(-2.5 * Av) * dust_d2g / 1d-2

  end function heating_PAH

end module prizmo_heating_PAH
