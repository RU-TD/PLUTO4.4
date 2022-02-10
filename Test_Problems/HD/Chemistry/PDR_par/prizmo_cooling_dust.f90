module prizmo_cooling_dust
contains

  ! ***********************
  function cooling_dust(n, Tgas, Tdust, jflux) result(cool)
    use prizmo_commons
    use prizmo_utils
    use prizmo_tdust
    use prizmo_dust_evaporation
    implicit none
    real*8,intent(in)::n(nmols), Tgas, Tdust, jflux(nphoto)
    real*8::cool, pre

    ! dust_mode = 3 means Tdust=Tgas, hence no cooling
    if(tdust_mode == 3) then
      cool = 0d0
    else
      pre = get_mu(n(:)) * get_ntot(n(:)) * proton_mass * d2g &
          / pi43 / dust_rho_bulk / dust_fact4

      !cool = (Tgas - Tdust) * 2d0 * sqrt(8d0 * kboltzmann * Tgas / pi / proton_mass) &
      !     * kboltzmann * pi * pre * dust_fact3 * 0.5
      cool = (femission(Tdust) - get_absorption_term(jflux(:))) * ev2erg * pre !&
      !         * get_dust_fevap(Tdust)
    end if

  end function cooling_dust

end module prizmo_cooling_dust
