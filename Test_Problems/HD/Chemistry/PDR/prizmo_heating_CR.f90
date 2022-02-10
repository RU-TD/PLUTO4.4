module prizmo_heating_CR

contains
  ! ***********************
  function heating_CR(n, Tgas) result(heat)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(nmols), Tgas
    real*8::heat

    heat = variable_crflux * (5.5d-12 * max(n(idx_H), 1d-40) &
         + 2.5d-11 * max(n(idx_H2), 1d-40))

  end function heating_CR
end module prizmo_heating_CR
