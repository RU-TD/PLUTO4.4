module prizmo_heating_CR
  use prizmo_commons
contains

  ! ***************
  function heating_CR(x) result(heat)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::heat

    heat = user_cr * (5.5d-12 * x(idx_H) + 2.5d-11 * x(idx_H2))

  end function heating_CR

end module prizmo_heating_CR
