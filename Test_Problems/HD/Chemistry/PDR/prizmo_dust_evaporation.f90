module prizmo_dust_evaporation
contains

! **************************
! get the fraction of active dust, i.e. that is not evaporated
function get_dust_fevap(Tdust) result(f)
  implicit none
  real*8,intent(in)::Tdust
  real*8::f

  f = 1d0 - 1d0 / (1d0 + exp((-Tdust + 2d3) * 5d-3))

end function get_dust_fevap

end module prizmo_dust_evaporation
