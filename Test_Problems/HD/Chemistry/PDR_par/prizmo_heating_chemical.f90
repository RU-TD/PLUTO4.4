module prizmo_heating_chemical
contains

  ! ***********************
  ! get chemical cooling (and heating), erg/cm3/s
  function heating_chemical(n, Tgas) result(heat)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(nmols), Tgas
    real*8::heat

    heat = 0d0

    !!BEGIN_HEATING_CHEMICAL
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! H + H -> H2 (H2 formation on dust heating, 4.5 eV)
    heat = heat + (7.209794e-12) * kall(289)*n(idx_H)*n(idx_H)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_HEATING_CHEMICAL

  end function heating_chemical

end module prizmo_heating_chemical
