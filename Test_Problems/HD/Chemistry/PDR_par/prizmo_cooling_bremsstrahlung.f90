module prizmo_cooling_bremsstrahlung
contains

  ! ***********************
  ! get bremsstrahlung cooling, erg/cm3/s
  function cooling_bremsstrahlung(n, Tgas) result(cool)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(nmols), Tgas
    real*8::cool, ions_sum_factor
    real*8,parameter::gaunt_factor = 1.5d0 !mean value

    !!BEGIN_COOLING_BREMSSTRAHLUNG
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ions_sum_factor = 0e0 &
        + 1.000000 * n(idx_CH2j) &
        + 1.000000 * n(idx_CHj) &
        + 1.000000 * n(idx_CH3j) &
        + 1.000000 * n(idx_Hej) &
        + 1.000000 * n(idx_Hj) &
        + 1.000000 * n(idx_Cj) &
        + 1.000000 * n(idx_Oj) &
        + 1.000000 * n(idx_H2j) &
        + 1.000000 * n(idx_COj) &
        + 1.000000 * n(idx_H3j) &
        + 1.000000 * n(idx_CH4j) &
        + 1.000000 * n(idx_OHj) &
        + 1.000000 * n(idx_CH5j) &
        + 1.000000 * n(idx_H2Oj) &
        + 1.000000 * n(idx_H3Oj) &
        + 1.000000 * n(idx_HCOj) &
        + 1.000000 * n(idx_O2j)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_COOLING_BREMSSTRAHLUNG

    cool = 1.42d-27 * gaunt_factor * sqrt(Tgas) * ions_sum_factor * n(idx_e)

    cool = max(cool, 0d0)

  end function cooling_bremsstrahlung

end module prizmo_cooling_bremsstrahlung
