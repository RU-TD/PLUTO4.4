module prizmo_cooling_chemical
contains

  ! ***********************
  ! get chemical cooling (and heating), erg/cm3/s
  function cooling_chemical(n, Tgas) result(cool)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(nmols), Tgas
    real*8::cool
    integer::i

    cool = 0d0

    !!BEGIN_COOLING_CHEMICAL
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! H+ + E -> H (cooling)
    cool = cool + kboltzmann * Tgas * kall(71)*n(idx_Hj)*n(idx_E)

    ! He+ + E -> He (cooling)
    cool = cool + kboltzmann * Tgas * kall(72)*n(idx_Hej)*n(idx_E)

    ! C+ + E -> C (cooling)
    cool = cool + kboltzmann * Tgas * kall(73)*n(idx_Cj)*n(idx_E)

    ! CH3+ + E -> CH3 (cooling)
    cool = cool + kboltzmann * Tgas * kall(74)*n(idx_CH3j)*n(idx_E)

    ! O+ + E -> O (cooling)
    cool = cool + kboltzmann * Tgas * kall(75)*n(idx_Oj)*n(idx_E)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_COOLING_CHEMICAL

  end function cooling_chemical

end module prizmo_cooling_chemical
