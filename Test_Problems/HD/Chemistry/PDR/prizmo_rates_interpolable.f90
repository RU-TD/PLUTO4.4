module prizmo_rates_interpolable
contains

  !*******************
  ! compute rates and store into commons kall(:)
  subroutine evaluate_interpolable(n, Tgas, Tdust)
    use prizmo_commons
    use prizmo_self_shielding
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas, Tdust
    real*8::inv_Tgas, sqrt_Tgas, Te, lnTe, invTe, t
    integer::i

    ! frequently used quantities
    inv_Tgas = 1d0 / Tgas
    sqrt_Tgas = sqrt(Tgas)
    Te = Tgas * 8.617343d-5  ! Tgas in eV (eV)
    lnTe = log(Te)  ! ln of Te (#)
    invTe = 1d0 / Te  ! inverse of T (1/eV)

    !!BEGIN_RATES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! H + CH -> C + H2
    kall(1) = 2.70e-11 * (Tgas / 3d2)**(0.38)

    ! H + CH3 -> CH2 + H2
    kall(2) = 1.00e-10 * exp(-7600.0 / Tgas)

    ! H + CH4 -> CH3 + H2
    kall(3) = 5.94e-13 * (Tgas / 3d2)**(3.00) * exp(-4045.0 / Tgas)

    ! H + OH -> O + H2
    kall(4) = 6.99e-14 * (Tgas / 3d2)**(2.80) * exp(-1950.0 / Tgas)

    ! H + H2O -> OH + H2
    kall(5) = 1.59e-11 * (Tgas / 3d2)**(1.20) * exp(-9610.0 / Tgas)

    ! H + CO -> OH + C
    kall(6) = 1.10e-10 * (Tgas / 3d2)**(0.50) * exp(-77700.0 / Tgas)

    ! H + O2 -> OH + O
    kall(7) = 2.61e-10 * exp(-8156.0 / Tgas)

    ! H2 + C -> CH + H
    kall(8) = 6.64e-10 * exp(-11700.0 / Tgas)

    ! H2 + CH -> CH2 + H
    kall(9) = 5.46e-10 * exp(-1943.0 / Tgas)

    ! H2 + CH2 -> CH3 + H
    kall(10) = 5.18e-11 * (Tgas / 3d2)**(0.17) * exp(-6400.0 / Tgas)

    ! H2 + CH3 -> CH4 + H
    kall(11) = 6.86e-14 * (Tgas / 3d2)**(2.74) * exp(-4740.0 / Tgas)

    ! H2 + O -> OH + H
    kall(12) = 3.14e-13 * (Tgas / 3d2)**(2.70) * exp(-3150.0 / Tgas)

    ! H2 + OH -> H2O + H
    kall(13) = 2.05e-12 * (Tgas / 3d2)**(1.52) * exp(-1736.0 / Tgas)

    ! H2 + O2 -> OH + OH
    kall(14) = 3.16e-10 * exp(-21890.0 / Tgas)

    ! C + CH2 -> CH + CH
    kall(15) = 2.69e-12 * exp(-23550.0 / Tgas)

    ! C + OH -> O + CH
    kall(16) = 2.25e-11 * (Tgas / 3d2)**(0.50) * exp(-14800.0 / Tgas)

    ! C + OH -> CO + H
    kall(17) = 1.10e-10 * (Tgas / 3d2)**(0.50)

    ! CH + O -> OH + C
    kall(18) = 2.52e-11 * exp(-2381.0 / Tgas)

    ! CH + CH4 -> CH3 + CH2
    kall(19) = 2.28e-11 * (Tgas / 3d2)**(0.70) * exp(-3000.0 / Tgas)

    ! CH2 + CH2 -> CH3 + CH
    kall(20) = 4.00e-10 * exp(-5000.0 / Tgas)

    ! CH2 + O -> OH + CH
    kall(21) = 4.98e-10 * exp(-6000.0 / Tgas)

    ! CH2 + CH4 -> CH3 + CH3
    kall(22) = 7.13e-12 * exp(-5050.0 / Tgas)

    ! CH2 + OH -> O + CH3
    kall(23) = 1.44e-11 * (Tgas / 3d2)**(0.50) * exp(-3000.0 / Tgas)

    ! CH2 + OH -> H2O + CH
    kall(24) = 1.44e-11 * (Tgas / 3d2)**(0.50) * exp(-3000.0 / Tgas)

    ! CH2 + O2 -> CO + H2O
    kall(25) = 2.48e-10 * (Tgas / 3d2)**(-3.30) * exp(-1443.0 / Tgas)

    ! CH3 + CH3 -> CH4 + CH2
    kall(26) = 7.13e-12 * exp(-5052.0 / Tgas)

    ! CH3 + OH -> CH4 + O
    kall(27) = 3.27e-14 * (Tgas / 3d2)**(2.20) * exp(-2240.0 / Tgas)

    ! CH3 + OH -> H2O + CH2
    kall(28) = 1.20e-10 * exp(-1400.0 / Tgas)

    ! CH3 + H2O -> OH + CH4
    kall(29) = 2.30e-15 * (Tgas / 3d2)**(3.47) * exp(-6681.0 / Tgas)

    ! O + CH4 -> OH + CH3
    kall(30) = 2.29e-12 * (Tgas / 3d2)**(2.20) * exp(-3820.0 / Tgas)

    ! O + OH -> O2 + H
    kall(31) = 4.34e-11 * (Tgas / 3d2)**(-0.50) * exp(-30.0 / Tgas)

    ! O + H2O -> OH + OH
    kall(32) = 1.85e-11 * (Tgas / 3d2)**(0.95) * exp(-8571.0 / Tgas)

    ! CH4 + OH -> H2O + CH3
    kall(33) = 3.77e-13 * (Tgas / 3d2)**(2.42) * exp(-1162.0 / Tgas)

    ! OH + OH -> H2O + O
    kall(34) = 1.65e-12 * (Tgas / 3d2)**(1.14) * exp(-50.0 / Tgas)

    ! H + CH2+ -> CH+ + H2
    kall(35) = 1.00e-09 * exp(-7080.0 / Tgas)

    ! H + CH3+ -> CH2+ + H2
    kall(36) = 7.00e-10 * exp(-10560.0 / Tgas)

    ! H2 + He+ -> He + H+ + H
    kall(37) = 3.70e-14 * exp(-35.0 / Tgas)

    ! H2 + C+ -> CH+ + H
    kall(38) = 1.00e-10 * exp(-4640.0 / Tgas)

    ! He+ + CO -> O+ + C + He
    kall(39) = 1.40e-16 * (Tgas / 3d2)**(-0.50)

    ! H+ + H2 -> H2+ + H
    kall(40) = 1.00e-10 * exp(-21200.0 / Tgas)

    ! H + He+ -> He + H+
    kall(41) = 4.85e-15 * (Tgas / 3d2)**(0.18)

    ! H + O+ -> O + H+
    kall(42) = 5.66e-10 * (Tgas / 3d2)**(0.36) * exp(+8.6 / Tgas)

    ! H+ + O -> O+ + H
    kall(43) = 7.31e-10 * (Tgas / 3d2)**(0.23) * exp(-225.9 / Tgas)

    ! He+ + C -> C+ + He
    kall(44) = 6.30e-15 * (Tgas / 3d2)**(0.75)

    ! O+ + CO -> CO+ + O
    kall(45) = 4.90e-12 * (Tgas / 3d2)**(0.50) * exp(-4580.0 / Tgas)

    ! H2+ + E -> H + H
    kall(46) = 1.60e-08 * (Tgas / 3d2)**(-0.43)

    ! H3+ + E -> H + H + H
    kall(47) = 7.50e-08 * (Tgas / 3d2)**(-0.30)

    ! H3+ + E -> H2 + H
    kall(48) = 2.50e-08 * (Tgas / 3d2)**(-0.30)

    ! CH+ + E -> C + H
    kall(49) = 1.50e-07 * (Tgas / 3d2)**(-0.42)

    ! CH2+ + E -> CH + H
    kall(50) = 1.60e-07 * (Tgas / 3d2)**(-0.60)

    ! CH2+ + E -> C + H + H
    kall(51) = 4.03e-07 * (Tgas / 3d2)**(-0.60)

    ! CH2+ + E -> C + H2
    kall(52) = 7.68e-08 * (Tgas / 3d2)**(-0.60)

    ! CH3+ + E -> CH2 + H
    kall(53) = 1.40e-07 * (Tgas / 3d2)**(-0.50)

    ! CH3+ + E -> CH + H2
    kall(54) = 4.90e-08 * (Tgas / 3d2)**(-0.50)

    ! CH3+ + E -> CH + H + H
    kall(55) = 5.60e-08 * (Tgas / 3d2)**(-0.50)

    ! CH4+ + E -> CH3 + H
    kall(56) = 1.75e-07 * (Tgas / 3d2)**(-0.50)

    ! CH4+ + E -> CH2 + H + H
    kall(57) = 1.75e-07 * (Tgas / 3d2)**(-0.50)

    ! OH+ + E -> O + H
    kall(58) = 3.75e-08 * (Tgas / 3d2)**(-0.50)

    ! CH5+ + E -> CH3 + H2
    kall(59) = 5.50e-07 * (Tgas / 3d2)**(-0.30)

    ! CH5+ + E -> CH4 + H
    kall(60) = 5.50e-07 * (Tgas / 3d2)**(-0.30)

    ! H2O+ + E -> O + H + H
    kall(61) = 2.45e-07 * (Tgas / 3d2)**(-0.50)

    ! H2O+ + E -> O + H2
    kall(62) = 3.60e-08 * (Tgas / 3d2)**(-0.50)

    ! H2O+ + E -> OH + H
    kall(63) = 7.92e-08 * (Tgas / 3d2)**(-0.50)

    ! H3O+ + E -> O + H2 + H
    kall(64) = 5.59e-09 * (Tgas / 3d2)**(-0.50)

    ! H3O+ + E -> OH + H + H
    kall(65) = 2.58e-07 * (Tgas / 3d2)**(-0.50)

    ! H3O+ + E -> OH + H2
    kall(66) = 6.45e-08 * (Tgas / 3d2)**(-0.50)

    ! H3O+ + E -> H2O + H
    kall(67) = 1.08e-07 * (Tgas / 3d2)**(-0.50)

    ! CO+ + E -> O + C
    kall(68) = 2.00e-07 * (Tgas / 3d2)**(-0.48)

    ! HCO+ + E -> CO + H
    kall(69) = 1.10e-07 * (Tgas / 3d2)**(-1.00)

    ! O2+ + E -> O + O
    kall(70) = 1.95e-07 * (Tgas / 3d2)**(-0.70)

    ! H+ + E -> H
    kall(71) = 3.50e-12 * (Tgas / 3d2)**(-0.75)

    ! He+ + E -> He
    kall(72) = 2.36e-12 * (Tgas / 3d2)**(-0.64)

    ! C+ + E -> C
    kall(73) = 4.67e-12 * (Tgas / 3d2)**(-0.60)

    ! CH3+ + E -> CH3
    kall(74) = 1.10e-10 * (Tgas / 3d2)**(-0.50)

    ! O+ + E -> O
    kall(75) = 3.24e-12 * (Tgas / 3d2)**(-0.66)

    ! H+ + H -> H2+
    kall(76) = 5.13e-19 * (Tgas / 3d2)**(1.85)

    ! H + O -> OH
    kall(77) = 9.90e-19 * (Tgas / 3d2)**(-0.38)

    ! H + OH -> H2O
    kall(78) = 5.26e-18 * (Tgas / 3d2)**(-5.22) * exp(-90.0 / Tgas)

    ! H2 + C+ -> CH2+
    kall(79) = 4.00e-16 * (Tgas / 3d2)**(-0.20)

    ! H2 + CH -> CH3
    kall(80) = 5.09e-18 * (Tgas / 3d2)**(-0.71) * exp(-11.6 / Tgas)

    ! H2 + CH3+ -> CH5+
    kall(81) = 1.30e-14 * (Tgas / 3d2)**(-1.00)

    ! O + O -> O2
    kall(82) = 4.90e-20 * (Tgas / 3d2)**(1.58)

    ! H + H2 -> H + H + H
    kall(83) = 4.67e-07 * (Tgas / 3d2)**(-1.00) * exp(-55000.0 / Tgas)

    ! H + CH -> C + H + H
    kall(84) = 6.00e-09 * exp(-40200.0 / Tgas)

    ! H + OH -> O + H + H
    kall(85) = 6.00e-09 * exp(-50900.0 / Tgas)

    ! H + H2O -> OH + H + H
    kall(86) = 5.80e-09 * exp(-52900.0 / Tgas)

    ! H + O2 -> O + O + H
    kall(87) = 6.00e-09 * exp(-52300.0 / Tgas)

    ! H2 + E -> H + H + E
    kall(88) = 3.22e-09 * (Tgas / 3d2)**(0.35) * exp(-102000.0 / Tgas)

    ! H2 + H2 -> H2 + H + H
    kall(89) = 1.00e-08 * exp(-84100.0 / Tgas)

    ! H2 + CH -> C + H2 + H
    kall(90) = 6.00e-09 * exp(-40200.0 / Tgas)

    ! H2 + OH -> O + H2 + H
    kall(91) = 6.00e-09 * exp(-50900.0 / Tgas)

    ! H2 + H2O -> OH + H2 + H
    kall(92) = 5.80e-09 * exp(-52900.0 / Tgas)

    ! H2 + O2 -> O + O + H2
    kall(93) = 6.00e-09 * exp(-52300.0 / Tgas)

    ! CH + O -> HCO+ + E
    kall(94) = 2.00e-11 * (Tgas / 3d2)**(0.44)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_RATES
  end subroutine evaluate_interpolable

end module prizmo_rates_interpolable
