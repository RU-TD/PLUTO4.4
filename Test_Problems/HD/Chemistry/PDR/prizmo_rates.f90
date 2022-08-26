module prizmo_rates
  use prizmo_commons
  use prizmo_shielding
  use prizmo_utils
contains

  ! ************************
  subroutine compute_rates(x, Tgas, Tdust)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust
    real*8::te, lnte, Tgas32, log_tgas, invTgas, sqrTgas32, invsqrTgas32, invsqrTgas
    real*8::invTgas32, Tgas14, sqrTgas, invTe, invsqrTe
    real*8::Tgas32_m06, Tgas32_m05, Tgas32_m03
    real*8::sticking, nu_debye

    ! temperature shortcuts
    te = Tgas * 8.617343d-5
    lnte = log(te)
    invTe = 1d0 / te
    invsqrTe = 1d0 / sqrt(te)
    Tgas32 = Tgas / 3d2
    Tgas14 = Tgas / 1d4
    invTgas32 = 1d0 / Tgas32
    log_tgas = log10(Tgas)
    invTgas = 1d0 / Tgas
    sqrTgas32 = sqrt(Tgas32)
    sqrTgas = sqrt(Tgas)
    invsqrTgas32 = 1d0 / sqrTgas32
    invsqrTgas = 1d0 / sqrTgas
    Tgas32_m06 = Tgas32**(-0.6)
    Tgas32_m05 = Tgas32**(-0.5)
    Tgas32_m03 = Tgas32**(-0.3)

    ! dust parameters
    sticking = 1d0
    nu_debye = 1d12  ! 1/s

    !! PREPROCESS_RATES

    ! H + CH -> C + H2
    kall(1) = 2.70e-11 * Tgas32**(0.38)
    
    ! H + CH2 -> CH + H2
    kall(2) = 6.64e-11
    
    ! H + CH3 -> CH2 + H2
    kall(3) = 1.00e-10 * exp(-7600.0*invTgas)
    
    ! H + CH4 -> CH3 + H2
    kall(4) = 5.94e-13 * Tgas32**(3.00) * exp(-4045.0*invTgas)
    
    ! H + OH -> O + H2
    kall(5) = 6.99e-14 * Tgas32**(2.80) * exp(-1950.0*invTgas)
    
    ! H + H2O -> OH + H2
    kall(6) = 1.59e-11 * Tgas32**(1.20) * exp(-9610.0*invTgas)
    
    ! H + CO -> OH + C
    kall(7) = 1.10e-10 * sqrTgas32 * exp(-77700.0*invTgas)
    
    ! H + O2 -> OH + O
    kall(8) = 2.61e-10 * exp(-8156.0*invTgas)
    
    ! H2 + C -> CH + H
    kall(9) = 6.64e-10 * exp(-11700.0*invTgas)
    
    ! H2 + CH -> CH2 + H
    kall(10) = 5.46e-10 * exp(-1943.0*invTgas)
    
    ! H2 + CH2 -> CH3 + H
    kall(11) = 5.18e-11 * Tgas32**(0.17) * exp(-6400.0*invTgas)
    
    ! H2 + CH3 -> CH4 + H
    kall(12) = 6.86e-14 * Tgas32**(2.74) * exp(-4740.0*invTgas)
    
    ! H2 + O -> OH + H
    kall(13) = 3.14e-13 * Tgas32**(2.70) * exp(-3150.0*invTgas)
    
    ! H2 + OH -> H2O + H
    kall(14) = 2.05e-12 * Tgas32**(1.52) * exp(-1736.0*invTgas)
    
    ! H2 + O2 -> OH + OH
    kall(15) = 3.16e-10 * exp(-21890.0*invTgas)
    
    ! C + CH2 -> CH + CH
    kall(16) = 2.69e-12 * exp(-23550.0*invTgas)
    
    ! C + OH -> O + CH
    kall(17) = 2.25e-11 * sqrTgas32 * exp(-14800.0*invTgas)
    
    ! C + OH -> CO + H
    kall(18) = 1.10e-10 * sqrTgas32
    
    ! C + O2 -> CO + O
    kall(19) = 3.30e-11
    
    ! CH + O -> OH + C
    kall(20) = 2.52e-11 * exp(-2381.0*invTgas)
    
    ! CH + O -> CO + H
    kall(21) = 6.60e-11
    
    ! CH + CH4 -> CH3 + CH2
    kall(22) = 2.28e-11 * Tgas32**(0.70) * exp(-3000.0*invTgas)
    
    ! CH + O2 -> CO + OH
    kall(23) = 2.60e-11
    
    ! CH2 + CH2 -> CH3 + CH
    kall(24) = 4.00e-10 * exp(-5000.0*invTgas)
    
    ! CH2 + O -> OH + CH
    kall(25) = 4.98e-10 * exp(-6000.0*invTgas)
    
    ! CH2 + O -> CO + H + H
    kall(26) = 1.33e-10
    
    ! CH2 + O -> CO + H2
    kall(27) = 8.00e-11
    
    ! CH2 + CH4 -> CH3 + CH3
    kall(28) = 7.13e-12 * exp(-5050.0*invTgas)
    
    ! CH2 + OH -> O + CH3
    kall(29) = 1.44e-11 * sqrTgas32 * exp(-3000.0*invTgas)
    
    ! CH2 + OH -> H2O + CH
    kall(30) = 1.44e-11 * sqrTgas32 * exp(-3000.0*invTgas)
    
    ! CH2 + O2 -> CO + H2O
    kall(31) = 2.48e-10 * Tgas32**(-3.30) * exp(-1443.0*invTgas)
    
    ! CH3 + CH3 -> CH4 + CH2
    kall(32) = 7.13e-12 * exp(-5052.0*invTgas)
    
    ! CH3 + OH -> CH4 + O
    kall(33) = 3.27e-14 * Tgas32**(2.20) * exp(-2240.0*invTgas)
    
    ! CH3 + OH -> H2O + CH2
    kall(34) = 1.20e-10 * exp(-1400.0*invTgas)
    
    ! CH3 + H2O -> OH + CH4
    kall(35) = 2.30e-15 * Tgas32**(3.47) * exp(-6681.0*invTgas)
    
    ! O + CH4 -> OH + CH3
    kall(36) = 2.29e-12 * Tgas32**(2.20) * exp(-3820.0*invTgas)
    
    ! O + OH -> O2 + H
    kall(37) = 4.34e-11 * Tgas32_m05 * exp(-30.0*invTgas)
    
    ! O + H2O -> OH + OH
    kall(38) = 1.85e-11 * Tgas32**(0.95) * exp(-8571.0*invTgas)
    
    ! CH4 + OH -> H2O + CH3
    kall(39) = 3.77e-13 * Tgas32**(2.42) * exp(-1162.0*invTgas)
    
    ! OH + OH -> H2O + O
    kall(40) = 1.65e-12 * Tgas32**(1.14) * exp(-50.0*invTgas)
    
    ! H + CH+ -> C+ + H2
    kall(41) = 7.50e-10
    
    ! H + CH2+ -> CH+ + H2
    kall(42) = 1.00e-09 * exp(-7080.0*invTgas)
    
    ! H + CH3+ -> CH2+ + H2
    kall(43) = 7.00e-10 * exp(-10560.0*invTgas)
    
    ! H2 + He+ -> He + H+ + H
    kall(44) = 3.70e-14 * exp(-35.0*invTgas)
    
    ! H2 + C+ -> CH+ + H
    kall(45) = 1.00e-10 * exp(-4640.0*invTgas)
    
    ! H2 + CH+ -> CH2+ + H
    kall(46) = 1.20e-09
    
    ! He+ + CO -> O+ + C + He
    kall(47) = 1.40e-16 * Tgas32_m05
    
    ! C+ + OH -> CO + H+
    kall(48) = 7.70e-10
    
    ! H+ + H2 -> H2+ + H
    kall(49) = 1.00e-10 * exp(-21200.0*invTgas)
    
    ! H + He+ -> He + H+
    kall(50) = 4.85e-15 * Tgas32**(0.18)
    
    ! H + O+ -> O + H+
    kall(51) = 5.66e-10 * Tgas32**(0.36) * exp(+8.6*invTgas)
    
    ! H+ + O -> O+ + H
    kall(52) = 7.31e-10 * Tgas32**(0.23) * exp(-225.9*invTgas)
    
    ! H+ + H2O -> H2O+ + H
    kall(53) = 6.90e-09
    
    ! H+ + O2 -> O2+ + H
    kall(54) = 2.00e-09
    
    ! H2 + He+ -> He + H2+
    kall(55) = 7.20e-15
    
    ! He+ + C -> C+ + He
    kall(56) = 6.30e-15 * Tgas32**(0.75)
    
    ! O+ + CO -> CO+ + O
    kall(57) = 4.90e-12 * sqrTgas32 * exp(-4580.0*invTgas)
    
    ! H2+ + E -> H + H
    kall(58) = 1.60e-08 * Tgas32**(-0.43)
    
    ! H3+ + E -> H + H + H
    kall(59) = 7.50e-08 * Tgas32_m03
    
    ! H3+ + E -> H2 + H
    kall(60) = 2.50e-08 * Tgas32_m03
    
    ! CH+ + E -> C + H
    kall(61) = 1.50e-07 * Tgas32**(-0.42)
    
    ! CH2+ + E -> CH + H
    kall(62) = 1.60e-07 * Tgas32_m06
    
    ! CH2+ + E -> C + H + H
    kall(63) = 4.03e-07 * Tgas32_m06
    
    ! CH2+ + E -> C + H2
    kall(64) = 7.68e-08 * Tgas32_m06
    
    ! CH3+ + E -> CH2 + H
    kall(65) = 1.40e-07 * Tgas32_m05
    
    ! CH3+ + E -> CH + H2
    kall(66) = 4.90e-08 * Tgas32_m05
    
    ! CH3+ + E -> CH + H + H
    kall(67) = 5.60e-08 * Tgas32_m05
    
    ! CH3+ + E -> C + H2 + H
    kall(68) = 1.05e-07 * Tgas32**(-0.5)
    
    ! CH4+ + E -> CH3 + H
    kall(69) = 1.75e-07 * Tgas32_m05
    
    ! CH4+ + E -> CH2 + H + H
    kall(70) = 1.75e-07 * Tgas32_m05
    
    ! OH+ + E -> O + H
    kall(71) = 3.75e-08 * Tgas32_m05
    
    ! CH5+ + E -> CH3 + H2
    kall(72) = 5.50e-07 * Tgas32_m03
    
    ! CH5+ + E -> CH4 + H
    kall(73) = 5.50e-07 * Tgas32_m03
    
    ! H2O+ + E -> O + H + H
    kall(74) = 2.45e-07 * Tgas32_m05
    
    ! H2O+ + E -> O + H2
    kall(75) = 3.60e-08 * Tgas32_m05
    
    ! H2O+ + E -> OH + H
    kall(76) = 7.92e-08 * Tgas32_m05
    
    ! H3O+ + E -> O + H2 + H
    kall(77) = 5.59e-09 * Tgas32_m05
    
    ! H3O+ + E -> OH + H + H
    kall(78) = 2.58e-07 * Tgas32_m05
    
    ! H3O+ + E -> OH + H2
    kall(79) = 6.45e-08 * Tgas32_m05
    
    ! H3O+ + E -> H2O + H
    kall(80) = 1.08e-07 * Tgas32_m05
    
    ! CO+ + E -> O + C
    kall(81) = 2.00e-07 * Tgas32**(-0.48)
    
    ! HCO+ + E -> CO + H
    kall(82) = 1.10e-07 * invTgas32
    
    ! O2+ + E -> O + O
    kall(83) = 1.95e-07 * Tgas32**(-0.70)
    
    ! CH3+ + E -> CH3
    kall(84) = 1.10e-10 * Tgas32_m05
    
    ! H+ + H -> H2+
    kall(85) = 5.13e-19 * Tgas32**(1.85)
    
    ! H + C -> CH
    kall(86) = 1.00e-17
    
    ! H + C+ -> CH+
    kall(87) = 1.70e-17
    
    ! H + O -> OH
    kall(88) = 9.90e-19 * Tgas32**(-0.38)
    
    ! H + OH -> H2O
    kall(89) = 5.26e-18 * Tgas32**(-5.22) * exp(-90.0*invTgas)
    
    ! H2 + C -> CH2
    kall(90) = 1.00e-17
    
    ! H2 + C+ -> CH2+
    kall(91) = 4.00e-16 * Tgas32**(-0.20)
    
    ! H2 + CH -> CH3
    kall(92) = 5.09e-18 * Tgas32**(-0.71) * exp(-11.6*invTgas)
    
    ! H2 + CH3+ -> CH5+
    kall(93) = 1.30e-14 * invTgas32
    
    ! C + O -> CO
    kall(94) = 2.10e-19
    
    ! C+ + O -> CO+
    kall(95) = 2.50e-18
    
    ! O + O -> O2
    kall(96) = 4.90e-20 * Tgas32**(1.58)
    
    ! H + H2 -> H + H + H
    kall(97) = 4.67e-07 * invTgas32 * exp(-55000.0*invTgas)
    
    ! H + CH -> C + H + H
    kall(98) = 6.00e-09 * exp(-40200.0*invTgas)
    
    ! H + OH -> O + H + H
    kall(99) = 6.00e-09 * exp(-50900.0*invTgas)
    
    ! H + H2O -> OH + H + H
    kall(100) = 5.80e-09 * exp(-52900.0*invTgas)
    
    ! H + O2 -> O + O + H
    kall(101) = 6.00e-09 * exp(-52300.0*invTgas)
    
    ! H2 + E -> H + H + E
    kall(102) = 3.22e-09 * Tgas32**(0.35) * exp(-102000.0*invTgas)
    
    ! H2 + H2 -> H2 + H + H
    kall(103) = 1.00e-08 * exp(-84100.0*invTgas)
    
    ! H2 + CH -> C + H2 + H
    kall(104) = 6.00e-09 * exp(-40200.0*invTgas)
    
    ! H2 + OH -> O + H2 + H
    kall(105) = 6.00e-09 * exp(-50900.0*invTgas)
    
    ! H2 + H2O -> OH + H2 + H
    kall(106) = 5.80e-09 * exp(-52900.0*invTgas)
    
    ! H2 + O2 -> O + O + H2
    kall(107) = 6.00e-09 * exp(-52300.0*invTgas)
    
    ! CH + O -> HCO+ + E
    kall(108) = 2.00e-11 * Tgas32**(0.44)
    
    ! H + H -> H2
    kall(109) = 2.121e-17 * d2g / 1d-2
    
    ! H+ + CH2 -> CH+ + H2
    kall(110) = 1.40e-09
    
    ! H+ + CH4 -> CH3+ + H2
    kall(111) = 2.30e-09
    
    ! H + CH4+ -> CH3+ + H2
    kall(112) = 1.00e-11
    
    ! H + CH5+ -> CH4+ + H2
    kall(113) = 2.00e-11
    
    ! H2+ + H2 -> H3+ + H
    kall(114) = 2.08e-09
    
    ! H2+ + C -> CH+ + H
    kall(115) = 2.40e-09
    
    ! H2+ + CH -> CH2+ + H
    kall(116) = 7.10e-10
    
    ! H2+ + CH2 -> CH3+ + H
    kall(117) = 1.00e-09
    
    ! H2 + CH2+ -> CH3+ + H
    kall(118) = 1.60e-09
    
    ! H2+ + O -> OH+ + H
    kall(119) = 1.50e-09
    
    ! H2 + O+ -> OH+ + H
    kall(120) = 1.70e-09
    
    ! H2+ + CH4 -> CH5+ + H
    kall(121) = 1.14e-10
    
    ! H2 + CH4+ -> CH5+ + H
    kall(122) = 3.30e-11
    
    ! H2+ + CH4 -> CH3+ + H2 + H
    kall(123) = 2.30e-09
    
    ! H2+ + OH -> H2O+ + H
    kall(124) = 7.60e-10
    
    ! H2 + OH+ -> H2O+ + H
    kall(125) = 1.01e-09
    
    ! H2+ + H2O -> H3O+ + H
    kall(126) = 3.40e-09
    
    ! H2 + H2O+ -> H3O+ + H
    kall(127) = 6.40e-10
    
    ! H2+ + CO -> HCO+ + H
    kall(128) = 2.16e-09
    
    ! H2 + CO+ -> HCO+ + H
    kall(129) = 1.80e-09
    
    ! H3+ + C -> CH+ + H2
    kall(130) = 2.00e-09
    
    ! H3+ + CH -> CH2+ + H2
    kall(131) = 1.20e-09
    
    ! H3+ + CH2 -> CH3+ + H2
    kall(132) = 1.70e-09
    
    ! H3+ + CH3 -> CH4+ + H2
    kall(133) = 2.10e-09
    
    ! H3+ + O -> OH+ + H2
    kall(134) = 8.00e-10
    
    ! H3+ + CH4 -> CH5+ + H2
    kall(135) = 2.40e-09
    
    ! H3+ + OH -> H2O+ + H2
    kall(136) = 1.30e-09
    
    ! H3+ + H2O -> H3O+ + H2
    kall(137) = 5.90e-09
    
    ! H3+ + CO -> HCO+ + H2
    kall(138) = 1.70e-09
    
    ! He+ + CH -> C+ + He + H
    kall(139) = 1.10e-09
    
    ! He+ + CH2 -> C+ + He + H2
    kall(140) = 7.50e-10
    
    ! He+ + CH2 -> CH+ + He + H
    kall(141) = 7.50e-10
    
    ! He+ + CH3 -> CH+ + He + H2
    kall(142) = 1.80e-09
    
    ! He+ + CH4 -> CH+ + He + H2 + H
    kall(143) = 2.40e-10
    
    ! He+ + CH4 -> CH2+ + He + H2
    kall(144) = 9.50e-10
    
    ! He+ + CH4 -> CH3 + He + H+
    kall(145) = 4.80e-10
    
    ! He+ + CH4 -> CH3+ + He + H
    kall(146) = 8.50e-11
    
    ! He+ + OH -> O+ + He + H
    kall(147) = 1.10e-09
    
    ! He+ + H2O -> OH + He + H+
    kall(148) = 2.04e-10
    
    ! He+ + H2O -> OH+ + He + H
    kall(149) = 2.86e-10
    
    ! He+ + CO -> O + C+ + He
    kall(150) = 1.60e-09
    
    ! He+ + O2 -> O+ + O + He
    kall(151) = 1.00e-09
    
    ! C + OH+ -> O + CH+
    kall(152) = 1.20e-09
    
    ! C+ + OH -> CO+ + H
    kall(153) = 7.70e-10
    
    ! C + CH5+ -> CH4 + CH+
    kall(154) = 1.20e-09
    
    ! C + H2O+ -> OH + CH+
    kall(155) = 1.10e-09
    
    ! C+ + H2O -> HCO+ + H
    kall(156) = 9.00e-10
    
    ! C + H3O+ -> HCO+ + H2
    kall(157) = 1.00e-11
    
    ! C + HCO+ -> CO + CH+
    kall(158) = 1.10e-09
    
    ! C+ + O2 -> CO+ + O
    kall(159) = 3.80e-10
    
    ! C+ + O2 -> CO + O+
    kall(160) = 6.20e-10
    
    ! C + O2+ -> CO+ + O
    kall(161) = 5.20e-11
    
    ! CH+ + O -> CO+ + H
    kall(162) = 3.50e-10
    
    ! CH + O+ -> CO+ + H
    kall(163) = 3.50e-10
    
    ! CH + OH+ -> O + CH2+
    kall(164) = 3.50e-10
    
    ! CH+ + OH -> CO+ + H2
    kall(165) = 7.50e-10
    
    ! CH + CH5+ -> CH4 + CH2+
    kall(166) = 6.90e-10
    
    ! CH+ + H2O -> H3O+ + C
    kall(167) = 5.80e-10
    
    ! CH + H2O+ -> OH + CH2+
    kall(168) = 3.40e-10
    
    ! CH+ + H2O -> HCO+ + H2
    kall(169) = 2.90e-09
    
    ! CH + H3O+ -> H2O + CH2+
    kall(170) = 6.80e-10
    
    ! CH + CO+ -> HCO+ + C
    kall(171) = 3.20e-10
    
    ! CH + HCO+ -> CO + CH2+
    kall(172) = 6.30e-10
    
    ! CH+ + O2 -> CO+ + OH
    kall(173) = 1.00e-11
    
    ! CH+ + O2 -> HCO+ + O
    kall(174) = 9.70e-10
    
    ! CH + O2+ -> HCO+ + O
    kall(175) = 3.10e-10
    
    ! CH2+ + O -> HCO+ + H
    kall(176) = 7.50e-10
    
    ! CH2 + OH+ -> O + CH3+
    kall(177) = 4.80e-10
    
    ! CH2 + CH5+ -> CH4 + CH3+
    kall(178) = 9.60e-10
    
    ! CH2 + H2O+ -> OH + CH3+
    kall(179) = 4.70e-10
    
    ! CH2 + H3O+ -> H2O + CH3+
    kall(180) = 9.40e-10
    
    ! CH2 + CO+ -> HCO+ + CH
    kall(181) = 4.30e-10
    
    ! CH2 + HCO+ -> CO + CH3+
    kall(182) = 8.60e-10
    
    ! CH2+ + O2 -> HCO+ + OH
    kall(183) = 9.10e-10
    
    ! CH3+ + O -> HCO+ + H2
    kall(184) = 4.00e-10
    
    ! O+ + CH4 -> OH + CH3+
    kall(185) = 1.10e-10
    
    ! O + CH4+ -> OH + CH3+
    kall(186) = 1.00e-09
    
    ! O+ + OH -> O2+ + H
    kall(187) = 3.60e-10
    
    ! O + OH+ -> O2+ + H
    kall(188) = 7.10e-10
    
    ! O + CH5+ -> H3O+ + CH2
    kall(189) = 2.20e-10
    
    ! O + H2O+ -> O2+ + H2
    kall(190) = 4.00e-11
    
    ! CH4+ + CH4 -> CH5+ + CH3
    kall(191) = 1.50e-09
    
    ! CH4 + OH+ -> CH5+ + O
    kall(192) = 1.95e-10
    
    ! CH4 + OH+ -> H3O+ + CH2
    kall(193) = 1.31e-09
    
    ! CH4+ + H2O -> H3O+ + CH3
    kall(194) = 2.60e-09
    
    ! CH4 + H2O+ -> H3O+ + CH3
    kall(195) = 1.40e-09
    
    ! CH4+ + CO -> HCO+ + CH3
    kall(196) = 1.40e-09
    
    ! CH4 + CO+ -> HCO+ + CH3
    kall(197) = 4.55e-10
    
    ! OH+ + OH -> H2O+ + O
    kall(198) = 7.00e-10
    
    ! OH + CH5+ -> H2O+ + CH4
    kall(199) = 7.00e-10
    
    ! OH+ + H2O -> H3O+ + O
    kall(200) = 1.30e-09
    
    ! OH + H2O+ -> H3O+ + O
    kall(201) = 6.90e-10
    
    ! OH+ + CO -> HCO+ + O
    kall(202) = 1.05e-09
    
    ! OH + CO+ -> HCO+ + O
    kall(203) = 3.10e-10
    
    ! OH + HCO+ -> CO + H2O+
    kall(204) = 6.20e-10
    
    ! CH5+ + H2O -> H3O+ + CH4
    kall(205) = 3.70e-09
    
    ! CH5+ + CO -> HCO+ + CH4
    kall(206) = 1.00e-09
    
    ! H2O+ + H2O -> H3O+ + OH
    kall(207) = 2.10e-09
    
    ! H2O+ + CO -> HCO+ + OH
    kall(208) = 5.00e-10
    
    ! H2O + CO+ -> HCO+ + OH
    kall(209) = 8.84e-10
    
    ! H2O + HCO+ -> CO + H3O+
    kall(210) = 2.50e-09
    
    ! H + H2+ -> H2 + H+
    kall(211) = 6.40e-10
    
    ! H+ + CH -> CH+ + H
    kall(212) = 1.90e-09
    
    ! H+ + CH2 -> CH2+ + H
    kall(213) = 1.40e-09
    
    ! H+ + CH3 -> CH3+ + H
    kall(214) = 3.40e-09
    
    ! H+ + CH4 -> CH4+ + H
    kall(215) = 1.50e-09
    
    ! H+ + OH -> OH+ + H
    kall(216) = 2.10e-09
    
    ! H + CO+ -> CO + H+
    kall(217) = 7.50e-10
    
    ! H2+ + CH -> CH+ + H2
    kall(218) = 7.10e-10
    
    ! H2+ + CH2 -> CH2+ + H2
    kall(219) = 1.00e-09
    
    ! H2+ + CH4 -> CH4+ + H2
    kall(220) = 1.40e-09
    
    ! H2+ + OH -> OH+ + H2
    kall(221) = 7.60e-10
    
    ! H2+ + H2O -> H2O+ + H2
    kall(222) = 3.90e-09
    
    ! H2+ + CO -> CO+ + H2
    kall(223) = 6.40e-10
    
    ! H2+ + O2 -> O2+ + H2
    kall(224) = 8.00e-10
    
    ! He+ + CH -> CH+ + He
    kall(225) = 5.00e-10
    
    ! He+ + CH4 -> CH4+ + He
    kall(226) = 5.10e-11
    
    ! He+ + H2O -> H2O+ + He
    kall(227) = 6.05e-11
    
    ! He+ + O2 -> O2+ + He
    kall(228) = 3.30e-11
    
    ! C+ + CH -> CH+ + C
    kall(229) = 3.80e-10
    
    ! C+ + CH2 -> CH2+ + C
    kall(230) = 5.20e-10
    
    ! C + CO+ -> CO + C+
    kall(231) = 1.10e-10
    
    ! C + O2+ -> O2 + C+
    kall(232) = 5.20e-11
    
    ! CH + O+ -> O + CH+
    kall(233) = 3.50e-10
    
    ! CH + OH+ -> OH + CH+
    kall(234) = 3.50e-10
    
    ! CH + H2O+ -> H2O + CH+
    kall(235) = 3.40e-10
    
    ! CH + CO+ -> CO + CH+
    kall(236) = 3.20e-10
    
    ! CH + O2+ -> O2 + CH+
    kall(237) = 3.10e-10
    
    ! CH2 + O+ -> O + CH2+
    kall(238) = 9.70e-10
    
    ! CH2 + OH+ -> OH + CH2+
    kall(239) = 4.80e-10
    
    ! CH2 + H2O+ -> H2O + CH2+
    kall(240) = 4.70e-10
    
    ! CH2 + CO+ -> CO + CH2+
    kall(241) = 4.30e-10
    
    ! CH2 + O2+ -> O2 + CH2+
    kall(242) = 4.30e-10
    
    ! O+ + CH4 -> CH4+ + O
    kall(243) = 8.90e-10
    
    ! O+ + OH -> OH+ + O
    kall(244) = 3.60e-10
    
    ! O+ + H2O -> H2O+ + O
    kall(245) = 3.20e-09
    
    ! O + CO+ -> CO + O+
    kall(246) = 1.40e-10
    
    ! O+ + O2 -> O2+ + O
    kall(247) = 1.90e-11
    
    ! CH4 + CO+ -> CO + CH4+
    kall(248) = 7.93e-10
    
    ! CH4+ + O2 -> O2+ + CH4
    kall(249) = 4.00e-10
    
    ! OH+ + H2O -> H2O+ + OH
    kall(250) = 1.59e-09
    
    ! OH + CO+ -> CO + OH+
    kall(251) = 3.10e-10
    
    ! OH+ + O2 -> O2+ + OH
    kall(252) = 5.90e-10
    
    ! H2O + CO+ -> CO + H2O+
    kall(253) = 1.72e-09
    
    ! H2O+ + O2 -> O2+ + H2O
    kall(254) = 4.60e-10
    
    ! CO+ + O2 -> O2+ + CO
    kall(255) = 1.20e-10
    
    ! H+ + E -> H
    kall(256) = 1.416215e-10 * invsqrTgas / ((1e0 + 5.636151e-01 * sqrTgas)**(0.25) * (1e0 + 1.192167e-03 * sqrTgas)**(1.75))
    
    ! He+ + E -> He
    kall(257) = 1.298521e-10 * invsqrTgas / ((1e0 + 2.536731e-01 * sqrTgas)**(0.31) * (1e0 + 1.649348e-04 * sqrTgas)**(1.69))
    
    ! He+ + E -> He
    kall(258) = 1.932416e-10 * invsqrTgas / ((1e0 + 4.841607e+00 * sqrTgas)**(0.21) * (1e0 + 4.623984e-04 * sqrTgas)**(1.79))
    
    ! He++ + E -> He+
    kall(259) = 5.788437e-10 * invsqrTgas / ((1e0 + 3.266858e-01 * sqrTgas)**(0.25) * (1e0 + 6.004084e-04 * sqrTgas)**(1.75))
    
    ! C+ + E -> C
    kall(260) = 4.700000e-13 / Tgas14**(0.62)
    
    ! C++ + E -> C+
    kall(261) = 2.300000e-12 / Tgas14**(0.65)
    
    ! C+++ + E -> C++
    kall(262) = 4.900000e-12 / Tgas14**(0.80)
    
    ! C++++ + E -> C+++
    kall(263) = 1.912274e-09 * invsqrTgas / ((1e0 + 4.465888e-02 * sqrTgas)**(0.48) * (1e0 + 2.600255e-04 * sqrTgas)**(1.52))
    
    ! O+ + E -> O
    kall(264) = 3.100000e-13 / Tgas14**(0.68)
    
    ! O++ + E -> O+
    kall(265) = 2.000000e-12 / Tgas14**(0.65)
    
    ! O+++ + E -> O++
    kall(266) = 5.100000e-12 / Tgas14**(0.66)
    
    ! O++++ + E -> O+++
    kall(267) = 9.600000e-12 / Tgas14**(0.67)
    
    ! H + E -> H+ + E + E
    kall(268) = 2.910000e-08 * 1e0 / (2.320000e-01 + 1.360000e+01 * invTe) * 2.767452e+00 * invTe**(0.39) * exp(-13.60 * invTe)
    
    ! He + E -> He+ + E + E
    kall(269) = 1.750000e-08 * 1e0 / (1.800000e-01 + 2.460000e+01 * invTe) * 3.067802e+00 * invTe**(0.35) * exp(-24.60 * invTe)
    
    ! He+ + E -> He++ + E + E
    kall(270) = 2.050000e-09 * (1e0 + 7.375636e+00*invsqrTe) / (2.650000e-01 + 5.440000e+01 * invTe) * 2.715812e+00 * invTe**(0.25) * exp(-54.40 * invTe)
    
    ! C + E -> C+ + E + E
    kall(271) = 6.850000e-08 * 1e0 / (1.930000e-01 + 1.130000e+01 * invTe) * 1.833452e+00 * invTe**(0.25) * exp(-11.30 * invTe)
    
    ! C+ + E -> C++ + E + E
    kall(272) = 1.860000e-08 * (1e0 + 4.939636e+00 * invsqrTe) / (2.860000e-01 + 2.440000e+01 * invTe) * 2.152651e+00 * invTe**(0.24) * exp(-24.40 * invTe)
    
    ! C++ + E -> C+++ + E + E
    kall(273) = 6.350000e-09 * (1e0 + 6.920983e+00 * invsqrTe) / (4.270000e-01 + 4.790000e+01 * invTe) * 2.253567e+00 * invTe**(0.21) * exp(-47.90 * invTe)
    
    ! C+++ + E -> C++++ + E + E
    kall(274) = 1.500000e-09 * (1e0 + 8.031189e+00 * invsqrTe) / (4.160000e-01 + 6.450000e+01 * invTe) * 1.718869e+00 * invTe**(0.13) * exp(-64.50 * invTe)
    
    ! O + E -> O+ + E + E
    kall(275) = 3.590000e-08 * 1e0 / (7.300000e-02 + 1.360000e+01 * invTe) * 2.428864e+00 * invTe**(0.34) * exp(-13.60 * invTe)
    
    ! O+ + E -> O++ + E + E
    kall(276) = 1.390000e-08 * (1e0 + 5.924525e+00 * invsqrTe) / (2.120000e-01 + 3.510000e+01 * invTe) * 2.187598e+00 * invTe**(0.22) * exp(-35.10 * invTe)
    
    ! O++ + E -> O+++ + E + E
    kall(277) = 9.310000e-09 * (1e0 + 7.409453e+00 * invsqrTe) / (2.700000e-01 + 5.490000e+01 * invTe) * 2.949066e+00 * invTe**(0.27) * exp(-54.90 * invTe)
    
    ! O+++ + E -> O++++ + E + E
    kall(278) = 1.020000e-08 * 1e0 / (6.140000e-01 + 7.740000e+01 * invTe) * 3.235639e+00 * invTe**(0.27) * exp(-77.40 * invTe)
    
    ! CO -> CO_DUST
    kall(279) = 1.325861e-03 * rho_dust * sticking * sqrTgas * 1.461242e+11 
    
    ! CO_DUST -> CO
    kall(280) = nu_debye * exp(-1.390000e+03 / Tdust)
    
    ! H2O -> H2O_DUST
    kall(281) = 1.325861e-03 * rho_dust * sticking * sqrTgas * 1.728965e+11 
    
    ! H2O_DUST -> H2O
    kall(282) = nu_debye * exp(-5.700000e+03 / Tdust)
    
    ! kall(283) = PHOTO, 13.60
    
    ! kall(284) = PHOTO, 24.59
    
    ! kall(285) = PHOTO, 54.42
    
    ! kall(286) = PHOTO, 54.49
    
    ! kall(287) = PHOTO, 77.41
    
    ! kall(288) = PHOTO, 113.9
    
    ! kall(289) = PHOTO, 138.1
    
    ! kall(290) = PHOTO, 11.2603
    
    ! kall(291) = PHOTO, 24.38
    
    ! kall(292) = PHOTO, 47.89
    
    ! kall(293) = PHOTO, 64.49
    
    ! kall(294) = PHOTO, 1e99
    
    ! kall(295) = PHOTO, 1e99
    
    ! kall(296) = PHOTO, 1e99
    
    ! kall(297) = PHOTO, 1e99
    
    ! kall(298) = PHOTO, 1e99
    
    ! kall(299) = PHOTO, 1e99
    
    ! kall(300) = PHOTO, 1e99
    
    ! kall(301) = PHOTO, 1e99
    
    ! kall(302) = PHOTO, 1e99
    
    ! kall(303) = PHOTO, 1e99
    
    ! kall(304) = PHOTO, 1e99
    
    ! kall(305) = PHOTO, 1e99
    
    ! kall(306) = PHOTO, 1e99
    
    ! kall(307) = PHOTO, 1e99
    
    ! kall(308) = PHOTO, 1e99
    
    ! kall(309) = PHOTO, 1e99
    
    ! kall(310) = PHOTO, 1e99
    
    ! H2 -> H+ + H + E
    kall(311) = 0.02e+00 * user_cr
    
    ! H2 -> H + H
    kall(312) = 0.10e+00 * user_cr
    
    ! H2 -> H2+ + E
    kall(313) = 0.88e+00 * user_cr
    
    ! H -> H+ + E
    kall(314) = 0.46e+00 * user_cr
    
    ! He -> He+ + E
    kall(315) = 0.50e+00 * user_cr

    !! PREPROCESS_END

  end subroutine compute_rates

end module prizmo_rates
