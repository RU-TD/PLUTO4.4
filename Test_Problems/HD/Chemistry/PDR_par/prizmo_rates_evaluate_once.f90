module prizmo_rates_evaluate_once
contains
  ! *****************
  subroutine init_evaluate_once()
    use prizmo_commons
    implicit none

    !!BEGIN_INIT_EVALUATE_ONCE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! H + CH2 -> CH + H2
    kall(95) = 6.64e-11

    ! C + O2 -> CO + O
    kall(96) = 3.30e-11

    ! CH + O -> CO + H
    kall(97) = 6.60e-11

    ! CH + O2 -> CO + OH
    kall(98) = 2.60e-11

    ! CH2 + O -> CO + H + H
    kall(99) = 1.33e-10

    ! CH2 + O -> CO + H2
    kall(100) = 8.00e-11

    ! H + CH+ -> C+ + H2
    kall(101) = 7.50e-10

    ! H2 + CH+ -> CH2+ + H
    kall(102) = 1.20e-09

    ! C+ + OH -> CO + H+
    kall(103) = 7.70e-10

    ! H+ + H2O -> H2O+ + H
    kall(104) = 6.90e-09

    ! H+ + O2 -> O2+ + H
    kall(105) = 2.00e-09

    ! H2 + He+ -> He + H2+
    kall(106) = 7.20e-15

    ! H + C -> CH
    kall(107) = 1.00e-17

    ! H + C+ -> CH+
    kall(108) = 1.70e-17

    ! H2 + C -> CH2
    kall(109) = 1.00e-17

    ! C + O -> CO
    kall(110) = 2.10e-19

    ! C+ + O -> CO+
    kall(111) = 2.50e-18

    ! H+ + CH2 -> CH+ + H2
    kall(112) = 1.40e-09

    ! H+ + CH4 -> CH3+ + H2
    kall(113) = 2.30e-09

    ! H + CH4+ -> CH3+ + H2
    kall(114) = 1.00e-11

    ! H + CH5+ -> CH4+ + H2
    kall(115) = 2.00e-11

    ! H2+ + H2 -> H3+ + H
    kall(116) = 2.08e-09

    ! H2+ + C -> CH+ + H
    kall(117) = 2.40e-09

    ! H2+ + CH -> CH2+ + H
    kall(118) = 7.10e-10

    ! H2+ + CH2 -> CH3+ + H
    kall(119) = 1.00e-09

    ! H2 + CH2+ -> CH3+ + H
    kall(120) = 1.60e-09

    ! H2+ + O -> OH+ + H
    kall(121) = 1.50e-09

    ! H2 + O+ -> OH+ + H
    kall(122) = 1.70e-09

    ! H2+ + CH4 -> CH5+ + H
    kall(123) = 1.14e-10

    ! H2 + CH4+ -> CH5+ + H
    kall(124) = 3.30e-11

    ! H2+ + CH4 -> CH3+ + H2 + H
    kall(125) = 2.30e-09

    ! H2+ + OH -> H2O+ + H
    kall(126) = 7.60e-10

    ! H2 + OH+ -> H2O+ + H
    kall(127) = 1.01e-09

    ! H2+ + H2O -> H3O+ + H
    kall(128) = 3.40e-09

    ! H2 + H2O+ -> H3O+ + H
    kall(129) = 6.40e-10

    ! H2+ + CO -> HCO+ + H
    kall(130) = 2.16e-09

    ! H2 + CO+ -> HCO+ + H
    kall(131) = 1.80e-09

    ! H3+ + C -> CH+ + H2
    kall(132) = 2.00e-09

    ! H3+ + CH -> CH2+ + H2
    kall(133) = 1.20e-09

    ! H3+ + CH2 -> CH3+ + H2
    kall(134) = 1.70e-09

    ! H3+ + CH3 -> CH4+ + H2
    kall(135) = 2.10e-09

    ! H3+ + O -> OH+ + H2
    kall(136) = 8.00e-10

    ! H3+ + CH4 -> CH5+ + H2
    kall(137) = 2.40e-09

    ! H3+ + OH -> H2O+ + H2
    kall(138) = 1.30e-09

    ! H3+ + H2O -> H3O+ + H2
    kall(139) = 5.90e-09

    ! H3+ + CO -> HCO+ + H2
    kall(140) = 1.70e-09

    ! He+ + CH -> C+ + He + H
    kall(141) = 1.10e-09

    ! He+ + CH2 -> C+ + He + H2
    kall(142) = 7.50e-10

    ! He+ + CH2 -> CH+ + He + H
    kall(143) = 7.50e-10

    ! He+ + CH3 -> CH+ + He + H2
    kall(144) = 1.80e-09

    ! He+ + CH4 -> CH+ + He + H2 + H
    kall(145) = 2.40e-10

    ! He+ + CH4 -> CH2+ + He + H2
    kall(146) = 9.50e-10

    ! He+ + CH4 -> CH3 + He + H+
    kall(147) = 4.80e-10

    ! He+ + CH4 -> CH3+ + He + H
    kall(148) = 8.50e-11

    ! He+ + OH -> O+ + He + H
    kall(149) = 1.10e-09

    ! He+ + H2O -> OH + He + H+
    kall(150) = 2.04e-10

    ! He+ + H2O -> OH+ + He + H
    kall(151) = 2.86e-10

    ! He+ + CO -> O + C+ + He
    kall(152) = 1.60e-09

    ! He+ + O2 -> O+ + O + He
    kall(153) = 1.00e-09

    ! C + OH+ -> O + CH+
    kall(154) = 1.20e-09

    ! C+ + OH -> CO+ + H
    kall(155) = 7.70e-10

    ! C + CH5+ -> CH4 + CH+
    kall(156) = 1.20e-09

    ! C + H2O+ -> OH + CH+
    kall(157) = 1.10e-09

    ! C+ + H2O -> HCO+ + H
    kall(158) = 9.00e-10

    ! C + H3O+ -> HCO+ + H2
    kall(159) = 1.00e-11

    ! C + HCO+ -> CO + CH+
    kall(160) = 1.10e-09

    ! C+ + O2 -> CO+ + O
    kall(161) = 3.80e-10

    ! C+ + O2 -> CO + O+
    kall(162) = 6.20e-10

    ! C + O2+ -> CO+ + O
    kall(163) = 5.20e-11

    ! CH+ + O -> CO+ + H
    kall(164) = 3.50e-10

    ! CH + O+ -> CO+ + H
    kall(165) = 3.50e-10

    ! CH + OH+ -> O + CH2+
    kall(166) = 3.50e-10

    ! CH+ + OH -> CO+ + H2
    kall(167) = 7.50e-10

    ! CH + CH5+ -> CH4 + CH2+
    kall(168) = 6.90e-10

    ! CH+ + H2O -> H3O+ + C
    kall(169) = 5.80e-10

    ! CH + H2O+ -> OH + CH2+
    kall(170) = 3.40e-10

    ! CH+ + H2O -> HCO+ + H2
    kall(171) = 2.90e-09

    ! CH + H3O+ -> H2O + CH2+
    kall(172) = 6.80e-10

    ! CH + CO+ -> HCO+ + C
    kall(173) = 3.20e-10

    ! CH + HCO+ -> CO + CH2+
    kall(174) = 6.30e-10

    ! CH+ + O2 -> CO+ + OH
    kall(175) = 1.00e-11

    ! CH+ + O2 -> HCO+ + O
    kall(176) = 9.70e-10

    ! CH + O2+ -> HCO+ + O
    kall(177) = 3.10e-10

    ! CH2+ + O -> HCO+ + H
    kall(178) = 7.50e-10

    ! CH2 + OH+ -> O + CH3+
    kall(179) = 4.80e-10

    ! CH2 + CH5+ -> CH4 + CH3+
    kall(180) = 9.60e-10

    ! CH2 + H2O+ -> OH + CH3+
    kall(181) = 4.70e-10

    ! CH2 + H3O+ -> H2O + CH3+
    kall(182) = 9.40e-10

    ! CH2 + CO+ -> HCO+ + CH
    kall(183) = 4.30e-10

    ! CH2 + HCO+ -> CO + CH3+
    kall(184) = 8.60e-10

    ! CH2+ + O2 -> HCO+ + OH
    kall(185) = 9.10e-10

    ! CH3+ + O -> HCO+ + H2
    kall(186) = 4.00e-10

    ! O+ + CH4 -> OH + CH3+
    kall(187) = 1.10e-10

    ! O + CH4+ -> OH + CH3+
    kall(188) = 1.00e-09

    ! O+ + OH -> O2+ + H
    kall(189) = 3.60e-10

    ! O + OH+ -> O2+ + H
    kall(190) = 7.10e-10

    ! O + CH5+ -> H3O+ + CH2
    kall(191) = 2.20e-10

    ! O + H2O+ -> O2+ + H2
    kall(192) = 4.00e-11

    ! CH4+ + CH4 -> CH5+ + CH3
    kall(193) = 1.50e-09

    ! CH4 + OH+ -> CH5+ + O
    kall(194) = 1.95e-10

    ! CH4 + OH+ -> H3O+ + CH2
    kall(195) = 1.31e-09

    ! CH4+ + H2O -> H3O+ + CH3
    kall(196) = 2.60e-09

    ! CH4 + H2O+ -> H3O+ + CH3
    kall(197) = 1.40e-09

    ! CH4+ + CO -> HCO+ + CH3
    kall(198) = 1.40e-09

    ! CH4 + CO+ -> HCO+ + CH3
    kall(199) = 4.55e-10

    ! OH+ + OH -> H2O+ + O
    kall(200) = 7.00e-10

    ! OH + CH5+ -> H2O+ + CH4
    kall(201) = 7.00e-10

    ! OH+ + H2O -> H3O+ + O
    kall(202) = 1.30e-09

    ! OH + H2O+ -> H3O+ + O
    kall(203) = 6.90e-10

    ! OH+ + CO -> HCO+ + O
    kall(204) = 1.05e-09

    ! OH + CO+ -> HCO+ + O
    kall(205) = 3.10e-10

    ! OH + HCO+ -> CO + H2O+
    kall(206) = 6.20e-10

    ! CH5+ + H2O -> H3O+ + CH4
    kall(207) = 3.70e-09

    ! CH5+ + CO -> HCO+ + CH4
    kall(208) = 1.00e-09

    ! H2O+ + H2O -> H3O+ + OH
    kall(209) = 2.10e-09

    ! H2O+ + CO -> HCO+ + OH
    kall(210) = 5.00e-10

    ! H2O + CO+ -> HCO+ + OH
    kall(211) = 8.84e-10

    ! H2O + HCO+ -> CO + H3O+
    kall(212) = 2.50e-09

    ! H + H2+ -> H2 + H+
    kall(213) = 6.40e-10

    ! H+ + CH -> CH+ + H
    kall(214) = 1.90e-09

    ! H+ + CH2 -> CH2+ + H
    kall(215) = 1.40e-09

    ! H+ + CH3 -> CH3+ + H
    kall(216) = 3.40e-09

    ! H+ + CH4 -> CH4+ + H
    kall(217) = 1.50e-09

    ! H+ + OH -> OH+ + H
    kall(218) = 2.10e-09

    ! H + CO+ -> CO + H+
    kall(219) = 7.50e-10

    ! H2+ + CH -> CH+ + H2
    kall(220) = 7.10e-10

    ! H2+ + CH2 -> CH2+ + H2
    kall(221) = 1.00e-09

    ! H2+ + CH4 -> CH4+ + H2
    kall(222) = 1.40e-09

    ! H2+ + OH -> OH+ + H2
    kall(223) = 7.60e-10

    ! H2+ + H2O -> H2O+ + H2
    kall(224) = 3.90e-09

    ! H2+ + CO -> CO+ + H2
    kall(225) = 6.40e-10

    ! H2+ + O2 -> O2+ + H2
    kall(226) = 8.00e-10

    ! He+ + CH -> CH+ + He
    kall(227) = 5.00e-10

    ! He+ + CH4 -> CH4+ + He
    kall(228) = 5.10e-11

    ! He+ + H2O -> H2O+ + He
    kall(229) = 6.05e-11

    ! He+ + O2 -> O2+ + He
    kall(230) = 3.30e-11

    ! C+ + CH -> CH+ + C
    kall(231) = 3.80e-10

    ! C+ + CH2 -> CH2+ + C
    kall(232) = 5.20e-10

    ! C + CO+ -> CO + C+
    kall(233) = 1.10e-10

    ! C + O2+ -> O2 + C+
    kall(234) = 5.20e-11

    ! CH + O+ -> O + CH+
    kall(235) = 3.50e-10

    ! CH + OH+ -> OH + CH+
    kall(236) = 3.50e-10

    ! CH + H2O+ -> H2O + CH+
    kall(237) = 3.40e-10

    ! CH + CO+ -> CO + CH+
    kall(238) = 3.20e-10

    ! CH + O2+ -> O2 + CH+
    kall(239) = 3.10e-10

    ! CH2 + O+ -> O + CH2+
    kall(240) = 9.70e-10

    ! CH2 + OH+ -> OH + CH2+
    kall(241) = 4.80e-10

    ! CH2 + H2O+ -> H2O + CH2+
    kall(242) = 4.70e-10

    ! CH2 + CO+ -> CO + CH2+
    kall(243) = 4.30e-10

    ! CH2 + O2+ -> O2 + CH2+
    kall(244) = 4.30e-10

    ! O+ + CH4 -> CH4+ + O
    kall(245) = 8.90e-10

    ! O+ + OH -> OH+ + O
    kall(246) = 3.60e-10

    ! O+ + H2O -> H2O+ + O
    kall(247) = 3.20e-09

    ! O + CO+ -> CO + O+
    kall(248) = 1.40e-10

    ! O+ + O2 -> O2+ + O
    kall(249) = 1.90e-11

    ! CH4 + CO+ -> CO + CH4+
    kall(250) = 7.93e-10

    ! CH4+ + O2 -> O2+ + CH4
    kall(251) = 4.00e-10

    ! OH+ + H2O -> H2O+ + OH
    kall(252) = 1.59e-09

    ! OH + CO+ -> CO + OH+
    kall(253) = 3.10e-10

    ! OH+ + O2 -> O2+ + OH
    kall(254) = 5.90e-10

    ! H2O + CO+ -> CO + H2O+
    kall(255) = 1.72e-09

    ! H2O+ + O2 -> O2+ + H2O
    kall(256) = 4.60e-10

    ! CO+ + O2 -> O2+ + CO
    kall(257) = 1.20e-10

    ! H2 -> H + H
    ! kall(258) = photochemistry, see rates_photo.f90

    ! CO -> C + O
    ! kall(259) = photochemistry, see rates_photo.f90

    ! H2+ -> H+ + H
    ! kall(260) = photochemistry, see rates_photo.f90

    ! H3+ -> H2+ + H
    ! kall(261) = photochemistry, see rates_photo.f90

    ! H3+ -> H2 + H+
    ! kall(262) = photochemistry, see rates_photo.f90

    ! C -> C+ + E
    ! kall(263) = photochemistry, see rates_photo.f90

    ! CH -> CH+ + E
    ! kall(264) = photochemistry, see rates_photo.f90

    ! CH -> C + H
    ! kall(265) = photochemistry, see rates_photo.f90

    ! CH+ -> C+ + H
    ! kall(266) = photochemistry, see rates_photo.f90

    ! CH2 -> CH2+ + E
    ! kall(267) = photochemistry, see rates_photo.f90

    ! CH2 -> CH + H
    ! kall(268) = photochemistry, see rates_photo.f90

    ! CH2+ -> CH+ + H
    ! kall(269) = photochemistry, see rates_photo.f90

    ! CH3 -> CH3+ + E
    ! kall(270) = photochemistry, see rates_photo.f90

    ! CH3 -> CH2 + H
    ! kall(271) = photochemistry, see rates_photo.f90

    ! CH3 -> CH + H2
    ! kall(272) = photochemistry, see rates_photo.f90

    ! CH4 -> CH3 + H
    ! kall(273) = photochemistry, see rates_photo.f90

    ! CH4 -> CH2 + H2
    ! kall(274) = photochemistry, see rates_photo.f90

    ! CH4 -> CH + H2 + H
    ! kall(275) = photochemistry, see rates_photo.f90

    ! OH -> OH+ + E
    ! kall(276) = photochemistry, see rates_photo.f90

    ! OH -> O + H
    ! kall(277) = photochemistry, see rates_photo.f90

    ! OH+ -> O + H+
    ! kall(278) = photochemistry, see rates_photo.f90

    ! H2O -> OH + H
    ! kall(279) = photochemistry, see rates_photo.f90

    ! H2O -> H2O+ + E
    ! kall(280) = photochemistry, see rates_photo.f90

    ! O2 -> O + O
    ! kall(281) = photochemistry, see rates_photo.f90

    ! O2 -> O2+ + E
    ! kall(282) = photochemistry, see rates_photo.f90

    ! H2 -> H+ + H + E
    kall(283) = 0.02e+00 * variable_crflux

    ! H2 -> H + H
    kall(284) = 0.10e+00 * variable_crflux

    ! H2 -> H2+ + E
    kall(285) = 0.88e+00 * variable_crflux

    ! H -> H+ + E
    kall(286) = 0.46e+00 * variable_crflux

    ! He -> He+ + E
    kall(287) = 0.50e+00 * variable_crflux

    ! load evaluate once rate into common variable for later use
    krate_evaluate_once(95:287) = kall(95:287)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_INIT_EVALUATE_ONCE

  end subroutine init_evaluate_once

end module prizmo_rates_evaluate_once
