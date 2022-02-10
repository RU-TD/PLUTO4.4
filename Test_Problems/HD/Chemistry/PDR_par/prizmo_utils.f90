module prizmo_utils
contains

  ! ****************************
  ! get electron number density from the other charged species
  function get_electrons(n) result(ne)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(nmols)
    real*8::ne

    !!BEGIN_GET_ELECTRONS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ne = 0d0 &
        + n(idx_CH2j) &
        + n(idx_CHj) &
        + n(idx_CH3j) &
        + n(idx_Hej) &
        + n(idx_Hj) &
        + n(idx_Cj) &
        + n(idx_Oj) &
        + n(idx_H2j) &
        + n(idx_COj) &
        + n(idx_H3j) &
        + n(idx_CH4j) &
        + n(idx_OHj) &
        + n(idx_CH5j) &
        + n(idx_H2Oj) &
        + n(idx_H3Oj) &
        + n(idx_HCOj) &
        + n(idx_O2j)
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_GET_ELECTRONS

  end function get_electrons

  ! **********************
  function get_species_name(idx) result(name)
    use prizmo_commons
    implicit none
    integer,intent(in):: idx
    character(len=max_character_len)::names(nmols), name

    names(:) = get_species_names()

    name = trim(names(idx))

  end function get_species_name

  ! **********************
  function get_species_names() result(names)
    use prizmo_commons
    implicit none
    character(len=max_character_len)::names(nmols)


    !!BEGIN_SPECIES_NAMES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    names(1) = "H"
    names(2) = "CH"
    names(3) = "C"
    names(4) = "H2"
    names(5) = "CH3"
    names(6) = "CH2"
    names(7) = "CH4"
    names(8) = "OH"
    names(9) = "O"
    names(10) = "H2O"
    names(11) = "CO"
    names(12) = "O2"
    names(13) = "CH2+"
    names(14) = "CH+"
    names(15) = "CH3+"
    names(16) = "He+"
    names(17) = "He"
    names(18) = "H+"
    names(19) = "C+"
    names(20) = "O+"
    names(21) = "H2+"
    names(22) = "CO+"
    names(23) = "E"
    names(24) = "H3+"
    names(25) = "CH4+"
    names(26) = "OH+"
    names(27) = "CH5+"
    names(28) = "H2O+"
    names(29) = "H3O+"
    names(30) = "HCO+"
    names(31) = "O2+"

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_SPECIES_NAMES

  end function get_species_names

  ! **********************
  ! get an array of string with reaction verbatim
  function get_reaction_names() result(rnames)
    use prizmo_commons
    implicit none
    character(len=max_character_len)::rnames(nrea)

    !!BEGIN_REACTION_NAMES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    rnames(1) = "H + CH -> C + H2"
    rnames(2) = "H + CH3 -> CH2 + H2"
    rnames(3) = "H + CH4 -> CH3 + H2"
    rnames(4) = "H + OH -> O + H2"
    rnames(5) = "H + H2O -> OH + H2"
    rnames(6) = "H + CO -> OH + C"
    rnames(7) = "H + O2 -> OH + O"
    rnames(8) = "H2 + C -> CH + H"
    rnames(9) = "H2 + CH -> CH2 + H"
    rnames(10) = "H2 + CH2 -> CH3 + H"
    rnames(11) = "H2 + CH3 -> CH4 + H"
    rnames(12) = "H2 + O -> OH + H"
    rnames(13) = "H2 + OH -> H2O + H"
    rnames(14) = "H2 + O2 -> OH + OH"
    rnames(15) = "C + CH2 -> CH + CH"
    rnames(16) = "C + OH -> O + CH"
    rnames(17) = "C + OH -> CO + H"
    rnames(18) = "CH + O -> OH + C"
    rnames(19) = "CH + CH4 -> CH3 + CH2"
    rnames(20) = "CH2 + CH2 -> CH3 + CH"
    rnames(21) = "CH2 + O -> OH + CH"
    rnames(22) = "CH2 + CH4 -> CH3 + CH3"
    rnames(23) = "CH2 + OH -> O + CH3"
    rnames(24) = "CH2 + OH -> H2O + CH"
    rnames(25) = "CH2 + O2 -> CO + H2O"
    rnames(26) = "CH3 + CH3 -> CH4 + CH2"
    rnames(27) = "CH3 + OH -> CH4 + O"
    rnames(28) = "CH3 + OH -> H2O + CH2"
    rnames(29) = "CH3 + H2O -> OH + CH4"
    rnames(30) = "O + CH4 -> OH + CH3"
    rnames(31) = "O + OH -> O2 + H"
    rnames(32) = "O + H2O -> OH + OH"
    rnames(33) = "CH4 + OH -> H2O + CH3"
    rnames(34) = "OH + OH -> H2O + O"
    rnames(35) = "H + CH2+ -> CH+ + H2"
    rnames(36) = "H + CH3+ -> CH2+ + H2"
    rnames(37) = "H2 + He+ -> He + H+ + H"
    rnames(38) = "H2 + C+ -> CH+ + H"
    rnames(39) = "He+ + CO -> O+ + C + He"
    rnames(40) = "H+ + H2 -> H2+ + H"
    rnames(41) = "H + He+ -> He + H+"
    rnames(42) = "H + O+ -> O + H+"
    rnames(43) = "H+ + O -> O+ + H"
    rnames(44) = "He+ + C -> C+ + He"
    rnames(45) = "O+ + CO -> CO+ + O"
    rnames(46) = "H2+ + E -> H + H"
    rnames(47) = "H3+ + E -> H + H + H"
    rnames(48) = "H3+ + E -> H2 + H"
    rnames(49) = "CH+ + E -> C + H"
    rnames(50) = "CH2+ + E -> CH + H"
    rnames(51) = "CH2+ + E -> C + H + H"
    rnames(52) = "CH2+ + E -> C + H2"
    rnames(53) = "CH3+ + E -> CH2 + H"
    rnames(54) = "CH3+ + E -> CH + H2"
    rnames(55) = "CH3+ + E -> CH + H + H"
    rnames(56) = "CH4+ + E -> CH3 + H"
    rnames(57) = "CH4+ + E -> CH2 + H + H"
    rnames(58) = "OH+ + E -> O + H"
    rnames(59) = "CH5+ + E -> CH3 + H2"
    rnames(60) = "CH5+ + E -> CH4 + H"
    rnames(61) = "H2O+ + E -> O + H + H"
    rnames(62) = "H2O+ + E -> O + H2"
    rnames(63) = "H2O+ + E -> OH + H"
    rnames(64) = "H3O+ + E -> O + H2 + H"
    rnames(65) = "H3O+ + E -> OH + H + H"
    rnames(66) = "H3O+ + E -> OH + H2"
    rnames(67) = "H3O+ + E -> H2O + H"
    rnames(68) = "CO+ + E -> O + C"
    rnames(69) = "HCO+ + E -> CO + H"
    rnames(70) = "O2+ + E -> O + O"
    rnames(71) = "H+ + E -> H"
    rnames(72) = "He+ + E -> He"
    rnames(73) = "C+ + E -> C"
    rnames(74) = "CH3+ + E -> CH3"
    rnames(75) = "O+ + E -> O"
    rnames(76) = "H+ + H -> H2+"
    rnames(77) = "H + O -> OH"
    rnames(78) = "H + OH -> H2O"
    rnames(79) = "H2 + C+ -> CH2+"
    rnames(80) = "H2 + CH -> CH3"
    rnames(81) = "H2 + CH3+ -> CH5+"
    rnames(82) = "O + O -> O2"
    rnames(83) = "H + H2 -> H + H + H"
    rnames(84) = "H + CH -> C + H + H"
    rnames(85) = "H + OH -> O + H + H"
    rnames(86) = "H + H2O -> OH + H + H"
    rnames(87) = "H + O2 -> O + O + H"
    rnames(88) = "H2 + E -> H + H + E"
    rnames(89) = "H2 + H2 -> H2 + H + H"
    rnames(90) = "H2 + CH -> C + H2 + H"
    rnames(91) = "H2 + OH -> O + H2 + H"
    rnames(92) = "H2 + H2O -> OH + H2 + H"
    rnames(93) = "H2 + O2 -> O + O + H2"
    rnames(94) = "CH + O -> HCO+ + E"
    rnames(95) = "H + CH2 -> CH + H2"
    rnames(96) = "C + O2 -> CO + O"
    rnames(97) = "CH + O -> CO + H"
    rnames(98) = "CH + O2 -> CO + OH"
    rnames(99) = "CH2 + O -> CO + H + H"
    rnames(100) = "CH2 + O -> CO + H2"
    rnames(101) = "H + CH+ -> C+ + H2"
    rnames(102) = "H2 + CH+ -> CH2+ + H"
    rnames(103) = "C+ + OH -> CO + H+"
    rnames(104) = "H+ + H2O -> H2O+ + H"
    rnames(105) = "H+ + O2 -> O2+ + H"
    rnames(106) = "H2 + He+ -> He + H2+"
    rnames(107) = "H + C -> CH"
    rnames(108) = "H + C+ -> CH+"
    rnames(109) = "H2 + C -> CH2"
    rnames(110) = "C + O -> CO"
    rnames(111) = "C+ + O -> CO+"
    rnames(112) = "H+ + CH2 -> CH+ + H2"
    rnames(113) = "H+ + CH4 -> CH3+ + H2"
    rnames(114) = "H + CH4+ -> CH3+ + H2"
    rnames(115) = "H + CH5+ -> CH4+ + H2"
    rnames(116) = "H2+ + H2 -> H3+ + H"
    rnames(117) = "H2+ + C -> CH+ + H"
    rnames(118) = "H2+ + CH -> CH2+ + H"
    rnames(119) = "H2+ + CH2 -> CH3+ + H"
    rnames(120) = "H2 + CH2+ -> CH3+ + H"
    rnames(121) = "H2+ + O -> OH+ + H"
    rnames(122) = "H2 + O+ -> OH+ + H"
    rnames(123) = "H2+ + CH4 -> CH5+ + H"
    rnames(124) = "H2 + CH4+ -> CH5+ + H"
    rnames(125) = "H2+ + CH4 -> CH3+ + H2 + H"
    rnames(126) = "H2+ + OH -> H2O+ + H"
    rnames(127) = "H2 + OH+ -> H2O+ + H"
    rnames(128) = "H2+ + H2O -> H3O+ + H"
    rnames(129) = "H2 + H2O+ -> H3O+ + H"
    rnames(130) = "H2+ + CO -> HCO+ + H"
    rnames(131) = "H2 + CO+ -> HCO+ + H"
    rnames(132) = "H3+ + C -> CH+ + H2"
    rnames(133) = "H3+ + CH -> CH2+ + H2"
    rnames(134) = "H3+ + CH2 -> CH3+ + H2"
    rnames(135) = "H3+ + CH3 -> CH4+ + H2"
    rnames(136) = "H3+ + O -> OH+ + H2"
    rnames(137) = "H3+ + CH4 -> CH5+ + H2"
    rnames(138) = "H3+ + OH -> H2O+ + H2"
    rnames(139) = "H3+ + H2O -> H3O+ + H2"
    rnames(140) = "H3+ + CO -> HCO+ + H2"
    rnames(141) = "He+ + CH -> C+ + He + H"
    rnames(142) = "He+ + CH2 -> C+ + He + H2"
    rnames(143) = "He+ + CH2 -> CH+ + He + H"
    rnames(144) = "He+ + CH3 -> CH+ + He + H2"
    rnames(145) = "He+ + CH4 -> CH+ + He + H2 + H"
    rnames(146) = "He+ + CH4 -> CH2+ + He + H2"
    rnames(147) = "He+ + CH4 -> CH3 + He + H+"
    rnames(148) = "He+ + CH4 -> CH3+ + He + H"
    rnames(149) = "He+ + OH -> O+ + He + H"
    rnames(150) = "He+ + H2O -> OH + He + H+"
    rnames(151) = "He+ + H2O -> OH+ + He + H"
    rnames(152) = "He+ + CO -> O + C+ + He"
    rnames(153) = "He+ + O2 -> O+ + O + He"
    rnames(154) = "C + OH+ -> O + CH+"
    rnames(155) = "C+ + OH -> CO+ + H"
    rnames(156) = "C + CH5+ -> CH4 + CH+"
    rnames(157) = "C + H2O+ -> OH + CH+"
    rnames(158) = "C+ + H2O -> HCO+ + H"
    rnames(159) = "C + H3O+ -> HCO+ + H2"
    rnames(160) = "C + HCO+ -> CO + CH+"
    rnames(161) = "C+ + O2 -> CO+ + O"
    rnames(162) = "C+ + O2 -> CO + O+"
    rnames(163) = "C + O2+ -> CO+ + O"
    rnames(164) = "CH+ + O -> CO+ + H"
    rnames(165) = "CH + O+ -> CO+ + H"
    rnames(166) = "CH + OH+ -> O + CH2+"
    rnames(167) = "CH+ + OH -> CO+ + H2"
    rnames(168) = "CH + CH5+ -> CH4 + CH2+"
    rnames(169) = "CH+ + H2O -> H3O+ + C"
    rnames(170) = "CH + H2O+ -> OH + CH2+"
    rnames(171) = "CH+ + H2O -> HCO+ + H2"
    rnames(172) = "CH + H3O+ -> H2O + CH2+"
    rnames(173) = "CH + CO+ -> HCO+ + C"
    rnames(174) = "CH + HCO+ -> CO + CH2+"
    rnames(175) = "CH+ + O2 -> CO+ + OH"
    rnames(176) = "CH+ + O2 -> HCO+ + O"
    rnames(177) = "CH + O2+ -> HCO+ + O"
    rnames(178) = "CH2+ + O -> HCO+ + H"
    rnames(179) = "CH2 + OH+ -> O + CH3+"
    rnames(180) = "CH2 + CH5+ -> CH4 + CH3+"
    rnames(181) = "CH2 + H2O+ -> OH + CH3+"
    rnames(182) = "CH2 + H3O+ -> H2O + CH3+"
    rnames(183) = "CH2 + CO+ -> HCO+ + CH"
    rnames(184) = "CH2 + HCO+ -> CO + CH3+"
    rnames(185) = "CH2+ + O2 -> HCO+ + OH"
    rnames(186) = "CH3+ + O -> HCO+ + H2"
    rnames(187) = "O+ + CH4 -> OH + CH3+"
    rnames(188) = "O + CH4+ -> OH + CH3+"
    rnames(189) = "O+ + OH -> O2+ + H"
    rnames(190) = "O + OH+ -> O2+ + H"
    rnames(191) = "O + CH5+ -> H3O+ + CH2"
    rnames(192) = "O + H2O+ -> O2+ + H2"
    rnames(193) = "CH4+ + CH4 -> CH5+ + CH3"
    rnames(194) = "CH4 + OH+ -> CH5+ + O"
    rnames(195) = "CH4 + OH+ -> H3O+ + CH2"
    rnames(196) = "CH4+ + H2O -> H3O+ + CH3"
    rnames(197) = "CH4 + H2O+ -> H3O+ + CH3"
    rnames(198) = "CH4+ + CO -> HCO+ + CH3"
    rnames(199) = "CH4 + CO+ -> HCO+ + CH3"
    rnames(200) = "OH+ + OH -> H2O+ + O"
    rnames(201) = "OH + CH5+ -> H2O+ + CH4"
    rnames(202) = "OH+ + H2O -> H3O+ + O"
    rnames(203) = "OH + H2O+ -> H3O+ + O"
    rnames(204) = "OH+ + CO -> HCO+ + O"
    rnames(205) = "OH + CO+ -> HCO+ + O"
    rnames(206) = "OH + HCO+ -> CO + H2O+"
    rnames(207) = "CH5+ + H2O -> H3O+ + CH4"
    rnames(208) = "CH5+ + CO -> HCO+ + CH4"
    rnames(209) = "H2O+ + H2O -> H3O+ + OH"
    rnames(210) = "H2O+ + CO -> HCO+ + OH"
    rnames(211) = "H2O + CO+ -> HCO+ + OH"
    rnames(212) = "H2O + HCO+ -> CO + H3O+"
    rnames(213) = "H + H2+ -> H2 + H+"
    rnames(214) = "H+ + CH -> CH+ + H"
    rnames(215) = "H+ + CH2 -> CH2+ + H"
    rnames(216) = "H+ + CH3 -> CH3+ + H"
    rnames(217) = "H+ + CH4 -> CH4+ + H"
    rnames(218) = "H+ + OH -> OH+ + H"
    rnames(219) = "H + CO+ -> CO + H+"
    rnames(220) = "H2+ + CH -> CH+ + H2"
    rnames(221) = "H2+ + CH2 -> CH2+ + H2"
    rnames(222) = "H2+ + CH4 -> CH4+ + H2"
    rnames(223) = "H2+ + OH -> OH+ + H2"
    rnames(224) = "H2+ + H2O -> H2O+ + H2"
    rnames(225) = "H2+ + CO -> CO+ + H2"
    rnames(226) = "H2+ + O2 -> O2+ + H2"
    rnames(227) = "He+ + CH -> CH+ + He"
    rnames(228) = "He+ + CH4 -> CH4+ + He"
    rnames(229) = "He+ + H2O -> H2O+ + He"
    rnames(230) = "He+ + O2 -> O2+ + He"
    rnames(231) = "C+ + CH -> CH+ + C"
    rnames(232) = "C+ + CH2 -> CH2+ + C"
    rnames(233) = "C + CO+ -> CO + C+"
    rnames(234) = "C + O2+ -> O2 + C+"
    rnames(235) = "CH + O+ -> O + CH+"
    rnames(236) = "CH + OH+ -> OH + CH+"
    rnames(237) = "CH + H2O+ -> H2O + CH+"
    rnames(238) = "CH + CO+ -> CO + CH+"
    rnames(239) = "CH + O2+ -> O2 + CH+"
    rnames(240) = "CH2 + O+ -> O + CH2+"
    rnames(241) = "CH2 + OH+ -> OH + CH2+"
    rnames(242) = "CH2 + H2O+ -> H2O + CH2+"
    rnames(243) = "CH2 + CO+ -> CO + CH2+"
    rnames(244) = "CH2 + O2+ -> O2 + CH2+"
    rnames(245) = "O+ + CH4 -> CH4+ + O"
    rnames(246) = "O+ + OH -> OH+ + O"
    rnames(247) = "O+ + H2O -> H2O+ + O"
    rnames(248) = "O + CO+ -> CO + O+"
    rnames(249) = "O+ + O2 -> O2+ + O"
    rnames(250) = "CH4 + CO+ -> CO + CH4+"
    rnames(251) = "CH4+ + O2 -> O2+ + CH4"
    rnames(252) = "OH+ + H2O -> H2O+ + OH"
    rnames(253) = "OH + CO+ -> CO + OH+"
    rnames(254) = "OH+ + O2 -> O2+ + OH"
    rnames(255) = "H2O + CO+ -> CO + H2O+"
    rnames(256) = "H2O+ + O2 -> O2+ + H2O"
    rnames(257) = "CO+ + O2 -> O2+ + CO"
    rnames(258) = "H2 -> H + H"
    rnames(259) = "CO -> C + O"
    rnames(260) = "H2+ -> H+ + H"
    rnames(261) = "H3+ -> H2+ + H"
    rnames(262) = "H3+ -> H2 + H+"
    rnames(263) = "C -> C+ + E"
    rnames(264) = "CH -> CH+ + E"
    rnames(265) = "CH -> C + H"
    rnames(266) = "CH+ -> C+ + H"
    rnames(267) = "CH2 -> CH2+ + E"
    rnames(268) = "CH2 -> CH + H"
    rnames(269) = "CH2+ -> CH+ + H"
    rnames(270) = "CH3 -> CH3+ + E"
    rnames(271) = "CH3 -> CH2 + H"
    rnames(272) = "CH3 -> CH + H2"
    rnames(273) = "CH4 -> CH3 + H"
    rnames(274) = "CH4 -> CH2 + H2"
    rnames(275) = "CH4 -> CH + H2 + H"
    rnames(276) = "OH -> OH+ + E"
    rnames(277) = "OH -> O + H"
    rnames(278) = "OH+ -> O + H+"
    rnames(279) = "H2O -> OH + H"
    rnames(280) = "H2O -> H2O+ + E"
    rnames(281) = "O2 -> O + O"
    rnames(282) = "O2 -> O2+ + E"
    rnames(283) = "H2 -> H+ + H + E"
    rnames(284) = "H2 -> H + H"
    rnames(285) = "H2 -> H2+ + E"
    rnames(286) = "H -> H+ + E"
    rnames(287) = "He -> He+ + E"
    rnames(288) = "CH3+ + E -> C + H2 + H"
    rnames(289) = "H + H -> H2"

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_REACTION_NAMES

  end function get_reaction_names

  ! *********************
  ! get species masses in g
  ! WARNING: this is an alias. use the common variable where possible
  function get_species_mass() result(masses)
    use prizmo_commons
    implicit none
    real*8::masses(nmols)

    masses(:) = mass(:)

  end function get_species_mass

  ! *********************************
  ! get total gas density, g/cm3
  function get_rho(n) result(rho)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(nmols)
    real*8::rho

    ! rho = sum(mass(:) * n(:))
    ! reduced version
    !!BEGIN_GET_RHO
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    rho =  &
        + mass(idx_H) * n(idx_H) &
        + mass(idx_CH) * n(idx_CH) &
        + mass(idx_C) * n(idx_C) &
        + mass(idx_H2) * n(idx_H2) &
        + mass(idx_CH3) * n(idx_CH3) &
        + mass(idx_CH2) * n(idx_CH2) &
        + mass(idx_CH4) * n(idx_CH4) &
        + mass(idx_OH) * n(idx_OH) &
        + mass(idx_O) * n(idx_O) &
        + mass(idx_H2O) * n(idx_H2O) &
        + mass(idx_CO) * n(idx_CO) &
        + mass(idx_O2) * n(idx_O2) &
        + mass(idx_CH2j) * n(idx_CH2j) &
        + mass(idx_CHj) * n(idx_CHj) &
        + mass(idx_CH3j) * n(idx_CH3j) &
        + mass(idx_Hej) * n(idx_Hej) &
        + mass(idx_He) * n(idx_He) &
        + mass(idx_Hj) * n(idx_Hj) &
        + mass(idx_Cj) * n(idx_Cj) &
        + mass(idx_Oj) * n(idx_Oj) &
        + mass(idx_H2j) * n(idx_H2j) &
        + mass(idx_COj) * n(idx_COj) &
        + mass(idx_E) * n(idx_E) &
        + mass(idx_H3j) * n(idx_H3j) &
        + mass(idx_CH4j) * n(idx_CH4j) &
        + mass(idx_OHj) * n(idx_OHj) &
        + mass(idx_CH5j) * n(idx_CH5j) &
        + mass(idx_H2Oj) * n(idx_H2Oj) &
        + mass(idx_H3Oj) * n(idx_H3Oj) &
        + mass(idx_HCOj) * n(idx_HCOj) &
        + mass(idx_O2j) * n(idx_O2j)
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_GET_RHO

  end function get_rho

  ! *********************
  ! get total density, cm-3
  function get_ntot(n) result(ntot)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(:)
    real*8::ntot

    !!BEGIN_GET_NTOT
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ntot = 0d0 &
        + n(idx_H) &
        + n(idx_CH) &
        + n(idx_C) &
        + n(idx_H2) &
        + n(idx_CH3) &
        + n(idx_CH2) &
        + n(idx_CH4) &
        + n(idx_OH) &
        + n(idx_O) &
        + n(idx_H2O) &
        + n(idx_CO) &
        + n(idx_O2) &
        + n(idx_CH2j) &
        + n(idx_CHj) &
        + n(idx_CH3j) &
        + n(idx_Hej) &
        + n(idx_He) &
        + n(idx_Hj) &
        + n(idx_Cj) &
        + n(idx_Oj) &
        + n(idx_H2j) &
        + n(idx_COj) &
        + n(idx_E) &
        + n(idx_H3j) &
        + n(idx_CH4j) &
        + n(idx_OHj) &
        + n(idx_CH5j) &
        + n(idx_H2Oj) &
        + n(idx_H3Oj) &
        + n(idx_HCOj) &
        + n(idx_O2j)
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_GET_NTOT
    !ntot = n(idx_H2) + n(idx_H) + n(idx_Hj)  + n(idx_H2j) + n(idx_He) + n(idx_Hej) + n(idx_E)

  end function get_ntot

  !!BEGIN_GET_X_NUCLEI
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! *********************
    ! get number of He nuclei, cm-3
    function get_Henuclei(n) result(nX)
        use prizmo_commons
        implicit none
        real*8,intent(in)::n(:)
        real*8::nX

        nX = 0d0 +n(idx_Hej) &
            + n(idx_He)

    end function get_Henuclei

    ! *********************
    ! get number of C nuclei, cm-3
    function get_Cnuclei(n) result(nX)
        use prizmo_commons
        implicit none
        real*8,intent(in)::n(:)
        real*8::nX

        nX = 0d0 +n(idx_CH) &
            + n(idx_C) &
            + n(idx_CH3) &
            + n(idx_CH2) &
            + n(idx_CH4) &
            + n(idx_CO) &
            + n(idx_CH2j) &
            + n(idx_CHj) &
            + n(idx_CH3j) &
            + n(idx_Cj) &
            + n(idx_COj) &
            + n(idx_CH4j) &
            + n(idx_CH5j) &
            + n(idx_HCOj)

    end function get_Cnuclei

    ! *********************
    ! get number of O nuclei, cm-3
    function get_Onuclei(n) result(nX)
        use prizmo_commons
        implicit none
        real*8,intent(in)::n(:)
        real*8::nX

        nX = 0d0 +n(idx_OH) &
            + n(idx_O) &
            + n(idx_H2O) &
            + n(idx_CO) &
            + 2 * n(idx_O2) &
            + n(idx_Oj) &
            + n(idx_COj) &
            + n(idx_OHj) &
            + n(idx_H2Oj) &
            + n(idx_H3Oj) &
            + n(idx_HCOj) &
            + 2 * n(idx_O2j)

    end function get_Onuclei

    ! *********************
    ! get number of H nuclei, cm-3
    function get_Hnuclei(n) result(nX)
        use prizmo_commons
        implicit none
        real*8,intent(in)::n(:)
        real*8::nX

        nX = 0d0 +n(idx_H) &
            + n(idx_CH) &
            + 2 * n(idx_H2) &
            + 3 * n(idx_CH3) &
            + 2 * n(idx_CH2) &
            + 4 * n(idx_CH4) &
            + n(idx_OH) &
            + 2 * n(idx_H2O) &
            + 2 * n(idx_CH2j) &
            + n(idx_CHj) &
            + 3 * n(idx_CH3j) &
            + n(idx_Hj) &
            + 2 * n(idx_H2j) &
            + 3 * n(idx_H3j) &
            + 4 * n(idx_CH4j) &
            + n(idx_OHj) &
            + 5 * n(idx_CH5j) &
            + 2 * n(idx_H2Oj) &
            + 3 * n(idx_H3Oj) &
            + n(idx_HCOj)

    end function get_Hnuclei

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_GET_X_NUCLEI

  ! **********************
  ! get adiabatic index
  function get_gamma_ad(n, Tgas) result(gamma_ad)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(:), Tgas
    real*8::gamma_ad

    gamma_ad = 7d0 / 5d0

  end function get_gamma_ad

  ! *******************
  ! mean molecular weight, no dimensions
  function get_mu(n) result(mu)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(:)
    real*8::mu

    mu = sum(n(:) * mass(:)) / sum(n(:)) / proton_mass

  end function get_mu

  ! *******************
  ! dimensional mean molecular weight, g
  function get_mu_g(n) result(mu)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(:)
    real*8::mu

    mu = sum(n(:) * mass(:)) / sum(n(:))

end function get_mu_g

end module prizmo_utils
