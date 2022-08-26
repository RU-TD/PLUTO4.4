module prizmo_ode
  use prizmo_commons
  use prizmo_flux
  use prizmo_heating
  use prizmo_cooling
  use prizmo_tdust
contains

  ! ***************************
  subroutine fex(neq, tt, yin, dy)
    integer::neq, i
    real*8::yin(neq), tt, dy(neq), y(neq)
    real*8::flux(nreactions), Tdust, jflux(nphoto), ntot, log_tgas, log_ngas
    real*8:: heat, cool, Tgas

    !y = max(yin, 0d0)
    y = yin
    Tgas = max(y(idx_tgas), 3d0)  ! FIXME

    jflux = jflux_common
    ntot = sum(max(y(1:nspecies), 1d-40))

    log_tgas = log10(Tgas)
    log_ngas = log10(ntot)
    Tdust = get_tdust(log_Tgas, log_ngas)

    flux(:) = get_flux(y(1:nspecies), Tgas, Tdust)

    if(minval(kall) < 0d0) then
      print *, "ERROR: negative kall!!!"
      do i=1,nreactions
        print *, i, kall(i)
      end do
      do i=1,nspecies
        print *, i, y(i)
      end do
      stop
    end if

    if(solve_thermo) then
      heat = heating(y(1:nspecies), Tgas, Tdust, jflux)
      cool = cooling(y(1:nspecies), Tgas, Tdust, jflux, flux)

      dy(idx_Tgas) = (gamma_ad - 1d0) * (heat - cool) / kboltzmann / ntot
    end if

    if(solve_chemistry) then
      !! PREPROCESS_ODE

      dy(idx_C) =  - flux(9) - flux(16) - flux(17) - flux(18) - flux(19) - flux(56) - flux(86) - flux(90) - flux(94) - flux(115) - flux(130) - flux(152) - flux(154) - flux(155) - flux(157) - flux(158) - flux(161) - flux(231) - flux(232) - flux(271) - flux(290) + flux(1) + flux(7) + flux(20) + flux(47) + flux(61) + flux(63) + flux(64) + flux(68) + flux(81) + flux(98) + flux(104) + flux(167) + flux(171) + flux(229) + flux(230) + flux(260) + flux(295) + flux(298)
      dy(idx_CH) =  - flux(1) - flux(10) - flux(20) - flux(21) - flux(22) - flux(23) - flux(92) - flux(98) - flux(104) - flux(108) - flux(116) - flux(131) - flux(139) - flux(163) - flux(164) - flux(166) - flux(168) - flux(170) - flux(171) - flux(172) - flux(175) - flux(212) - flux(218) - flux(225) - flux(229) - flux(233) - flux(234) - flux(235) - flux(236) - flux(237) - flux(297) - flux(298) + flux(2) + flux(9) + flux(16) + flux(16) + flux(17) + flux(24) + flux(25) + flux(30) + flux(62) + flux(66) + flux(67) + flux(86) + flux(181) + flux(300)
      dy(idx_CH2) =  - flux(2) - flux(11) - flux(16) - flux(24) - flux(24) - flux(25) - flux(26) - flux(27) - flux(28) - flux(29) - flux(30) - flux(31) - flux(110) - flux(117) - flux(132) - flux(140) - flux(141) - flux(177) - flux(178) - flux(179) - flux(180) - flux(181) - flux(182) - flux(213) - flux(219) - flux(230) - flux(238) - flux(239) - flux(240) - flux(241) - flux(242) - flux(300) + flux(3) + flux(10) + flux(22) + flux(32) + flux(34) + flux(65) + flux(70) + flux(90) + flux(189) + flux(193) + flux(303)
      dy(idx_CH2j) =  - flux(42) - flux(62) - flux(63) - flux(64) - flux(118) - flux(176) - flux(183) - flux(301) + flux(43) + flux(46) + flux(91) + flux(116) + flux(131) + flux(144) + flux(164) + flux(166) + flux(168) + flux(170) + flux(172) + flux(213) + flux(219) + flux(230) + flux(238) + flux(239) + flux(240) + flux(241) + flux(242)
      dy(idx_CH3) =  - flux(3) - flux(12) - flux(32) - flux(32) - flux(33) - flux(34) - flux(35) - flux(133) - flux(142) - flux(214) - flux(302) - flux(303) + flux(4) + flux(11) + flux(22) + flux(24) + flux(28) + flux(28) + flux(29) + flux(36) + flux(39) + flux(69) + flux(72) + flux(84) + flux(92) + flux(145) + flux(191) + flux(194) + flux(195) + flux(196) + flux(197) + flux(304)
      dy(idx_CH3j) =  - flux(43) - flux(65) - flux(66) - flux(67) - flux(68) - flux(84) - flux(93) - flux(184) + flux(111) + flux(112) + flux(117) + flux(118) + flux(123) + flux(132) + flux(146) + flux(177) + flux(178) + flux(179) + flux(180) + flux(182) + flux(185) + flux(186) + flux(214) + flux(302)
      dy(idx_CH4) =  - flux(4) - flux(22) - flux(28) - flux(36) - flux(39) - flux(111) - flux(121) - flux(123) - flux(135) - flux(143) - flux(144) - flux(145) - flux(146) - flux(185) - flux(191) - flux(192) - flux(193) - flux(195) - flux(197) - flux(215) - flux(220) - flux(226) - flux(243) - flux(248) - flux(304) + flux(12) + flux(32) + flux(33) + flux(35) + flux(73) + flux(154) + flux(166) + flux(178) + flux(199) + flux(205) + flux(206) + flux(249)
      dy(idx_CH4j) =  - flux(69) - flux(70) - flux(112) - flux(122) - flux(186) - flux(191) - flux(194) - flux(196) - flux(249) + flux(113) + flux(133) + flux(215) + flux(220) + flux(226) + flux(243) + flux(248)
      dy(idx_CH5j) =  - flux(72) - flux(73) - flux(113) - flux(154) - flux(166) - flux(178) - flux(189) - flux(199) - flux(205) - flux(206) + flux(93) + flux(121) + flux(122) + flux(135) + flux(191) + flux(192)
      dy(idx_CHj) =  - flux(41) - flux(46) - flux(61) - flux(162) - flux(165) - flux(167) - flux(169) - flux(173) - flux(174) - flux(299) + flux(42) + flux(45) + flux(87) + flux(110) + flux(115) + flux(130) + flux(141) + flux(142) + flux(143) + flux(152) + flux(154) + flux(155) + flux(158) + flux(212) + flux(218) + flux(225) + flux(229) + flux(233) + flux(234) + flux(235) + flux(236) + flux(237) + flux(297) + flux(301)
      dy(idx_CO) =  - flux(7) - flux(47) - flux(57) - flux(128) - flux(138) - flux(150) - flux(196) - flux(202) - flux(206) - flux(208) - flux(223) - flux(279) - flux(295) + flux(18) + flux(19) + flux(21) + flux(23) + flux(26) + flux(27) + flux(31) + flux(48) + flux(82) + flux(94) + flux(158) + flux(160) + flux(172) + flux(182) + flux(204) + flux(210) + flux(217) + flux(231) + flux(236) + flux(241) + flux(246) + flux(248) + flux(251) + flux(253) + flux(255) + flux(280)
      dy(idx_CO_DUST) =  - flux(280) + flux(279)
      dy(idx_COj) =  - flux(81) - flux(129) - flux(171) - flux(181) - flux(197) - flux(203) - flux(209) - flux(217) - flux(231) - flux(236) - flux(241) - flux(246) - flux(248) - flux(251) - flux(253) - flux(255) + flux(57) + flux(95) + flux(153) + flux(159) + flux(161) + flux(162) + flux(163) + flux(165) + flux(173) + flux(223)
      dy(idx_Cj) =  - flux(45) - flux(48) - flux(87) - flux(91) - flux(95) - flux(153) - flux(156) - flux(159) - flux(160) - flux(229) - flux(230) - flux(260) - flux(272) - flux(291) + flux(41) + flux(56) + flux(139) + flux(140) + flux(150) + flux(231) + flux(232) + flux(261) + flux(271) + flux(290) + flux(299)
      dy(idx_Cjj) =  - flux(261) - flux(273) - flux(292) + flux(262) + flux(272) + flux(291)
      dy(idx_Cjjj) =  - flux(262) - flux(274) - flux(293) + flux(263) + flux(273) + flux(292)
      dy(idx_Cjjjj) =  - flux(263) + flux(274) + flux(293)
      dy(idx_E) =  - flux(58) - flux(59) - flux(60) - flux(61) - flux(62) - flux(63) - flux(64) - flux(65) - flux(66) - flux(67) - flux(68) - flux(69) - flux(70) - flux(71) - flux(72) - flux(73) - flux(74) - flux(75) - flux(76) - flux(77) - flux(78) - flux(79) - flux(80) - flux(81) - flux(82) - flux(83) - flux(84) - flux(102) - flux(256) - flux(257) - flux(258) - flux(259) - flux(260) - flux(261) - flux(262) - flux(263) - flux(264) - flux(265) - flux(266) - flux(267) - flux(268) - flux(269) - flux(270) - flux(271) - flux(272) - flux(273) - flux(274) - flux(275) - flux(276) - flux(277) - flux(278) + flux(102) + flux(108) + flux(268) + flux(268) + flux(269) + flux(269) + flux(270) + flux(270) + flux(271) + flux(271) + flux(272) + flux(272) + flux(273) + flux(273) + flux(274) + flux(274) + flux(275) + flux(275) + flux(276) + flux(276) + flux(277) + flux(277) + flux(278) + flux(278) + flux(283) + flux(284) + flux(285) + flux(286) + flux(287) + flux(288) + flux(289) + flux(290) + flux(291) + flux(292) + flux(293) + flux(297) + flux(302) + flux(308) + flux(310) + flux(311) + flux(313) + flux(314) + flux(315)
      dy(idx_H) =  - flux(1) - flux(2) - flux(3) - flux(4) - flux(5) - flux(6) - flux(7) - flux(8) - flux(41) - flux(42) - flux(43) - flux(50) - flux(51) - flux(85) - flux(86) - flux(87) - flux(88) - flux(89) - flux(97) - flux(98) - flux(99) - flux(100) - flux(101) - flux(109) - flux(109) - flux(112) - flux(113) - flux(211) - flux(217) - flux(268) - flux(283) - flux(314) + flux(9) + flux(10) + flux(11) + flux(12) + flux(13) + flux(14) + flux(18) + flux(21) + flux(26) + flux(26) + flux(37) + flux(44) + flux(45) + flux(46) + flux(49) + flux(52) + flux(53) + flux(54) + flux(58) + flux(58) + flux(59) + flux(59) + flux(59) + flux(60) + flux(61) + flux(62) + flux(63) + flux(63) + flux(65) + flux(67) + flux(67) + flux(68) + flux(69) + flux(70) + flux(70) + flux(71) + flux(73) + flux(74) + flux(74) + flux(76) + flux(77) + flux(78) + flux(78) + flux(80) + flux(82) + flux(97) + flux(97) + flux(97) + flux(98) + flux(98) + flux(99) + flux(99) + flux(100) + flux(100) + flux(101) + flux(102) + flux(102) + flux(103) + flux(103) + flux(104) + flux(105) + flux(106) + flux(114) + flux(115) + flux(116) + flux(117) + flux(118) + flux(119) + flux(120) + flux(121) + flux(122) + flux(123) + flux(124) + flux(125) + flux(126) + flux(127) + flux(128) + flux(129) + flux(139) + flux(141) + flux(143) + flux(146) + flux(147) + flux(149) + flux(153) + flux(156) + flux(162) + flux(163) + flux(176) + flux(187) + flux(188) + flux(212) + flux(213) + flux(214) + flux(215) + flux(216) + flux(256) + flux(294) + flux(294) + flux(296) + flux(298) + flux(299) + flux(300) + flux(301) + flux(303) + flux(304) + flux(305) + flux(307) + flux(311) + flux(312) + flux(312)
      dy(idx_H2) =  - flux(9) - flux(10) - flux(11) - flux(12) - flux(13) - flux(14) - flux(15) - flux(44) - flux(45) - flux(46) - flux(49) - flux(55) - flux(90) - flux(91) - flux(92) - flux(93) - flux(97) - flux(102) - flux(103) - flux(103) - flux(104) - flux(105) - flux(106) - flux(107) - flux(114) - flux(118) - flux(120) - flux(122) - flux(125) - flux(127) - flux(129) - flux(294) - flux(311) - flux(312) - flux(313) + flux(1) + flux(2) + flux(3) + flux(4) + flux(5) + flux(6) + flux(27) + flux(41) + flux(42) + flux(43) + flux(60) + flux(64) + flux(66) + flux(68) + flux(72) + flux(75) + flux(77) + flux(79) + flux(103) + flux(104) + flux(105) + flux(106) + flux(107) + flux(109) + flux(110) + flux(111) + flux(112) + flux(113) + flux(123) + flux(130) + flux(131) + flux(132) + flux(133) + flux(134) + flux(135) + flux(136) + flux(137) + flux(138) + flux(140) + flux(142) + flux(143) + flux(144) + flux(157) + flux(165) + flux(169) + flux(184) + flux(190) + flux(211) + flux(218) + flux(219) + flux(220) + flux(221) + flux(222) + flux(223) + flux(224)
      dy(idx_H2O) =  - flux(6) - flux(35) - flux(38) - flux(53) - flux(100) - flux(106) - flux(126) - flux(137) - flux(148) - flux(149) - flux(156) - flux(167) - flux(169) - flux(194) - flux(200) - flux(205) - flux(207) - flux(209) - flux(210) - flux(222) - flux(227) - flux(245) - flux(250) - flux(253) - flux(281) - flux(307) - flux(308) + flux(14) + flux(30) + flux(31) + flux(34) + flux(39) + flux(40) + flux(80) + flux(89) + flux(170) + flux(180) + flux(235) + flux(240) + flux(254) + flux(282)
      dy(idx_H2O_DUST) =  - flux(282) + flux(281)
      dy(idx_H2Oj) =  - flux(74) - flux(75) - flux(76) - flux(127) - flux(155) - flux(168) - flux(179) - flux(190) - flux(195) - flux(201) - flux(207) - flux(208) - flux(235) - flux(240) - flux(254) + flux(53) + flux(124) + flux(125) + flux(136) + flux(198) + flux(199) + flux(204) + flux(222) + flux(227) + flux(245) + flux(250) + flux(253) + flux(308)
      dy(idx_H2j) =  - flux(58) - flux(114) - flux(115) - flux(116) - flux(117) - flux(119) - flux(121) - flux(123) - flux(124) - flux(126) - flux(128) - flux(211) - flux(218) - flux(219) - flux(220) - flux(221) - flux(222) - flux(223) - flux(224) - flux(296) + flux(49) + flux(55) + flux(85) + flux(313)
      dy(idx_H3Oj) =  - flux(77) - flux(78) - flux(79) - flux(80) - flux(157) - flux(170) - flux(180) + flux(126) + flux(127) + flux(137) + flux(167) + flux(189) + flux(193) + flux(194) + flux(195) + flux(200) + flux(201) + flux(205) + flux(207) + flux(210)
      dy(idx_H3j) =  - flux(59) - flux(60) - flux(130) - flux(131) - flux(132) - flux(133) - flux(134) - flux(135) - flux(136) - flux(137) - flux(138) + flux(114)
      dy(idx_HCOj) =  - flux(82) - flux(158) - flux(172) - flux(182) - flux(204) - flux(210) + flux(108) + flux(128) + flux(129) + flux(138) + flux(156) + flux(157) + flux(169) + flux(171) + flux(174) + flux(175) + flux(176) + flux(181) + flux(183) + flux(184) + flux(196) + flux(197) + flux(202) + flux(203) + flux(206) + flux(208) + flux(209)
      dy(idx_He) =  - flux(269) - flux(284) - flux(315) + flux(44) + flux(47) + flux(50) + flux(55) + flux(56) + flux(139) + flux(140) + flux(141) + flux(142) + flux(143) + flux(144) + flux(145) + flux(146) + flux(147) + flux(148) + flux(149) + flux(150) + flux(151) + flux(225) + flux(226) + flux(227) + flux(228) + flux(257) + flux(258)
      dy(idx_Hej) =  - flux(44) - flux(47) - flux(50) - flux(55) - flux(56) - flux(139) - flux(140) - flux(141) - flux(142) - flux(143) - flux(144) - flux(145) - flux(146) - flux(147) - flux(148) - flux(149) - flux(150) - flux(151) - flux(225) - flux(226) - flux(227) - flux(228) - flux(257) - flux(258) - flux(270) - flux(285) + flux(259) + flux(269) + flux(284) + flux(315)
      dy(idx_Hejj) =  - flux(259) + flux(270) + flux(285)
      dy(idx_Hj) =  - flux(49) - flux(52) - flux(53) - flux(54) - flux(85) - flux(110) - flux(111) - flux(212) - flux(213) - flux(214) - flux(215) - flux(216) - flux(256) + flux(44) + flux(48) + flux(50) + flux(51) + flux(145) + flux(148) + flux(211) + flux(217) + flux(268) + flux(283) + flux(296) + flux(306) + flux(311) + flux(314)
      dy(idx_O) =  - flux(13) - flux(20) - flux(21) - flux(25) - flux(26) - flux(27) - flux(36) - flux(37) - flux(38) - flux(52) - flux(88) - flux(94) - flux(95) - flux(96) - flux(96) - flux(108) - flux(119) - flux(134) - flux(162) - flux(176) - flux(184) - flux(186) - flux(188) - flux(189) - flux(190) - flux(246) - flux(275) - flux(286) + flux(5) + flux(8) + flux(17) + flux(19) + flux(29) + flux(33) + flux(40) + flux(51) + flux(57) + flux(71) + flux(74) + flux(75) + flux(77) + flux(81) + flux(83) + flux(83) + flux(99) + flux(101) + flux(101) + flux(105) + flux(107) + flux(107) + flux(150) + flux(151) + flux(152) + flux(159) + flux(161) + flux(164) + flux(174) + flux(175) + flux(177) + flux(192) + flux(198) + flux(200) + flux(201) + flux(202) + flux(203) + flux(233) + flux(238) + flux(243) + flux(244) + flux(245) + flux(247) + flux(264) + flux(295) + flux(305) + flux(306) + flux(309) + flux(309)
      dy(idx_O2) =  - flux(8) - flux(15) - flux(19) - flux(23) - flux(31) - flux(54) - flux(101) - flux(107) - flux(151) - flux(159) - flux(160) - flux(173) - flux(174) - flux(183) - flux(224) - flux(228) - flux(247) - flux(249) - flux(252) - flux(254) - flux(255) - flux(309) - flux(310) + flux(37) + flux(96) + flux(232) + flux(237) + flux(242)
      dy(idx_O2j) =  - flux(83) - flux(161) - flux(175) - flux(232) - flux(237) - flux(242) + flux(54) + flux(187) + flux(188) + flux(190) + flux(224) + flux(228) + flux(247) + flux(249) + flux(252) + flux(254) + flux(255) + flux(310)
      dy(idx_OH) =  - flux(5) - flux(14) - flux(17) - flux(18) - flux(29) - flux(30) - flux(33) - flux(34) - flux(37) - flux(39) - flux(40) - flux(40) - flux(48) - flux(89) - flux(99) - flux(105) - flux(124) - flux(136) - flux(147) - flux(153) - flux(165) - flux(187) - flux(198) - flux(199) - flux(201) - flux(203) - flux(204) - flux(216) - flux(221) - flux(244) - flux(251) - flux(305) + flux(6) + flux(7) + flux(8) + flux(13) + flux(15) + flux(15) + flux(20) + flux(23) + flux(25) + flux(35) + flux(36) + flux(38) + flux(38) + flux(76) + flux(78) + flux(79) + flux(88) + flux(100) + flux(106) + flux(148) + flux(155) + flux(168) + flux(173) + flux(179) + flux(183) + flux(185) + flux(186) + flux(207) + flux(208) + flux(209) + flux(234) + flux(239) + flux(250) + flux(252) + flux(307)
      dy(idx_OHj) =  - flux(71) - flux(125) - flux(152) - flux(164) - flux(177) - flux(188) - flux(192) - flux(193) - flux(198) - flux(200) - flux(202) - flux(234) - flux(239) - flux(250) - flux(252) - flux(306) + flux(119) + flux(120) + flux(134) + flux(149) + flux(216) + flux(221) + flux(244) + flux(251)
      dy(idx_Oj) =  - flux(51) - flux(57) - flux(120) - flux(163) - flux(185) - flux(187) - flux(233) - flux(238) - flux(243) - flux(244) - flux(245) - flux(247) - flux(264) - flux(276) - flux(287) + flux(47) + flux(52) + flux(147) + flux(151) + flux(160) + flux(246) + flux(265) + flux(275) + flux(286)
      dy(idx_Ojj) =  - flux(265) - flux(277) - flux(288) + flux(266) + flux(276) + flux(287)
      dy(idx_Ojjj) =  - flux(266) - flux(278) - flux(289) + flux(267) + flux(277) + flux(288)
      dy(idx_Ojjjj) =  - flux(267) + flux(278) + flux(289)

      !! PREPROCESS_END
    end if

  end subroutine fex

  ! ***************************
  ! Jacobian, pd(i,j)=df(i)/dx(j), see DLSODES documentation
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use prizmo_commons
    implicit none
    integer::neq, j, ian, jan
    real*8::tt, n(neq), pdj(neq)

  end subroutine jes

  ! *******************
  function get_fex(x, Tgas) result(dy)
    use prizmo_commons
    implicit none
    integer,parameter::neq=nspecies+1
    real*8,intent(in)::x(nspecies), Tgas
    real*8::dy(neq), tt, y(neq)

    tt = 0d0

    y(1:nspecies) = x
    y(idx_Tgas) = Tgas

    call fex(neq, tt, y, dy)

  end function get_fex

end module prizmo_ode
