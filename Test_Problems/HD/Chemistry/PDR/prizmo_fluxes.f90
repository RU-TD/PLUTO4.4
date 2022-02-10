module prizmo_fluxes
  !!BEGIN_MAX_SPECIES_NUMBER
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! max number of indexes, reactants and products
    integer, parameter::max_indexes = 6

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_MAX_SPECIES_NUMBER
contains

  ! ***********************
  subroutine print_sorted_fluxes(n, Tgas, jflux, nbest)
    use prizmo_commons
    use prizmo_utils
    implicit none
    integer,intent(in)::nbest
    integer::idx(nrea), i
    real*8,intent(in)::n(nmols), Tgas, jflux(nphoto)
    real*8::fluxes(nrea)
    character(len=max_character_len)::rnames(nrea)

    fluxes(:) = get_fluxes(n(:), Tgas, jflux(:))
    idx(:) = sort_fluxes(fluxes(:))
    rnames(:) = get_reaction_names()

    print '(2a5,a40,2a17)', "#", "idx", "reaction", "flux", "krate"
    do i=1,nbest
      write(*, '(2I5,a40,2E17.8e3)') i, idx(i), rnames(idx(i)), fluxes(idx(i)), kall(idx(i))
    end do

  end subroutine print_sorted_fluxes

  ! ***********************
  subroutine print_sorted_fluxes_species(n, Tgas, jflux, nbest, species)
    use prizmo_commons
    use prizmo_utils
    implicit none
    integer,intent(in)::nbest, species(:)
    integer::idx(nrea), i
    real*8,intent(in)::n(nmols), Tgas, jflux(nphoto)
    real*8::fluxes(nrea)
    character(len=max_character_len)::rnames(nrea)

    fluxes(:) = get_fluxes(n(:), Tgas, jflux(:))
    do i=1,nrea
      if(.not. reaction_has_any_reactants(i, species)) then
        fluxes(i) = 0d0
      end if
    end do
    idx(:) = sort_fluxes(fluxes(:))
    rnames(:) = get_reaction_names()

    print '(2a5,a40,2a17)', "#", "idx", "reaction", "flux", "krate"
    do i=1,min(nbest, nrea)
      write(*, '(2I5,a40,2E17.8e3)') i, idx(i), rnames(idx(i)), fluxes(idx(i)), kall(idx(i))
    end do

  end subroutine print_sorted_fluxes_species

  ! *********************
  function get_sorted_fluxes(n, Tgas, jflux) result(fluxes_sorted)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(nmols), Tgas, jflux(nphoto)
    real*8::fluxes_sorted(nrea)
    real*8::fluxes(nrea)
    integer::idx(nrea), i

    fluxes(:) = get_fluxes(n(:), Tgas, jflux(:))
    idx(:) = sort_fluxes(fluxes(:))

    do i=1,nrea
      fluxes_sorted(i) = fluxes(idx(i))
    end do

  end function get_sorted_fluxes

  ! *********************
  function get_sorted_fluxes_species(n, Tgas, jflux, species) result(fluxes_sorted)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(nmols), Tgas, jflux(nphoto)
    real*8::fluxes_sorted(nrea)
    real*8::fluxes(nrea)
    integer::idx(nrea), i
    integer,intent(in)::species(:)

    fluxes(:) = get_fluxes(n(:), Tgas, jflux(:))
    do i=1,nrea
      if(.not. reaction_has_any_reactants(i, species)) then
        fluxes(i) = 0d0
      end if
    end do
    idx(:) = sort_fluxes(fluxes(:))

    do i=1,nrea
      fluxes_sorted(i) = fluxes(idx(i))
    end do

  end function get_sorted_fluxes_species


  ! **********************
  function sort_fluxes(fluxes_in) result(idx)
    use prizmo_commons
    implicit none
    real*8,intent(in)::fluxes_in(nrea)
    real*8::tmp, fluxes(nrea)
    integer::idx(nrea), i, itmp
    logical::swap

    fluxes(:) = fluxes_in(:)

    ! fill indexes
    do i=1,nrea
      idx(i) = i
    end do

    ! loop until there is nothing to swap
    do
      ! swap flag
      swap = .false.
      ! loop on elements to swap if not sorted
      do i=2,nrea
        ! if not sorted swap
        if(fluxes(i) > fluxes(i-1)) then
          swap = .true.
          ! swap fluxes
          tmp = fluxes(i)
          fluxes(i) = fluxes(i-1)
          fluxes(i-1) = tmp
          ! swap indeces
          itmp = idx(i)
          idx(i) = idx(i-1)
          idx(i-1) = itmp
        end if
      end do
      if(swap.eqv..false.) exit
    end do


  end function sort_fluxes

  ! **********************
  function get_fluxes(n, Tgas, jflux) result(fluxes)
    use prizmo_commons
    use prizmo_tdust
    use prizmo_rates
    use prizmo_rates_photo
    implicit none
    real*8,intent(in)::n(nmols), Tgas
    real*8::fluxes(nrea), Tdust, jflux(nphoto)

    Tdust = get_Tdust(n(:), Tgas, jflux(:))

    ! compute non-photochemical rates
    call compute_rates(n(:), Tgas, Tdust, jflux(:))

    ! compute photochemical rates
    call compute_rates_photo(n(:), Tgas, jflux(:))

    !!BEGIN_FLUXES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    fluxes(1) = n(idx_H)*n(idx_CH)
    fluxes(2) = n(idx_H)*n(idx_CH3)
    fluxes(3) = n(idx_H)*n(idx_CH4)
    fluxes(4) = n(idx_H)*n(idx_OH)
    fluxes(5) = n(idx_H)*n(idx_H2O)
    fluxes(6) = n(idx_H)*n(idx_CO)
    fluxes(7) = n(idx_H)*n(idx_O2)
    fluxes(8) = n(idx_H2)*n(idx_C)
    fluxes(9) = n(idx_H2)*n(idx_CH)
    fluxes(10) = n(idx_H2)*n(idx_CH2)
    fluxes(11) = n(idx_H2)*n(idx_CH3)
    fluxes(12) = n(idx_H2)*n(idx_O)
    fluxes(13) = n(idx_H2)*n(idx_OH)
    fluxes(14) = n(idx_H2)*n(idx_O2)
    fluxes(15) = n(idx_C)*n(idx_CH2)
    fluxes(16) = n(idx_C)*n(idx_OH)
    fluxes(17) = n(idx_C)*n(idx_OH)
    fluxes(18) = n(idx_CH)*n(idx_O)
    fluxes(19) = n(idx_CH)*n(idx_CH4)
    fluxes(20) = n(idx_CH2)*n(idx_CH2)
    fluxes(21) = n(idx_CH2)*n(idx_O)
    fluxes(22) = n(idx_CH2)*n(idx_CH4)
    fluxes(23) = n(idx_CH2)*n(idx_OH)
    fluxes(24) = n(idx_CH2)*n(idx_OH)
    fluxes(25) = n(idx_CH2)*n(idx_O2)
    fluxes(26) = n(idx_CH3)*n(idx_CH3)
    fluxes(27) = n(idx_CH3)*n(idx_OH)
    fluxes(28) = n(idx_CH3)*n(idx_OH)
    fluxes(29) = n(idx_CH3)*n(idx_H2O)
    fluxes(30) = n(idx_O)*n(idx_CH4)
    fluxes(31) = n(idx_O)*n(idx_OH)
    fluxes(32) = n(idx_O)*n(idx_H2O)
    fluxes(33) = n(idx_CH4)*n(idx_OH)
    fluxes(34) = n(idx_OH)*n(idx_OH)
    fluxes(35) = n(idx_H)*n(idx_CH2j)
    fluxes(36) = n(idx_H)*n(idx_CH3j)
    fluxes(37) = n(idx_H2)*n(idx_Hej)
    fluxes(38) = n(idx_H2)*n(idx_Cj)
    fluxes(39) = n(idx_Hej)*n(idx_CO)
    fluxes(40) = n(idx_Hj)*n(idx_H2)
    fluxes(41) = n(idx_H)*n(idx_Hej)
    fluxes(42) = n(idx_H)*n(idx_Oj)
    fluxes(43) = n(idx_Hj)*n(idx_O)
    fluxes(44) = n(idx_Hej)*n(idx_C)
    fluxes(45) = n(idx_Oj)*n(idx_CO)
    fluxes(46) = n(idx_H2j)*n(idx_E)
    fluxes(47) = n(idx_H3j)*n(idx_E)
    fluxes(48) = n(idx_H3j)*n(idx_E)
    fluxes(49) = n(idx_CHj)*n(idx_E)
    fluxes(50) = n(idx_CH2j)*n(idx_E)
    fluxes(51) = n(idx_CH2j)*n(idx_E)
    fluxes(52) = n(idx_CH2j)*n(idx_E)
    fluxes(53) = n(idx_CH3j)*n(idx_E)
    fluxes(54) = n(idx_CH3j)*n(idx_E)
    fluxes(55) = n(idx_CH3j)*n(idx_E)
    fluxes(56) = n(idx_CH4j)*n(idx_E)
    fluxes(57) = n(idx_CH4j)*n(idx_E)
    fluxes(58) = n(idx_OHj)*n(idx_E)
    fluxes(59) = n(idx_CH5j)*n(idx_E)
    fluxes(60) = n(idx_CH5j)*n(idx_E)
    fluxes(61) = n(idx_H2Oj)*n(idx_E)
    fluxes(62) = n(idx_H2Oj)*n(idx_E)
    fluxes(63) = n(idx_H2Oj)*n(idx_E)
    fluxes(64) = n(idx_H3Oj)*n(idx_E)
    fluxes(65) = n(idx_H3Oj)*n(idx_E)
    fluxes(66) = n(idx_H3Oj)*n(idx_E)
    fluxes(67) = n(idx_H3Oj)*n(idx_E)
    fluxes(68) = n(idx_COj)*n(idx_E)
    fluxes(69) = n(idx_HCOj)*n(idx_E)
    fluxes(70) = n(idx_O2j)*n(idx_E)
    fluxes(71) = n(idx_Hj)*n(idx_E)
    fluxes(72) = n(idx_Hej)*n(idx_E)
    fluxes(73) = n(idx_Cj)*n(idx_E)
    fluxes(74) = n(idx_CH3j)*n(idx_E)
    fluxes(75) = n(idx_Oj)*n(idx_E)
    fluxes(76) = n(idx_Hj)*n(idx_H)
    fluxes(77) = n(idx_H)*n(idx_O)
    fluxes(78) = n(idx_H)*n(idx_OH)
    fluxes(79) = n(idx_H2)*n(idx_Cj)
    fluxes(80) = n(idx_H2)*n(idx_CH)
    fluxes(81) = n(idx_H2)*n(idx_CH3j)
    fluxes(82) = n(idx_O)*n(idx_O)
    fluxes(83) = n(idx_H)*n(idx_H2)
    fluxes(84) = n(idx_H)*n(idx_CH)
    fluxes(85) = n(idx_H)*n(idx_OH)
    fluxes(86) = n(idx_H)*n(idx_H2O)
    fluxes(87) = n(idx_H)*n(idx_O2)
    fluxes(88) = n(idx_H2)*n(idx_E)
    fluxes(89) = n(idx_H2)*n(idx_H2)
    fluxes(90) = n(idx_H2)*n(idx_CH)
    fluxes(91) = n(idx_H2)*n(idx_OH)
    fluxes(92) = n(idx_H2)*n(idx_H2O)
    fluxes(93) = n(idx_H2)*n(idx_O2)
    fluxes(94) = n(idx_CH)*n(idx_O)
    fluxes(95) = n(idx_H)*n(idx_CH2)
    fluxes(96) = n(idx_C)*n(idx_O2)
    fluxes(97) = n(idx_CH)*n(idx_O)
    fluxes(98) = n(idx_CH)*n(idx_O2)
    fluxes(99) = n(idx_CH2)*n(idx_O)
    fluxes(100) = n(idx_CH2)*n(idx_O)
    fluxes(101) = n(idx_H)*n(idx_CHj)
    fluxes(102) = n(idx_H2)*n(idx_CHj)
    fluxes(103) = n(idx_Cj)*n(idx_OH)
    fluxes(104) = n(idx_Hj)*n(idx_H2O)
    fluxes(105) = n(idx_Hj)*n(idx_O2)
    fluxes(106) = n(idx_H2)*n(idx_Hej)
    fluxes(107) = n(idx_H)*n(idx_C)
    fluxes(108) = n(idx_H)*n(idx_Cj)
    fluxes(109) = n(idx_H2)*n(idx_C)
    fluxes(110) = n(idx_C)*n(idx_O)
    fluxes(111) = n(idx_Cj)*n(idx_O)
    fluxes(112) = n(idx_Hj)*n(idx_CH2)
    fluxes(113) = n(idx_Hj)*n(idx_CH4)
    fluxes(114) = n(idx_H)*n(idx_CH4j)
    fluxes(115) = n(idx_H)*n(idx_CH5j)
    fluxes(116) = n(idx_H2j)*n(idx_H2)
    fluxes(117) = n(idx_H2j)*n(idx_C)
    fluxes(118) = n(idx_H2j)*n(idx_CH)
    fluxes(119) = n(idx_H2j)*n(idx_CH2)
    fluxes(120) = n(idx_H2)*n(idx_CH2j)
    fluxes(121) = n(idx_H2j)*n(idx_O)
    fluxes(122) = n(idx_H2)*n(idx_Oj)
    fluxes(123) = n(idx_H2j)*n(idx_CH4)
    fluxes(124) = n(idx_H2)*n(idx_CH4j)
    fluxes(125) = n(idx_H2j)*n(idx_CH4)
    fluxes(126) = n(idx_H2j)*n(idx_OH)
    fluxes(127) = n(idx_H2)*n(idx_OHj)
    fluxes(128) = n(idx_H2j)*n(idx_H2O)
    fluxes(129) = n(idx_H2)*n(idx_H2Oj)
    fluxes(130) = n(idx_H2j)*n(idx_CO)
    fluxes(131) = n(idx_H2)*n(idx_COj)
    fluxes(132) = n(idx_H3j)*n(idx_C)
    fluxes(133) = n(idx_H3j)*n(idx_CH)
    fluxes(134) = n(idx_H3j)*n(idx_CH2)
    fluxes(135) = n(idx_H3j)*n(idx_CH3)
    fluxes(136) = n(idx_H3j)*n(idx_O)
    fluxes(137) = n(idx_H3j)*n(idx_CH4)
    fluxes(138) = n(idx_H3j)*n(idx_OH)
    fluxes(139) = n(idx_H3j)*n(idx_H2O)
    fluxes(140) = n(idx_H3j)*n(idx_CO)
    fluxes(141) = n(idx_Hej)*n(idx_CH)
    fluxes(142) = n(idx_Hej)*n(idx_CH2)
    fluxes(143) = n(idx_Hej)*n(idx_CH2)
    fluxes(144) = n(idx_Hej)*n(idx_CH3)
    fluxes(145) = n(idx_Hej)*n(idx_CH4)
    fluxes(146) = n(idx_Hej)*n(idx_CH4)
    fluxes(147) = n(idx_Hej)*n(idx_CH4)
    fluxes(148) = n(idx_Hej)*n(idx_CH4)
    fluxes(149) = n(idx_Hej)*n(idx_OH)
    fluxes(150) = n(idx_Hej)*n(idx_H2O)
    fluxes(151) = n(idx_Hej)*n(idx_H2O)
    fluxes(152) = n(idx_Hej)*n(idx_CO)
    fluxes(153) = n(idx_Hej)*n(idx_O2)
    fluxes(154) = n(idx_C)*n(idx_OHj)
    fluxes(155) = n(idx_Cj)*n(idx_OH)
    fluxes(156) = n(idx_C)*n(idx_CH5j)
    fluxes(157) = n(idx_C)*n(idx_H2Oj)
    fluxes(158) = n(idx_Cj)*n(idx_H2O)
    fluxes(159) = n(idx_C)*n(idx_H3Oj)
    fluxes(160) = n(idx_C)*n(idx_HCOj)
    fluxes(161) = n(idx_Cj)*n(idx_O2)
    fluxes(162) = n(idx_Cj)*n(idx_O2)
    fluxes(163) = n(idx_C)*n(idx_O2j)
    fluxes(164) = n(idx_CHj)*n(idx_O)
    fluxes(165) = n(idx_CH)*n(idx_Oj)
    fluxes(166) = n(idx_CH)*n(idx_OHj)
    fluxes(167) = n(idx_CHj)*n(idx_OH)
    fluxes(168) = n(idx_CH)*n(idx_CH5j)
    fluxes(169) = n(idx_CHj)*n(idx_H2O)
    fluxes(170) = n(idx_CH)*n(idx_H2Oj)
    fluxes(171) = n(idx_CHj)*n(idx_H2O)
    fluxes(172) = n(idx_CH)*n(idx_H3Oj)
    fluxes(173) = n(idx_CH)*n(idx_COj)
    fluxes(174) = n(idx_CH)*n(idx_HCOj)
    fluxes(175) = n(idx_CHj)*n(idx_O2)
    fluxes(176) = n(idx_CHj)*n(idx_O2)
    fluxes(177) = n(idx_CH)*n(idx_O2j)
    fluxes(178) = n(idx_CH2j)*n(idx_O)
    fluxes(179) = n(idx_CH2)*n(idx_OHj)
    fluxes(180) = n(idx_CH2)*n(idx_CH5j)
    fluxes(181) = n(idx_CH2)*n(idx_H2Oj)
    fluxes(182) = n(idx_CH2)*n(idx_H3Oj)
    fluxes(183) = n(idx_CH2)*n(idx_COj)
    fluxes(184) = n(idx_CH2)*n(idx_HCOj)
    fluxes(185) = n(idx_CH2j)*n(idx_O2)
    fluxes(186) = n(idx_CH3j)*n(idx_O)
    fluxes(187) = n(idx_Oj)*n(idx_CH4)
    fluxes(188) = n(idx_O)*n(idx_CH4j)
    fluxes(189) = n(idx_Oj)*n(idx_OH)
    fluxes(190) = n(idx_O)*n(idx_OHj)
    fluxes(191) = n(idx_O)*n(idx_CH5j)
    fluxes(192) = n(idx_O)*n(idx_H2Oj)
    fluxes(193) = n(idx_CH4j)*n(idx_CH4)
    fluxes(194) = n(idx_CH4)*n(idx_OHj)
    fluxes(195) = n(idx_CH4)*n(idx_OHj)
    fluxes(196) = n(idx_CH4j)*n(idx_H2O)
    fluxes(197) = n(idx_CH4)*n(idx_H2Oj)
    fluxes(198) = n(idx_CH4j)*n(idx_CO)
    fluxes(199) = n(idx_CH4)*n(idx_COj)
    fluxes(200) = n(idx_OHj)*n(idx_OH)
    fluxes(201) = n(idx_OH)*n(idx_CH5j)
    fluxes(202) = n(idx_OHj)*n(idx_H2O)
    fluxes(203) = n(idx_OH)*n(idx_H2Oj)
    fluxes(204) = n(idx_OHj)*n(idx_CO)
    fluxes(205) = n(idx_OH)*n(idx_COj)
    fluxes(206) = n(idx_OH)*n(idx_HCOj)
    fluxes(207) = n(idx_CH5j)*n(idx_H2O)
    fluxes(208) = n(idx_CH5j)*n(idx_CO)
    fluxes(209) = n(idx_H2Oj)*n(idx_H2O)
    fluxes(210) = n(idx_H2Oj)*n(idx_CO)
    fluxes(211) = n(idx_H2O)*n(idx_COj)
    fluxes(212) = n(idx_H2O)*n(idx_HCOj)
    fluxes(213) = n(idx_H)*n(idx_H2j)
    fluxes(214) = n(idx_Hj)*n(idx_CH)
    fluxes(215) = n(idx_Hj)*n(idx_CH2)
    fluxes(216) = n(idx_Hj)*n(idx_CH3)
    fluxes(217) = n(idx_Hj)*n(idx_CH4)
    fluxes(218) = n(idx_Hj)*n(idx_OH)
    fluxes(219) = n(idx_H)*n(idx_COj)
    fluxes(220) = n(idx_H2j)*n(idx_CH)
    fluxes(221) = n(idx_H2j)*n(idx_CH2)
    fluxes(222) = n(idx_H2j)*n(idx_CH4)
    fluxes(223) = n(idx_H2j)*n(idx_OH)
    fluxes(224) = n(idx_H2j)*n(idx_H2O)
    fluxes(225) = n(idx_H2j)*n(idx_CO)
    fluxes(226) = n(idx_H2j)*n(idx_O2)
    fluxes(227) = n(idx_Hej)*n(idx_CH)
    fluxes(228) = n(idx_Hej)*n(idx_CH4)
    fluxes(229) = n(idx_Hej)*n(idx_H2O)
    fluxes(230) = n(idx_Hej)*n(idx_O2)
    fluxes(231) = n(idx_Cj)*n(idx_CH)
    fluxes(232) = n(idx_Cj)*n(idx_CH2)
    fluxes(233) = n(idx_C)*n(idx_COj)
    fluxes(234) = n(idx_C)*n(idx_O2j)
    fluxes(235) = n(idx_CH)*n(idx_Oj)
    fluxes(236) = n(idx_CH)*n(idx_OHj)
    fluxes(237) = n(idx_CH)*n(idx_H2Oj)
    fluxes(238) = n(idx_CH)*n(idx_COj)
    fluxes(239) = n(idx_CH)*n(idx_O2j)
    fluxes(240) = n(idx_CH2)*n(idx_Oj)
    fluxes(241) = n(idx_CH2)*n(idx_OHj)
    fluxes(242) = n(idx_CH2)*n(idx_H2Oj)
    fluxes(243) = n(idx_CH2)*n(idx_COj)
    fluxes(244) = n(idx_CH2)*n(idx_O2j)
    fluxes(245) = n(idx_Oj)*n(idx_CH4)
    fluxes(246) = n(idx_Oj)*n(idx_OH)
    fluxes(247) = n(idx_Oj)*n(idx_H2O)
    fluxes(248) = n(idx_O)*n(idx_COj)
    fluxes(249) = n(idx_Oj)*n(idx_O2)
    fluxes(250) = n(idx_CH4)*n(idx_COj)
    fluxes(251) = n(idx_CH4j)*n(idx_O2)
    fluxes(252) = n(idx_OHj)*n(idx_H2O)
    fluxes(253) = n(idx_OH)*n(idx_COj)
    fluxes(254) = n(idx_OHj)*n(idx_O2)
    fluxes(255) = n(idx_H2O)*n(idx_COj)
    fluxes(256) = n(idx_H2Oj)*n(idx_O2)
    fluxes(257) = n(idx_COj)*n(idx_O2)
    fluxes(258) = n(idx_H2)
    fluxes(259) = n(idx_CO)
    fluxes(260) = n(idx_H2j)
    fluxes(261) = n(idx_H3j)
    fluxes(262) = n(idx_H3j)
    fluxes(263) = n(idx_C)
    fluxes(264) = n(idx_CH)
    fluxes(265) = n(idx_CH)
    fluxes(266) = n(idx_CHj)
    fluxes(267) = n(idx_CH2)
    fluxes(268) = n(idx_CH2)
    fluxes(269) = n(idx_CH2j)
    fluxes(270) = n(idx_CH3)
    fluxes(271) = n(idx_CH3)
    fluxes(272) = n(idx_CH3)
    fluxes(273) = n(idx_CH4)
    fluxes(274) = n(idx_CH4)
    fluxes(275) = n(idx_CH4)
    fluxes(276) = n(idx_OH)
    fluxes(277) = n(idx_OH)
    fluxes(278) = n(idx_OHj)
    fluxes(279) = n(idx_H2O)
    fluxes(280) = n(idx_H2O)
    fluxes(281) = n(idx_O2)
    fluxes(282) = n(idx_O2)
    fluxes(283) = n(idx_H2)
    fluxes(284) = n(idx_H2)
    fluxes(285) = n(idx_H2)
    fluxes(286) = n(idx_H)
    fluxes(287) = n(idx_He)
    fluxes(288) = n(idx_CH3j)*n(idx_E)
    fluxes(289) = n(idx_H)*n(idx_H)

    ! multiply rate coefficient (vectorized)
    fluxes(:) = fluxes(:) * kall(:)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_FLUXES

  end function get_fluxes

  ! **********************
  ! check if the reaction with index idx_reaction has any of the species listed in species(:)
  function reaction_has_any_reactants(idx_reaction, species) result(has_any)
    use prizmo_commons
    implicit none
    integer,intent(in)::species(:), idx_reaction
    logical::has_any
    integer::i, j, indexes(nrea, max_indexes)

    ! default is not found
    has_any = .false.

    indexes(:, :) = get_species_indexes()

    ! loop on indexes
    do i=1,size(species)
      ! loop on reactants
      do j=1,max_indexes
        ! if found return with true
        if(indexes(idx_reaction, j) == species(i)) then
          has_any = .true.
          return
        end if
      end do
    end do

  end function reaction_has_any_reactants

  ! **********************
  function get_species_indexes() result(indexes)
    use prizmo_commons
    implicit none
    real*8::indexes(nrea, max_indexes)

    !!BEGIN_SPECIES_INDEXES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    indexes(1, :) = (/idx_H, idx_CH, idx_C, idx_H2, -1, -1/)
    indexes(2, :) = (/idx_H, idx_CH3, idx_CH2, idx_H2, -1, -1/)
    indexes(3, :) = (/idx_H, idx_CH4, idx_CH3, idx_H2, -1, -1/)
    indexes(4, :) = (/idx_H, idx_OH, idx_O, idx_H2, -1, -1/)
    indexes(5, :) = (/idx_H, idx_H2O, idx_OH, idx_H2, -1, -1/)
    indexes(6, :) = (/idx_H, idx_CO, idx_OH, idx_C, -1, -1/)
    indexes(7, :) = (/idx_H, idx_O2, idx_OH, idx_O, -1, -1/)
    indexes(8, :) = (/idx_H2, idx_C, idx_CH, idx_H, -1, -1/)
    indexes(9, :) = (/idx_H2, idx_CH, idx_CH2, idx_H, -1, -1/)
    indexes(10, :) = (/idx_H2, idx_CH2, idx_CH3, idx_H, -1, -1/)
    indexes(11, :) = (/idx_H2, idx_CH3, idx_CH4, idx_H, -1, -1/)
    indexes(12, :) = (/idx_H2, idx_O, idx_OH, idx_H, -1, -1/)
    indexes(13, :) = (/idx_H2, idx_OH, idx_H2O, idx_H, -1, -1/)
    indexes(14, :) = (/idx_H2, idx_O2, idx_OH, idx_OH, -1, -1/)
    indexes(15, :) = (/idx_C, idx_CH2, idx_CH, idx_CH, -1, -1/)
    indexes(16, :) = (/idx_C, idx_OH, idx_O, idx_CH, -1, -1/)
    indexes(17, :) = (/idx_C, idx_OH, idx_CO, idx_H, -1, -1/)
    indexes(18, :) = (/idx_CH, idx_O, idx_OH, idx_C, -1, -1/)
    indexes(19, :) = (/idx_CH, idx_CH4, idx_CH3, idx_CH2, -1, -1/)
    indexes(20, :) = (/idx_CH2, idx_CH2, idx_CH3, idx_CH, -1, -1/)
    indexes(21, :) = (/idx_CH2, idx_O, idx_OH, idx_CH, -1, -1/)
    indexes(22, :) = (/idx_CH2, idx_CH4, idx_CH3, idx_CH3, -1, -1/)
    indexes(23, :) = (/idx_CH2, idx_OH, idx_O, idx_CH3, -1, -1/)
    indexes(24, :) = (/idx_CH2, idx_OH, idx_H2O, idx_CH, -1, -1/)
    indexes(25, :) = (/idx_CH2, idx_O2, idx_CO, idx_H2O, -1, -1/)
    indexes(26, :) = (/idx_CH3, idx_CH3, idx_CH4, idx_CH2, -1, -1/)
    indexes(27, :) = (/idx_CH3, idx_OH, idx_CH4, idx_O, -1, -1/)
    indexes(28, :) = (/idx_CH3, idx_OH, idx_H2O, idx_CH2, -1, -1/)
    indexes(29, :) = (/idx_CH3, idx_H2O, idx_OH, idx_CH4, -1, -1/)
    indexes(30, :) = (/idx_O, idx_CH4, idx_OH, idx_CH3, -1, -1/)
    indexes(31, :) = (/idx_O, idx_OH, idx_O2, idx_H, -1, -1/)
    indexes(32, :) = (/idx_O, idx_H2O, idx_OH, idx_OH, -1, -1/)
    indexes(33, :) = (/idx_CH4, idx_OH, idx_H2O, idx_CH3, -1, -1/)
    indexes(34, :) = (/idx_OH, idx_OH, idx_H2O, idx_O, -1, -1/)
    indexes(35, :) = (/idx_H, idx_CH2j, idx_CHj, idx_H2, -1, -1/)
    indexes(36, :) = (/idx_H, idx_CH3j, idx_CH2j, idx_H2, -1, -1/)
    indexes(37, :) = (/idx_H2, idx_Hej, idx_He, idx_Hj, idx_H, -1/)
    indexes(38, :) = (/idx_H2, idx_Cj, idx_CHj, idx_H, -1, -1/)
    indexes(39, :) = (/idx_Hej, idx_CO, idx_Oj, idx_C, idx_He, -1/)
    indexes(40, :) = (/idx_Hj, idx_H2, idx_H2j, idx_H, -1, -1/)
    indexes(41, :) = (/idx_H, idx_Hej, idx_He, idx_Hj, -1, -1/)
    indexes(42, :) = (/idx_H, idx_Oj, idx_O, idx_Hj, -1, -1/)
    indexes(43, :) = (/idx_Hj, idx_O, idx_Oj, idx_H, -1, -1/)
    indexes(44, :) = (/idx_Hej, idx_C, idx_Cj, idx_He, -1, -1/)
    indexes(45, :) = (/idx_Oj, idx_CO, idx_COj, idx_O, -1, -1/)
    indexes(46, :) = (/idx_H2j, idx_E, idx_H, idx_H, -1, -1/)
    indexes(47, :) = (/idx_H3j, idx_E, idx_H, idx_H, idx_H, -1/)
    indexes(48, :) = (/idx_H3j, idx_E, idx_H2, idx_H, -1, -1/)
    indexes(49, :) = (/idx_CHj, idx_E, idx_C, idx_H, -1, -1/)
    indexes(50, :) = (/idx_CH2j, idx_E, idx_CH, idx_H, -1, -1/)
    indexes(51, :) = (/idx_CH2j, idx_E, idx_C, idx_H, idx_H, -1/)
    indexes(52, :) = (/idx_CH2j, idx_E, idx_C, idx_H2, -1, -1/)
    indexes(53, :) = (/idx_CH3j, idx_E, idx_CH2, idx_H, -1, -1/)
    indexes(54, :) = (/idx_CH3j, idx_E, idx_CH, idx_H2, -1, -1/)
    indexes(55, :) = (/idx_CH3j, idx_E, idx_CH, idx_H, idx_H, -1/)
    indexes(56, :) = (/idx_CH4j, idx_E, idx_CH3, idx_H, -1, -1/)
    indexes(57, :) = (/idx_CH4j, idx_E, idx_CH2, idx_H, idx_H, -1/)
    indexes(58, :) = (/idx_OHj, idx_E, idx_O, idx_H, -1, -1/)
    indexes(59, :) = (/idx_CH5j, idx_E, idx_CH3, idx_H2, -1, -1/)
    indexes(60, :) = (/idx_CH5j, idx_E, idx_CH4, idx_H, -1, -1/)
    indexes(61, :) = (/idx_H2Oj, idx_E, idx_O, idx_H, idx_H, -1/)
    indexes(62, :) = (/idx_H2Oj, idx_E, idx_O, idx_H2, -1, -1/)
    indexes(63, :) = (/idx_H2Oj, idx_E, idx_OH, idx_H, -1, -1/)
    indexes(64, :) = (/idx_H3Oj, idx_E, idx_O, idx_H2, idx_H, -1/)
    indexes(65, :) = (/idx_H3Oj, idx_E, idx_OH, idx_H, idx_H, -1/)
    indexes(66, :) = (/idx_H3Oj, idx_E, idx_OH, idx_H2, -1, -1/)
    indexes(67, :) = (/idx_H3Oj, idx_E, idx_H2O, idx_H, -1, -1/)
    indexes(68, :) = (/idx_COj, idx_E, idx_O, idx_C, -1, -1/)
    indexes(69, :) = (/idx_HCOj, idx_E, idx_CO, idx_H, -1, -1/)
    indexes(70, :) = (/idx_O2j, idx_E, idx_O, idx_O, -1, -1/)
    indexes(71, :) = (/idx_Hj, idx_E, idx_H, -1, -1, -1/)
    indexes(72, :) = (/idx_Hej, idx_E, idx_He, -1, -1, -1/)
    indexes(73, :) = (/idx_Cj, idx_E, idx_C, -1, -1, -1/)
    indexes(74, :) = (/idx_CH3j, idx_E, idx_CH3, -1, -1, -1/)
    indexes(75, :) = (/idx_Oj, idx_E, idx_O, -1, -1, -1/)
    indexes(76, :) = (/idx_Hj, idx_H, idx_H2j, -1, -1, -1/)
    indexes(77, :) = (/idx_H, idx_O, idx_OH, -1, -1, -1/)
    indexes(78, :) = (/idx_H, idx_OH, idx_H2O, -1, -1, -1/)
    indexes(79, :) = (/idx_H2, idx_Cj, idx_CH2j, -1, -1, -1/)
    indexes(80, :) = (/idx_H2, idx_CH, idx_CH3, -1, -1, -1/)
    indexes(81, :) = (/idx_H2, idx_CH3j, idx_CH5j, -1, -1, -1/)
    indexes(82, :) = (/idx_O, idx_O, idx_O2, -1, -1, -1/)
    indexes(83, :) = (/idx_H, idx_H2, idx_H, idx_H, idx_H, -1/)
    indexes(84, :) = (/idx_H, idx_CH, idx_C, idx_H, idx_H, -1/)
    indexes(85, :) = (/idx_H, idx_OH, idx_O, idx_H, idx_H, -1/)
    indexes(86, :) = (/idx_H, idx_H2O, idx_OH, idx_H, idx_H, -1/)
    indexes(87, :) = (/idx_H, idx_O2, idx_O, idx_O, idx_H, -1/)
    indexes(88, :) = (/idx_H2, idx_E, idx_H, idx_H, idx_E, -1/)
    indexes(89, :) = (/idx_H2, idx_H2, idx_H2, idx_H, idx_H, -1/)
    indexes(90, :) = (/idx_H2, idx_CH, idx_C, idx_H2, idx_H, -1/)
    indexes(91, :) = (/idx_H2, idx_OH, idx_O, idx_H2, idx_H, -1/)
    indexes(92, :) = (/idx_H2, idx_H2O, idx_OH, idx_H2, idx_H, -1/)
    indexes(93, :) = (/idx_H2, idx_O2, idx_O, idx_O, idx_H2, -1/)
    indexes(94, :) = (/idx_CH, idx_O, idx_HCOj, idx_E, -1, -1/)
    indexes(95, :) = (/idx_H, idx_CH2, idx_CH, idx_H2, -1, -1/)
    indexes(96, :) = (/idx_C, idx_O2, idx_CO, idx_O, -1, -1/)
    indexes(97, :) = (/idx_CH, idx_O, idx_CO, idx_H, -1, -1/)
    indexes(98, :) = (/idx_CH, idx_O2, idx_CO, idx_OH, -1, -1/)
    indexes(99, :) = (/idx_CH2, idx_O, idx_CO, idx_H, idx_H, -1/)
    indexes(100, :) = (/idx_CH2, idx_O, idx_CO, idx_H2, -1, -1/)
    indexes(101, :) = (/idx_H, idx_CHj, idx_Cj, idx_H2, -1, -1/)
    indexes(102, :) = (/idx_H2, idx_CHj, idx_CH2j, idx_H, -1, -1/)
    indexes(103, :) = (/idx_Cj, idx_OH, idx_CO, idx_Hj, -1, -1/)
    indexes(104, :) = (/idx_Hj, idx_H2O, idx_H2Oj, idx_H, -1, -1/)
    indexes(105, :) = (/idx_Hj, idx_O2, idx_O2j, idx_H, -1, -1/)
    indexes(106, :) = (/idx_H2, idx_Hej, idx_He, idx_H2j, -1, -1/)
    indexes(107, :) = (/idx_H, idx_C, idx_CH, -1, -1, -1/)
    indexes(108, :) = (/idx_H, idx_Cj, idx_CHj, -1, -1, -1/)
    indexes(109, :) = (/idx_H2, idx_C, idx_CH2, -1, -1, -1/)
    indexes(110, :) = (/idx_C, idx_O, idx_CO, -1, -1, -1/)
    indexes(111, :) = (/idx_Cj, idx_O, idx_COj, -1, -1, -1/)
    indexes(112, :) = (/idx_Hj, idx_CH2, idx_CHj, idx_H2, -1, -1/)
    indexes(113, :) = (/idx_Hj, idx_CH4, idx_CH3j, idx_H2, -1, -1/)
    indexes(114, :) = (/idx_H, idx_CH4j, idx_CH3j, idx_H2, -1, -1/)
    indexes(115, :) = (/idx_H, idx_CH5j, idx_CH4j, idx_H2, -1, -1/)
    indexes(116, :) = (/idx_H2j, idx_H2, idx_H3j, idx_H, -1, -1/)
    indexes(117, :) = (/idx_H2j, idx_C, idx_CHj, idx_H, -1, -1/)
    indexes(118, :) = (/idx_H2j, idx_CH, idx_CH2j, idx_H, -1, -1/)
    indexes(119, :) = (/idx_H2j, idx_CH2, idx_CH3j, idx_H, -1, -1/)
    indexes(120, :) = (/idx_H2, idx_CH2j, idx_CH3j, idx_H, -1, -1/)
    indexes(121, :) = (/idx_H2j, idx_O, idx_OHj, idx_H, -1, -1/)
    indexes(122, :) = (/idx_H2, idx_Oj, idx_OHj, idx_H, -1, -1/)
    indexes(123, :) = (/idx_H2j, idx_CH4, idx_CH5j, idx_H, -1, -1/)
    indexes(124, :) = (/idx_H2, idx_CH4j, idx_CH5j, idx_H, -1, -1/)
    indexes(125, :) = (/idx_H2j, idx_CH4, idx_CH3j, idx_H2, idx_H, -1/)
    indexes(126, :) = (/idx_H2j, idx_OH, idx_H2Oj, idx_H, -1, -1/)
    indexes(127, :) = (/idx_H2, idx_OHj, idx_H2Oj, idx_H, -1, -1/)
    indexes(128, :) = (/idx_H2j, idx_H2O, idx_H3Oj, idx_H, -1, -1/)
    indexes(129, :) = (/idx_H2, idx_H2Oj, idx_H3Oj, idx_H, -1, -1/)
    indexes(130, :) = (/idx_H2j, idx_CO, idx_HCOj, idx_H, -1, -1/)
    indexes(131, :) = (/idx_H2, idx_COj, idx_HCOj, idx_H, -1, -1/)
    indexes(132, :) = (/idx_H3j, idx_C, idx_CHj, idx_H2, -1, -1/)
    indexes(133, :) = (/idx_H3j, idx_CH, idx_CH2j, idx_H2, -1, -1/)
    indexes(134, :) = (/idx_H3j, idx_CH2, idx_CH3j, idx_H2, -1, -1/)
    indexes(135, :) = (/idx_H3j, idx_CH3, idx_CH4j, idx_H2, -1, -1/)
    indexes(136, :) = (/idx_H3j, idx_O, idx_OHj, idx_H2, -1, -1/)
    indexes(137, :) = (/idx_H3j, idx_CH4, idx_CH5j, idx_H2, -1, -1/)
    indexes(138, :) = (/idx_H3j, idx_OH, idx_H2Oj, idx_H2, -1, -1/)
    indexes(139, :) = (/idx_H3j, idx_H2O, idx_H3Oj, idx_H2, -1, -1/)
    indexes(140, :) = (/idx_H3j, idx_CO, idx_HCOj, idx_H2, -1, -1/)
    indexes(141, :) = (/idx_Hej, idx_CH, idx_Cj, idx_He, idx_H, -1/)
    indexes(142, :) = (/idx_Hej, idx_CH2, idx_Cj, idx_He, idx_H2, -1/)
    indexes(143, :) = (/idx_Hej, idx_CH2, idx_CHj, idx_He, idx_H, -1/)
    indexes(144, :) = (/idx_Hej, idx_CH3, idx_CHj, idx_He, idx_H2, -1/)
    indexes(145, :) = (/idx_Hej, idx_CH4, idx_CHj, idx_He, idx_H2, idx_H/)
    indexes(146, :) = (/idx_Hej, idx_CH4, idx_CH2j, idx_He, idx_H2, -1/)
    indexes(147, :) = (/idx_Hej, idx_CH4, idx_CH3, idx_He, idx_Hj, -1/)
    indexes(148, :) = (/idx_Hej, idx_CH4, idx_CH3j, idx_He, idx_H, -1/)
    indexes(149, :) = (/idx_Hej, idx_OH, idx_Oj, idx_He, idx_H, -1/)
    indexes(150, :) = (/idx_Hej, idx_H2O, idx_OH, idx_He, idx_Hj, -1/)
    indexes(151, :) = (/idx_Hej, idx_H2O, idx_OHj, idx_He, idx_H, -1/)
    indexes(152, :) = (/idx_Hej, idx_CO, idx_O, idx_Cj, idx_He, -1/)
    indexes(153, :) = (/idx_Hej, idx_O2, idx_Oj, idx_O, idx_He, -1/)
    indexes(154, :) = (/idx_C, idx_OHj, idx_O, idx_CHj, -1, -1/)
    indexes(155, :) = (/idx_Cj, idx_OH, idx_COj, idx_H, -1, -1/)
    indexes(156, :) = (/idx_C, idx_CH5j, idx_CH4, idx_CHj, -1, -1/)
    indexes(157, :) = (/idx_C, idx_H2Oj, idx_OH, idx_CHj, -1, -1/)
    indexes(158, :) = (/idx_Cj, idx_H2O, idx_HCOj, idx_H, -1, -1/)
    indexes(159, :) = (/idx_C, idx_H3Oj, idx_HCOj, idx_H2, -1, -1/)
    indexes(160, :) = (/idx_C, idx_HCOj, idx_CO, idx_CHj, -1, -1/)
    indexes(161, :) = (/idx_Cj, idx_O2, idx_COj, idx_O, -1, -1/)
    indexes(162, :) = (/idx_Cj, idx_O2, idx_CO, idx_Oj, -1, -1/)
    indexes(163, :) = (/idx_C, idx_O2j, idx_COj, idx_O, -1, -1/)
    indexes(164, :) = (/idx_CHj, idx_O, idx_COj, idx_H, -1, -1/)
    indexes(165, :) = (/idx_CH, idx_Oj, idx_COj, idx_H, -1, -1/)
    indexes(166, :) = (/idx_CH, idx_OHj, idx_O, idx_CH2j, -1, -1/)
    indexes(167, :) = (/idx_CHj, idx_OH, idx_COj, idx_H2, -1, -1/)
    indexes(168, :) = (/idx_CH, idx_CH5j, idx_CH4, idx_CH2j, -1, -1/)
    indexes(169, :) = (/idx_CHj, idx_H2O, idx_H3Oj, idx_C, -1, -1/)
    indexes(170, :) = (/idx_CH, idx_H2Oj, idx_OH, idx_CH2j, -1, -1/)
    indexes(171, :) = (/idx_CHj, idx_H2O, idx_HCOj, idx_H2, -1, -1/)
    indexes(172, :) = (/idx_CH, idx_H3Oj, idx_H2O, idx_CH2j, -1, -1/)
    indexes(173, :) = (/idx_CH, idx_COj, idx_HCOj, idx_C, -1, -1/)
    indexes(174, :) = (/idx_CH, idx_HCOj, idx_CO, idx_CH2j, -1, -1/)
    indexes(175, :) = (/idx_CHj, idx_O2, idx_COj, idx_OH, -1, -1/)
    indexes(176, :) = (/idx_CHj, idx_O2, idx_HCOj, idx_O, -1, -1/)
    indexes(177, :) = (/idx_CH, idx_O2j, idx_HCOj, idx_O, -1, -1/)
    indexes(178, :) = (/idx_CH2j, idx_O, idx_HCOj, idx_H, -1, -1/)
    indexes(179, :) = (/idx_CH2, idx_OHj, idx_O, idx_CH3j, -1, -1/)
    indexes(180, :) = (/idx_CH2, idx_CH5j, idx_CH4, idx_CH3j, -1, -1/)
    indexes(181, :) = (/idx_CH2, idx_H2Oj, idx_OH, idx_CH3j, -1, -1/)
    indexes(182, :) = (/idx_CH2, idx_H3Oj, idx_H2O, idx_CH3j, -1, -1/)
    indexes(183, :) = (/idx_CH2, idx_COj, idx_HCOj, idx_CH, -1, -1/)
    indexes(184, :) = (/idx_CH2, idx_HCOj, idx_CO, idx_CH3j, -1, -1/)
    indexes(185, :) = (/idx_CH2j, idx_O2, idx_HCOj, idx_OH, -1, -1/)
    indexes(186, :) = (/idx_CH3j, idx_O, idx_HCOj, idx_H2, -1, -1/)
    indexes(187, :) = (/idx_Oj, idx_CH4, idx_OH, idx_CH3j, -1, -1/)
    indexes(188, :) = (/idx_O, idx_CH4j, idx_OH, idx_CH3j, -1, -1/)
    indexes(189, :) = (/idx_Oj, idx_OH, idx_O2j, idx_H, -1, -1/)
    indexes(190, :) = (/idx_O, idx_OHj, idx_O2j, idx_H, -1, -1/)
    indexes(191, :) = (/idx_O, idx_CH5j, idx_H3Oj, idx_CH2, -1, -1/)
    indexes(192, :) = (/idx_O, idx_H2Oj, idx_O2j, idx_H2, -1, -1/)
    indexes(193, :) = (/idx_CH4j, idx_CH4, idx_CH5j, idx_CH3, -1, -1/)
    indexes(194, :) = (/idx_CH4, idx_OHj, idx_CH5j, idx_O, -1, -1/)
    indexes(195, :) = (/idx_CH4, idx_OHj, idx_H3Oj, idx_CH2, -1, -1/)
    indexes(196, :) = (/idx_CH4j, idx_H2O, idx_H3Oj, idx_CH3, -1, -1/)
    indexes(197, :) = (/idx_CH4, idx_H2Oj, idx_H3Oj, idx_CH3, -1, -1/)
    indexes(198, :) = (/idx_CH4j, idx_CO, idx_HCOj, idx_CH3, -1, -1/)
    indexes(199, :) = (/idx_CH4, idx_COj, idx_HCOj, idx_CH3, -1, -1/)
    indexes(200, :) = (/idx_OHj, idx_OH, idx_H2Oj, idx_O, -1, -1/)
    indexes(201, :) = (/idx_OH, idx_CH5j, idx_H2Oj, idx_CH4, -1, -1/)
    indexes(202, :) = (/idx_OHj, idx_H2O, idx_H3Oj, idx_O, -1, -1/)
    indexes(203, :) = (/idx_OH, idx_H2Oj, idx_H3Oj, idx_O, -1, -1/)
    indexes(204, :) = (/idx_OHj, idx_CO, idx_HCOj, idx_O, -1, -1/)
    indexes(205, :) = (/idx_OH, idx_COj, idx_HCOj, idx_O, -1, -1/)
    indexes(206, :) = (/idx_OH, idx_HCOj, idx_CO, idx_H2Oj, -1, -1/)
    indexes(207, :) = (/idx_CH5j, idx_H2O, idx_H3Oj, idx_CH4, -1, -1/)
    indexes(208, :) = (/idx_CH5j, idx_CO, idx_HCOj, idx_CH4, -1, -1/)
    indexes(209, :) = (/idx_H2Oj, idx_H2O, idx_H3Oj, idx_OH, -1, -1/)
    indexes(210, :) = (/idx_H2Oj, idx_CO, idx_HCOj, idx_OH, -1, -1/)
    indexes(211, :) = (/idx_H2O, idx_COj, idx_HCOj, idx_OH, -1, -1/)
    indexes(212, :) = (/idx_H2O, idx_HCOj, idx_CO, idx_H3Oj, -1, -1/)
    indexes(213, :) = (/idx_H, idx_H2j, idx_H2, idx_Hj, -1, -1/)
    indexes(214, :) = (/idx_Hj, idx_CH, idx_CHj, idx_H, -1, -1/)
    indexes(215, :) = (/idx_Hj, idx_CH2, idx_CH2j, idx_H, -1, -1/)
    indexes(216, :) = (/idx_Hj, idx_CH3, idx_CH3j, idx_H, -1, -1/)
    indexes(217, :) = (/idx_Hj, idx_CH4, idx_CH4j, idx_H, -1, -1/)
    indexes(218, :) = (/idx_Hj, idx_OH, idx_OHj, idx_H, -1, -1/)
    indexes(219, :) = (/idx_H, idx_COj, idx_CO, idx_Hj, -1, -1/)
    indexes(220, :) = (/idx_H2j, idx_CH, idx_CHj, idx_H2, -1, -1/)
    indexes(221, :) = (/idx_H2j, idx_CH2, idx_CH2j, idx_H2, -1, -1/)
    indexes(222, :) = (/idx_H2j, idx_CH4, idx_CH4j, idx_H2, -1, -1/)
    indexes(223, :) = (/idx_H2j, idx_OH, idx_OHj, idx_H2, -1, -1/)
    indexes(224, :) = (/idx_H2j, idx_H2O, idx_H2Oj, idx_H2, -1, -1/)
    indexes(225, :) = (/idx_H2j, idx_CO, idx_COj, idx_H2, -1, -1/)
    indexes(226, :) = (/idx_H2j, idx_O2, idx_O2j, idx_H2, -1, -1/)
    indexes(227, :) = (/idx_Hej, idx_CH, idx_CHj, idx_He, -1, -1/)
    indexes(228, :) = (/idx_Hej, idx_CH4, idx_CH4j, idx_He, -1, -1/)
    indexes(229, :) = (/idx_Hej, idx_H2O, idx_H2Oj, idx_He, -1, -1/)
    indexes(230, :) = (/idx_Hej, idx_O2, idx_O2j, idx_He, -1, -1/)
    indexes(231, :) = (/idx_Cj, idx_CH, idx_CHj, idx_C, -1, -1/)
    indexes(232, :) = (/idx_Cj, idx_CH2, idx_CH2j, idx_C, -1, -1/)
    indexes(233, :) = (/idx_C, idx_COj, idx_CO, idx_Cj, -1, -1/)
    indexes(234, :) = (/idx_C, idx_O2j, idx_O2, idx_Cj, -1, -1/)
    indexes(235, :) = (/idx_CH, idx_Oj, idx_O, idx_CHj, -1, -1/)
    indexes(236, :) = (/idx_CH, idx_OHj, idx_OH, idx_CHj, -1, -1/)
    indexes(237, :) = (/idx_CH, idx_H2Oj, idx_H2O, idx_CHj, -1, -1/)
    indexes(238, :) = (/idx_CH, idx_COj, idx_CO, idx_CHj, -1, -1/)
    indexes(239, :) = (/idx_CH, idx_O2j, idx_O2, idx_CHj, -1, -1/)
    indexes(240, :) = (/idx_CH2, idx_Oj, idx_O, idx_CH2j, -1, -1/)
    indexes(241, :) = (/idx_CH2, idx_OHj, idx_OH, idx_CH2j, -1, -1/)
    indexes(242, :) = (/idx_CH2, idx_H2Oj, idx_H2O, idx_CH2j, -1, -1/)
    indexes(243, :) = (/idx_CH2, idx_COj, idx_CO, idx_CH2j, -1, -1/)
    indexes(244, :) = (/idx_CH2, idx_O2j, idx_O2, idx_CH2j, -1, -1/)
    indexes(245, :) = (/idx_Oj, idx_CH4, idx_CH4j, idx_O, -1, -1/)
    indexes(246, :) = (/idx_Oj, idx_OH, idx_OHj, idx_O, -1, -1/)
    indexes(247, :) = (/idx_Oj, idx_H2O, idx_H2Oj, idx_O, -1, -1/)
    indexes(248, :) = (/idx_O, idx_COj, idx_CO, idx_Oj, -1, -1/)
    indexes(249, :) = (/idx_Oj, idx_O2, idx_O2j, idx_O, -1, -1/)
    indexes(250, :) = (/idx_CH4, idx_COj, idx_CO, idx_CH4j, -1, -1/)
    indexes(251, :) = (/idx_CH4j, idx_O2, idx_O2j, idx_CH4, -1, -1/)
    indexes(252, :) = (/idx_OHj, idx_H2O, idx_H2Oj, idx_OH, -1, -1/)
    indexes(253, :) = (/idx_OH, idx_COj, idx_CO, idx_OHj, -1, -1/)
    indexes(254, :) = (/idx_OHj, idx_O2, idx_O2j, idx_OH, -1, -1/)
    indexes(255, :) = (/idx_H2O, idx_COj, idx_CO, idx_H2Oj, -1, -1/)
    indexes(256, :) = (/idx_H2Oj, idx_O2, idx_O2j, idx_H2O, -1, -1/)
    indexes(257, :) = (/idx_COj, idx_O2, idx_O2j, idx_CO, -1, -1/)
    indexes(258, :) = (/idx_H2, idx_H, idx_H, -1, -1, -1/)
    indexes(259, :) = (/idx_CO, idx_C, idx_O, -1, -1, -1/)
    indexes(260, :) = (/idx_H2j, idx_Hj, idx_H, -1, -1, -1/)
    indexes(261, :) = (/idx_H3j, idx_H2j, idx_H, -1, -1, -1/)
    indexes(262, :) = (/idx_H3j, idx_H2, idx_Hj, -1, -1, -1/)
    indexes(263, :) = (/idx_C, idx_Cj, idx_E, -1, -1, -1/)
    indexes(264, :) = (/idx_CH, idx_CHj, idx_E, -1, -1, -1/)
    indexes(265, :) = (/idx_CH, idx_C, idx_H, -1, -1, -1/)
    indexes(266, :) = (/idx_CHj, idx_Cj, idx_H, -1, -1, -1/)
    indexes(267, :) = (/idx_CH2, idx_CH2j, idx_E, -1, -1, -1/)
    indexes(268, :) = (/idx_CH2, idx_CH, idx_H, -1, -1, -1/)
    indexes(269, :) = (/idx_CH2j, idx_CHj, idx_H, -1, -1, -1/)
    indexes(270, :) = (/idx_CH3, idx_CH3j, idx_E, -1, -1, -1/)
    indexes(271, :) = (/idx_CH3, idx_CH2, idx_H, -1, -1, -1/)
    indexes(272, :) = (/idx_CH3, idx_CH, idx_H2, -1, -1, -1/)
    indexes(273, :) = (/idx_CH4, idx_CH3, idx_H, -1, -1, -1/)
    indexes(274, :) = (/idx_CH4, idx_CH2, idx_H2, -1, -1, -1/)
    indexes(275, :) = (/idx_CH4, idx_CH, idx_H2, idx_H, -1, -1/)
    indexes(276, :) = (/idx_OH, idx_OHj, idx_E, -1, -1, -1/)
    indexes(277, :) = (/idx_OH, idx_O, idx_H, -1, -1, -1/)
    indexes(278, :) = (/idx_OHj, idx_O, idx_Hj, -1, -1, -1/)
    indexes(279, :) = (/idx_H2O, idx_OH, idx_H, -1, -1, -1/)
    indexes(280, :) = (/idx_H2O, idx_H2Oj, idx_E, -1, -1, -1/)
    indexes(281, :) = (/idx_O2, idx_O, idx_O, -1, -1, -1/)
    indexes(282, :) = (/idx_O2, idx_O2j, idx_E, -1, -1, -1/)
    indexes(283, :) = (/idx_H2, idx_Hj, idx_H, idx_E, -1, -1/)
    indexes(284, :) = (/idx_H2, idx_H, idx_H, -1, -1, -1/)
    indexes(285, :) = (/idx_H2, idx_H2j, idx_E, -1, -1, -1/)
    indexes(286, :) = (/idx_H, idx_Hj, idx_E, -1, -1, -1/)
    indexes(287, :) = (/idx_He, idx_Hej, idx_E, -1, -1, -1/)
    indexes(288, :) = (/idx_CH3j, idx_E, idx_C, idx_H2, idx_H, -1/)
    indexes(289, :) = (/idx_H, idx_H, idx_H2, -1, -1, -1/)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_SPECIES_INDEXES

  end function get_species_indexes


end module prizmo_fluxes
