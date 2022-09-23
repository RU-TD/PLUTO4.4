module prizmo_cooling_chemical
  use prizmo_commons
  use prizmo_utils
contains

  ! ***************
  function cooling_chemical(x, fluxes, Tgas) result(cool)
    implicit none
    real*8,intent(in)::x(nspecies), fluxes(nreactions), Tgas
    real*8::coola(nreactions), cool

    coola = get_cooling_chemical_array(x, fluxes, Tgas)

    ! call ranker(abs(coola), 5)  ! DEBUG

    cool = sum(coola)

  end function cooling_chemical

  ! ***********************
  function get_cooling_chemical_array(x, fluxes, Tgas) result(coola)
    implicit none
    real*8,intent(in)::x(nspecies), fluxes(nreactions), Tgas
    real*8::coola(nreactions)

    coola = 0d0

    !! PREPROCESS_RECOMBINATION

    ! H2+ + E -> H + H
    coola(58) = fluxes(58)
    
    ! H3+ + E -> H + H + H
    coola(59) = fluxes(59)
    
    ! H3+ + E -> H2 + H
    coola(60) = fluxes(60)
    
    ! CH+ + E -> C + H
    coola(61) = fluxes(61)
    
    ! CH2+ + E -> CH + H
    coola(62) = fluxes(62)
    
    ! CH2+ + E -> C + H + H
    coola(63) = fluxes(63)
    
    ! CH2+ + E -> C + H2
    coola(64) = fluxes(64)
    
    ! CH3+ + E -> CH2 + H
    coola(65) = fluxes(65)
    
    ! CH3+ + E -> CH + H2
    coola(66) = fluxes(66)
    
    ! CH3+ + E -> CH + H + H
    coola(67) = fluxes(67)
    
    ! CH3+ + E -> C + H2 + H
    coola(68) = fluxes(68)
    
    ! CH4+ + E -> CH3 + H
    coola(69) = fluxes(69)
    
    ! CH4+ + E -> CH2 + H + H
    coola(70) = fluxes(70)
    
    ! OH+ + E -> O + H
    coola(71) = fluxes(71)
    
    ! CH5+ + E -> CH3 + H2
    coola(72) = fluxes(72)
    
    ! CH5+ + E -> CH4 + H
    coola(73) = fluxes(73)
    
    ! H2O+ + E -> O + H + H
    coola(74) = fluxes(74)
    
    ! H2O+ + E -> O + H2
    coola(75) = fluxes(75)
    
    ! H2O+ + E -> OH + H
    coola(76) = fluxes(76)
    
    ! H3O+ + E -> O + H2 + H
    coola(77) = fluxes(77)
    
    ! H3O+ + E -> OH + H + H
    coola(78) = fluxes(78)
    
    ! H3O+ + E -> OH + H2
    coola(79) = fluxes(79)
    
    ! H3O+ + E -> H2O + H
    coola(80) = fluxes(80)
    
    ! CO+ + E -> O + C
    coola(81) = fluxes(81)
    
    ! HCO+ + E -> CO + H
    coola(82) = fluxes(82)
    
    ! O2+ + E -> O + O
    coola(83) = fluxes(83)
    
    ! CH3+ + E -> CH3
    coola(84) = fluxes(84)
    
    ! H+ + E -> H
    coola(256) = fluxes(256)
    
    ! He+ + E -> He
    coola(257) = fluxes(257)
    
    ! He+ + E -> He
    coola(258) = fluxes(258)
    
    ! He++ + E -> He+
    coola(259) = fluxes(259)
    
    ! C+ + E -> C
    coola(260) = fluxes(260)
    
    ! C++ + E -> C+
    coola(261) = fluxes(261)
    
    ! C+++ + E -> C++
    coola(262) = fluxes(262)
    
    ! C++++ + E -> C+++
    coola(263) = fluxes(263)
    
    ! O+ + E -> O
    coola(264) = fluxes(264)
    
    ! O++ + E -> O+
    coola(265) = fluxes(265)
    
    ! O+++ + E -> O++
    coola(266) = fluxes(266)
    
    ! O++++ + E -> O+++
    coola(267) = fluxes(267)

    !! PREPROCESS_END

    coola = coola * kboltzmann * Tgas

    !! PREPROCESS_CHEMICAL

    ! H + H2 -> H + H + H
    coola(97) = (7.177754e-12) * fluxes(97)
    
    ! H2 + H2 -> H2 + H + H
    coola(103) = (7.177754e-12) * fluxes(103)
    
    ! H + H -> H2
    coola(109) = (-6.729145e-12) * fluxes(109)

    !! PREPROCESS_END

  end function get_cooling_chemical_array

end module prizmo_cooling_chemical
