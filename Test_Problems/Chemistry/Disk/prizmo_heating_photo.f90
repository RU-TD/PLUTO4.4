module prizmo_heating_photo
  use prizmo_commons
contains

  ! ***************
  function heating_photo(x, Tgas, Tdust, jflux, ntot) result(heat)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto), ntot
    real*8::heat

    heat = 0d0

    !! PREPROCESS_PHOTOHEATING

    ! H -> H+ + E
    heat = heat + kall_heat(283) * x(idx_H)
    ! He -> He+ + E
    heat = heat + kall_heat(284) * x(idx_He)
    ! He+ -> He++ + E
    heat = heat + kall_heat(285) * x(idx_Hej)
    ! O -> O+ + E
    heat = heat + kall_heat(286) * x(idx_O)
    ! O+ -> O++ + E
    heat = heat + kall_heat(287) * x(idx_Oj)
    ! O++ -> O+++ + E
    heat = heat + kall_heat(288) * x(idx_Ojj)
    ! O+++ -> O++++ + E
    heat = heat + kall_heat(289) * x(idx_Ojjj)
    ! C -> C+ + E
    heat = heat + kall_heat(290) * x(idx_C)
    ! C+ -> C++ + E
    heat = heat + kall_heat(291) * x(idx_Cj)
    ! C++ -> C+++ + E
    heat = heat + kall_heat(292) * x(idx_Cjj)
    ! C+++ -> C++++ + E
    heat = heat + kall_heat(293) * x(idx_Cjjj)
    
    ! CO -> C + O
    heat = heat + kall_heat(295) * x(idx_CO)
    ! H2+ -> H+ + H
    heat = heat + kall_heat(296) * x(idx_H2j)
    ! CH -> CH+ + E
    heat = heat + kall_heat(297) * x(idx_CH)
    ! CH -> C + H
    heat = heat + kall_heat(298) * x(idx_CH)
    ! CH+ -> C+ + H
    heat = heat + kall_heat(299) * x(idx_CHj)
    ! CH2 -> CH + H
    heat = heat + kall_heat(300) * x(idx_CH2)
    ! CH2+ -> CH+ + H
    heat = heat + kall_heat(301) * x(idx_CH2j)
    ! CH3 -> CH3+ + E
    heat = heat + kall_heat(302) * x(idx_CH3)
    ! CH3 -> CH2 + H
    heat = heat + kall_heat(303) * x(idx_CH3)
    ! CH4 -> CH3 + H
    heat = heat + kall_heat(304) * x(idx_CH4)
    ! OH -> O + H
    heat = heat + kall_heat(305) * x(idx_OH)
    ! OH+ -> O + H+
    heat = heat + kall_heat(306) * x(idx_OHj)
    ! H2O -> OH + H
    heat = heat + kall_heat(307) * x(idx_H2O)
    ! H2O -> H2O+ + E
    heat = heat + kall_heat(308) * x(idx_H2O)
    ! O2 -> O + O
    heat = heat + kall_heat(309) * x(idx_O2)
    ! O2 -> O2+ + E
    heat = heat + kall_heat(310) * x(idx_O2)

    !! PREPROCESS_END

  end function heating_photo

end module prizmo_heating_photo
