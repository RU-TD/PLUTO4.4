module prizmo_rates_heating
  use prizmo_commons
contains

  ! ************************
  subroutine compute_photorates_heating(x, Tgas, jflux)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, jflux(nphoto)
    real*8::f(nphoto)

    kall_heat = 0d0

    !! PREPROCESS_PHOTOHEATING_RATE

    ! H -> H+ + E
    f(:) = photo_xsecs(:, 283) * jflux * max(energy - energy_threshold(283), 0d0) / energy / hplanck
    kall_heat(283) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! He -> He+ + E
    f(:) = photo_xsecs(:, 284) * jflux * max(energy - energy_threshold(284), 0d0) / energy / hplanck
    kall_heat(284) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! He+ -> He++ + E
    f(:) = photo_xsecs(:, 285) * jflux * max(energy - energy_threshold(285), 0d0) / energy / hplanck
    kall_heat(285) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! O -> O+ + E
    f(:) = photo_xsecs(:, 286) * jflux * max(energy - energy_threshold(286), 0d0) / energy / hplanck
    kall_heat(286) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! O+ -> O++ + E
    f(:) = photo_xsecs(:, 287) * jflux * max(energy - energy_threshold(287), 0d0) / energy / hplanck
    kall_heat(287) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! O++ -> O+++ + E
    f(:) = photo_xsecs(:, 288) * jflux * max(energy - energy_threshold(288), 0d0) / energy / hplanck
    kall_heat(288) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! O+++ -> O++++ + E
    f(:) = photo_xsecs(:, 289) * jflux * max(energy - energy_threshold(289), 0d0) / energy / hplanck
    kall_heat(289) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! C -> C+ + E
    f(:) = photo_xsecs(:, 290) * jflux * max(energy - energy_threshold(290), 0d0) / energy / hplanck
    kall_heat(290) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! C+ -> C++ + E
    f(:) = photo_xsecs(:, 291) * jflux * max(energy - energy_threshold(291), 0d0) / energy / hplanck
    kall_heat(291) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! C++ -> C+++ + E
    f(:) = photo_xsecs(:, 292) * jflux * max(energy - energy_threshold(292), 0d0) / energy / hplanck
    kall_heat(292) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! C+++ -> C++++ + E
    f(:) = photo_xsecs(:, 293) * jflux * max(energy - energy_threshold(293), 0d0) / energy / hplanck
    kall_heat(293) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! H2 -> H + H
    ! skipped
    ! CO -> C + O
    f(:) = photo_xsecs(:, 295) * jflux * max(energy - energy_threshold(295), 0d0) / energy / hplanck
    kall_heat(295) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! H2+ -> H+ + H
    f(:) = photo_xsecs(:, 296) * jflux * max(energy - energy_threshold(296), 0d0) / energy / hplanck
    kall_heat(296) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! CH -> CH+ + E
    f(:) = photo_xsecs(:, 297) * jflux * max(energy - energy_threshold(297), 0d0) / energy / hplanck
    kall_heat(297) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! CH -> C + H
    f(:) = photo_xsecs(:, 298) * jflux * max(energy - energy_threshold(298), 0d0) / energy / hplanck
    kall_heat(298) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! CH+ -> C+ + H
    f(:) = photo_xsecs(:, 299) * jflux * max(energy - energy_threshold(299), 0d0) / energy / hplanck
    kall_heat(299) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! CH2 -> CH + H
    f(:) = photo_xsecs(:, 300) * jflux * max(energy - energy_threshold(300), 0d0) / energy / hplanck
    kall_heat(300) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! CH2+ -> CH+ + H
    f(:) = photo_xsecs(:, 301) * jflux * max(energy - energy_threshold(301), 0d0) / energy / hplanck
    kall_heat(301) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! CH3 -> CH3+ + E
    f(:) = photo_xsecs(:, 302) * jflux * max(energy - energy_threshold(302), 0d0) / energy / hplanck
    kall_heat(302) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! CH3 -> CH2 + H
    f(:) = photo_xsecs(:, 303) * jflux * max(energy - energy_threshold(303), 0d0) / energy / hplanck
    kall_heat(303) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! CH4 -> CH3 + H
    f(:) = photo_xsecs(:, 304) * jflux * max(energy - energy_threshold(304), 0d0) / energy / hplanck
    kall_heat(304) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! OH -> O + H
    f(:) = photo_xsecs(:, 305) * jflux * max(energy - energy_threshold(305), 0d0) / energy / hplanck
    kall_heat(305) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! OH+ -> O + H+
    f(:) = photo_xsecs(:, 306) * jflux * max(energy - energy_threshold(306), 0d0) / energy / hplanck
    kall_heat(306) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! H2O -> OH + H
    f(:) = photo_xsecs(:, 307) * jflux * max(energy - energy_threshold(307), 0d0) / energy / hplanck
    kall_heat(307) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! H2O -> H2O+ + E
    f(:) = photo_xsecs(:, 308) * jflux * max(energy - energy_threshold(308), 0d0) / energy / hplanck
    kall_heat(308) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! O2 -> O + O
    f(:) = photo_xsecs(:, 309) * jflux * max(energy - energy_threshold(309), 0d0) / energy / hplanck
    kall_heat(309) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    ! O2 -> O2+ + E
    f(:) = photo_xsecs(:, 310) * jflux * max(energy - energy_threshold(310), 0d0) / energy / hplanck
    kall_heat(310) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0

    !! PREPROCESS_END

  end subroutine compute_photorates_heating

end module prizmo_rates_heating
