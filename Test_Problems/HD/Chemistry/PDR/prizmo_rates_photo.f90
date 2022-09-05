module prizmo_rates_photo
  use prizmo_commons
  use prizmo_shielding
contains

  ! ************************
  subroutine compute_photorates(x, Tgas, jflux)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, jflux(nphoto)
    real*8::f(nphoto), log_tgas, log_NH2, log_NCO

    log_tgas = log10(Tgas)
    log_NH2 = log10(radial_Ncol_H2 + 1d-40)
    log_NCO = log10(radial_Ncol_CO + 1d-40)

    !! PREPROCESS_PHOTORATES

    ! H -> H+ + E
    f(:) = photo_xsecs(:, 283) * jflux / energy / hplanck
    kall(283) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! He -> He+ + E
    f(:) = photo_xsecs(:, 284) * jflux / energy / hplanck
    kall(284) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! He+ -> He++ + E
    f(:) = photo_xsecs(:, 285) * jflux / energy / hplanck
    kall(285) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! O -> O+ + E
    f(:) = photo_xsecs(:, 286) * jflux / energy / hplanck
    kall(286) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! O+ -> O++ + E
    f(:) = photo_xsecs(:, 287) * jflux / energy / hplanck
    kall(287) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! O++ -> O+++ + E
    f(:) = photo_xsecs(:, 288) * jflux / energy / hplanck
    kall(288) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! O+++ -> O++++ + E
    f(:) = photo_xsecs(:, 289) * jflux / energy / hplanck
    kall(289) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! C -> C+ + E
    f(:) = photo_xsecs(:, 290) * jflux / energy / hplanck
    kall(290) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! C+ -> C++ + E
    f(:) = photo_xsecs(:, 291) * jflux / energy / hplanck
    kall(291) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! C++ -> C+++ + E
    f(:) = photo_xsecs(:, 292) * jflux / energy / hplanck
    kall(292) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! C+++ -> C++++ + E
    f(:) = photo_xsecs(:, 293) * jflux / energy / hplanck
    kall(293) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! H2 -> H + H
    f(:) = photo_xsecs(:, 294) * jflux / energy / hplanck
    kall(294) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0  * shielding_H2(log_NH2, log_tgas)
    
    ! CO -> C + O
    f(:) = photo_xsecs(:, 295) * jflux / energy / hplanck
    kall(295) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0  * shielding_CO(log_NH2, log_NCO)
    
    ! H2+ -> H+ + H
    f(:) = photo_xsecs(:, 296) * jflux / energy / hplanck
    kall(296) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! CH -> CH+ + E
    f(:) = photo_xsecs(:, 297) * jflux / energy / hplanck
    kall(297) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! CH -> C + H
    f(:) = photo_xsecs(:, 298) * jflux / energy / hplanck
    kall(298) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! CH+ -> C+ + H
    f(:) = photo_xsecs(:, 299) * jflux / energy / hplanck
    kall(299) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! CH2 -> CH + H
    f(:) = photo_xsecs(:, 300) * jflux / energy / hplanck
    kall(300) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! CH2+ -> CH+ + H
    f(:) = photo_xsecs(:, 301) * jflux / energy / hplanck
    kall(301) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! CH3 -> CH3+ + E
    f(:) = photo_xsecs(:, 302) * jflux / energy / hplanck
    kall(302) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! CH3 -> CH2 + H
    f(:) = photo_xsecs(:, 303) * jflux / energy / hplanck
    kall(303) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! CH4 -> CH3 + H
    f(:) = photo_xsecs(:, 304) * jflux / energy / hplanck
    kall(304) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! OH -> O + H
    f(:) = photo_xsecs(:, 305) * jflux / energy / hplanck
    kall(305) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! OH+ -> O + H+
    f(:) = photo_xsecs(:, 306) * jflux / energy / hplanck
    kall(306) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! H2O -> OH + H
    f(:) = photo_xsecs(:, 307) * jflux / energy / hplanck
    kall(307) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! H2O -> H2O+ + E
    f(:) = photo_xsecs(:, 308) * jflux / energy / hplanck
    kall(308) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! O2 -> O + O
    f(:) = photo_xsecs(:, 309) * jflux / energy / hplanck
    kall(309) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 
    
    ! O2 -> O2+ + E
    f(:) = photo_xsecs(:, 310) * jflux / energy / hplanck
    kall(310) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0

    !! PREPROCESS_END

  end subroutine compute_photorates

end module prizmo_rates_photo
