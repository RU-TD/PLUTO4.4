module prizmo_rates_not_interpolable
contains

  subroutine evaluate_not_interpolable(n, Tgas, Tdust, jflux)
    use prizmo_commons
    use prizmo_utils
    use prizmo_self_shielding
    use prizmo_rates_photo
    use prizmo_H2_dust
    implicit none
    real*8,intent(in)::n(nmols), Tgas, Tdust, jflux(nphoto)
    real*8::ntot, mu, inv_Tdust, invsqrt32, prev, sqrt_Tgas
    real*8:: rho_dust, stick, pre_scale, pre_freeze, ndns, pre_2body

    ! get total number density
    ntot = get_ntot(n(:))

    ! mean molecular weight
    mu = get_mu(n(:))

    ! frequently used quantities
    inv_Tdust = 1d0 / Tdust
    invsqrt32 = (Tgas / 3e2)**(-0.5)
    sqrt_Tgas = sqrt(Tgas)
    prev = sqrt(8e0 * kboltzmann / pi)

    ! dust mass density, g/cm3
    rho_dust = ntot * mu * dust_d2g * proton_mass

    ! sticking coefficient
    stick = 1d0 / (1d0 + 4d-2 * sqrt(Tgas + Tdust) + 2d-3 * Tgas + 8d-6 * Tgas**2)

    ! dust content scaling factor
    pre_scale = rho_dust * dust_fact3 / dust_fact4 / pi43 / dust_rho_bulk

    ! freezeout prefactor, will be multiplied by 1/sqrt(mass_species)
    pre_freeze = stick * prev * sqrt_Tgas * pre_scale

    ! site density factor
    ndns = 4d0 * pi / dust_app2 * pre_scale

    ! 2body prefactor
    pre_2body = 1d0 / ndns

    !!BEGIN_NOT_INTERPOLABLE_RATES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! CH3+ + E -> C + H2 + H
    kall(288) = 1.05e-07*invsqrt32

    ! H + H -> H2
    kall(289) = get_H2_dust(Tgas, Tdust, pre_scale, n(idx_H))

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_NOT_INTERPOLABLE_RATES

  end subroutine evaluate_not_interpolable

end module prizmo_rates_not_interpolable
