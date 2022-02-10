module lkrm_rates
contains

  !*******************
  !compute rates and store into commons kall(:)
  subroutine compute_rates(n, Tgas, Tdust, variable_crflux, &
       dust_amin, dust_amax, dust_d2g, dust_pexp)
    use lkrm_commons, only: nmols, kall, mass, pi43
    implicit none
    real*8,intent(in)::n(nmols), Tgas, Tdust
    real*8,intent(in)::variable_crflux
    real*8,intent(in)::dust_amin, dust_amax, dust_d2g, dust_pexp
    real*8::inv_Tgas, sqrt_Tgas, pre_freeze, stick, p3, p4, rho_dust
    real*8::pre_scale, pre_2body, ntot, mu
    real*8,parameter::debye_nu=1d12 !Debye frequency, 1/s
    real*8,parameter::rho0=3d0 !bulk density, g/cm3
    real*8,parameter::app2=(3e-8)**2 !sites separation area, cm2

    ! total number density 1/cm3
    ntot = sum(n(:))

    ! mean molecular weight, g
    mu = sum(n(:)*mass(:)) / ntot

    !frequently used quantities
    inv_Tgas = 1d0/Tgas
    sqrt_Tgas = sqrt(Tgas)

    !dust mass density, g/cm3 (note: mu in g)
    rho_dust = ntot*mu*dust_d2g

    p3 = dust_pexp + 3d0
    p4 = dust_pexp + 4d0

    !sticking coefficient
    stick = 1d0/(1d0 + 4d-2*sqrt(Tgas+Tdust) + 2d-3*Tgas + 8d-6*Tgas**2)

    !dust content scaling factor
    pre_scale = rho_dust / rho0 * p4 / p3 * (dust_amax**p3-dust_amin**p3) &
         / (dust_amax**p4-dust_amin**p4)

    !freezeout prefactor
    pre_freeze = stick * sqrt_Tgas / pi43 * pre_scale

    !2body prefactor
    pre_2body = 3d0 / app2 * pre_scale

    !default rate
    kall(:) = 0d0

    !!BEGIN_RATES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2018-03-08 13:34:08
    ! CHANGESET: 3f0f625
    ! URL: https://GiovanniPicogna@bitbucket.org/tgrassi/mocassin_xray_chem.git
    ! BY: picogna@picogna-laptop

    ! H2 -> H + H
    kall(1) = 1e-12

    ! H + H -> H2
    kall(2) = 5e-13*Tgas

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_RATES

  end subroutine compute_rates
end module lkrm_rates

