module lkrm_main

  !!BEGIN_USER_COMMONS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2018-03-08 13:34:08
    ! CHANGESET: 3f0f625
    ! URL: https://GiovanniPicogna@bitbucket.org/tgrassi/mocassin_xray_chem.git
    ! BY: picogna@picogna-laptop

    ! number of species
    integer,parameter::lkrm_nmols=2

    ! number of reactions
    integer,parameter::lkrm_nrea=2

    ! number of energy bins
    integer,parameter::lkrm_nphoto=10000

    ! species indexes
    integer,parameter::lkrm_idx_H2=1
    integer,parameter::lkrm_idx_H=2

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_USER_COMMONS

contains
  !************************
  !evolve chemistry for a time-step dt (s)
  ! n(:) are species number densities
  subroutine lkrm(n, Tgas, jflux, dt)
    use lkrm_commons, only: nmols, nphoto
    use lkrm_ode, only: dochem
    use lkrm_rates, only: compute_rates
    implicit none
    real*8,intent(inout)::n(nmols)
    real*8,intent(in)::dt, Tgas, jflux(nphoto)

    real*8::variable_crflux, dust_amin, dust_amax, dust_d2g
    real*8::dust_pexp, Tdust

    ! cosmic rays ionization rate, 1/s
    variable_crflux = 5d-17

    ! dust temperature, K
    Tdust = Tgas

    ! dust variables
    dust_amin = 1e-7 ! min size, cm
    dust_amax = 1e-5 ! max size, cm
    dust_d2g = 1d-2 ! dust to gas mass ratio
    dust_pexp = -3.5 ! dust power-law exponent

    ! compute non-photochemmical rates
    call compute_rates(n(:), Tgas, Tdust, variable_crflux, &
         dust_amin, dust_amax, dust_d2g, dust_pexp)

    call dochem(n(:), Tgas, jflux(:), dt)

  end subroutine lkrm

  ! *******************
  ! initialize LKRM
  subroutine lkrm_init()
    use lkrm_xsecs, only: load_xsecs_all
    implicit none

    call load_xsecs_all()

  end subroutine lkrm_init

  ! *******************
  ! initialize LKRM, C interface
  subroutine lkrm_init_c() bind(C)
    use lkrm_xsecs, only: load_xsecs_all
    use iso_c_binding

    call load_xsecs_all()

  end subroutine lkrm_init_c

  ! *******************
  ! do chemistry, C interface
  subroutine lkrm_c(n, Tgas, dt) bind(C)
    use iso_c_binding
    use lkrm_commons, only: nmols, nphoto
    implicit none
    real(C_DOUBLE),intent(inout)::n(nmols)
    real(C_DOUBLE),intent(in)::dt, Tgas
    real*8::jflux(nphoto)

    ! set default radiation flux to zero
    jflux(:) = 0d0

    ! call chemistry
    call lkrm(n(:), Tgas, jflux(:), dt)

  end subroutine lkrm_c

  ! *********************
  ! call chemistry using mass fractions (x),
  ! mass density (rho, g/cm3), temperature (Tgas, K)
  ! and time-step (dt, s)
  subroutine lkrm_rho_c(x, rho, Tgas, dt) bind(C)
    use iso_c_binding
    use lkrm_commons, only: nmols, nphoto, mass
    implicit none
    real(C_DOUBLE),intent(inout)::x(nmols)
    real(C_DOUBLE),intent(in)::dt, Tgas, rho
    real*8::jflux(nphoto), n(nmols)

    ! set default radiation flux to zero
    jflux(:) = 0d0

    ! convert mass fractions to number density
    n(:) = x(:) * rho / mass(:)

    ! call chemistry
    call lkrm(n(:), Tgas, jflux(:), dt)

    ! convert number density to mass fractions
    x(:) = n(:) * mass(:) / rho

  end subroutine lkrm_rho_c

end module lkrm_main
